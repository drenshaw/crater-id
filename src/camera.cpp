#include <iostream>
#include <iomanip>
#include <vector>
#include <opencv2/viz/types.hpp>

#include "math_utils.h"
#include "camera.h"
// TODO: Ensure that the "rotation" matrix is a transformation matrix, or transpose it

// This is the base constructor
// We will store the transformation as an active one since Eigen seems to prefer it
Camera::Camera( const double dx,
                const double dy,
                const double up,
                const double vp,
                const double skew,
                const cv::Size2i& image_size,
                const Eigen::Isometry3d& transform) :  
                  dx_(dx), dy_(dy), skew_(skew), up_(up), vp_(vp),
                  image_size_(image_size) {
  // TODO: does this work? Or do we need to invert the whole transform?                        
  state_.linear() = transform.rotation().inverse();
}
Camera::Camera( const double dx,
                const double dy,
                const double up,
                const double vp,
                const double skew,
                const cv::Size2i& image_size,
                const Eigen::Quaterniond& attitude,
                const Eigen::Vector3d& position) :  
                Camera( dx, dy, up, vp, skew, 
                        image_size, 
                        attitude * Eigen::Translation3d(position)) {}

Camera::Camera() : Camera(1000, 1000, 
                          1296.5, 1024.5, 
                          0.0, 
                          cv::Size2i(2592, 2048)) {

}

Camera::Camera( const Eigen::Matrix3d& intrinsic, 
                const cv::Size2i& image_size,
                const Eigen::Isometry3d& transform) {

}

// These are delegating constructors
Camera::Camera( const Eigen::Matrix3d& intrinsic_mtx, 
                const cv::Size2i& image_size,
                const Eigen::Quaterniond& attitude,
                const Eigen::Vector3d& position) : 
                  image_size_(image_size),
                  intrinsic_matrix_(intrinsic_mtx),
                  state_(attitude * Eigen::Translation3d(position)) {
  // Note the Eigen quaternion has the scalar value last (i.e., v0, v1, v2, s) 
  // when printing, BUUUUUUT when inputting the values, scalar part is first
  // However, you can access the elements using quat.vec() and quat.w()
  // https://tools.glowbuzzer.com/rotationconverter for visualizations
  std::array<double, CAMERA_INTRINSIC_PARAM> params;
  getCameraIntrinsicParams(intrinsic_mtx, params);
  this->dx_  = params.at(0);
  this->dy_  = params.at(1);
  this->up_  = params.at(2);
  this->vp_  = params.at(3);
  this->skew_= params.at(4);
}
Camera::Camera( const std::vector<double>& intrinsics, 
                const cv::Size2i& image_size,
                const Eigen::Quaterniond& attitude,
                const Eigen::Vector3d& position) : 
  Camera( intrinsics.at(0),
          intrinsics.at(1),
          intrinsics.at(2),
          intrinsics.at(3),
          intrinsics.at(4),
          image_size,
          attitude, position) {
}
Camera::Camera(const Eigen::Matrix3d& intrinsic,
                const cv::Size2i& image_size,
               const Eigen::Matrix4d& extrinsic) :
  Camera( intrinsic, 
          image_size,
          extrinsic.topLeftCorner(3,3), 
          extrinsic.topRightCorner(3,1)) {
}
Camera::Camera( const Eigen::Matrix3d& intrinsic, 
                const cv::Size2i& image_size,
                const Eigen::Matrix3d& attitude,
                const Eigen::Vector3d& translation) :
                Camera( intrinsic,
                        image_size,
                        Eigen::Quaterniond(attitude),
                        translation) {
}

Eigen::Isometry3d Camera::getTransformation() const {
  return this->state_;
}

Eigen::Vector3d Camera::getPosition() const {
  return this->getTransformation().translation();
}

Eigen::Vector3d Camera::getLocation() const {
  return this->getPosition();
}

Eigen::Quaterniond Camera::getAttitude() const {
  Eigen::Quaterniond quat;
  this->getAttitude(quat);
  return quat.inverse();
}

void Camera::getAttitude(Eigen::Quaterniond& quat) const {
  quat = Eigen::Quaterniond(this->state_.rotation());
}

Eigen::Quaterniond Camera::getAttitudeQuaternion() const {
  return this->getAttitude();
}

Eigen::Matrix3d Camera::getAttitudeMatrix() const {
  return this->getAttitude().toRotationMatrix();
}

Eigen::Vector3d Camera::getAttitudeEuler( Eigen::Index a0, 
                                          Eigen::Index a1, 
                                          Eigen::Index a2) const {
  return this->getAttitudeMatrix().eulerAngles(a0, a1, a2);
}

Eigen::Matrix3d Camera::getIntrinsicMatrix() const {
  Eigen::Matrix3d intrinsic_matrix;
  intrinsic_matrix << dx_, skew_, up_,
                        0,   dy_, vp_,
                        0,     0,   1;
  return intrinsic_matrix;
}

Eigen::Matrix3d Camera::getInverseIntrinsicMatrix() const {
  Eigen::Matrix3d inv_intrinsic;
  inv_intrinsic <<  1/dx_, -skew_/(dx_*dy_), (skew_*vp_ - dy_*up_)/(dx_*dy_),
                    0,                1/dy_,                        -vp_/dy_,
                    0,                    0,                               1;
  return inv_intrinsic;
}

Eigen::Isometry3d Camera::getHomogeneousExtrinsicMatrix() const {
  return this->state_.inverse();
}

Eigen::MatrixXd Camera::getExtrinsicMatrix() const {
  Eigen::Quaterniond att = this->getAttitude();
  Eigen::Vector3d pos = this->getPosition();
  Eigen::MatrixXd extrinsic(3,4);
  extrinsic.topLeftCorner(3,3) = Eigen::Matrix3d::Identity();
  extrinsic.topRightCorner(3,1) = -pos;
  extrinsic = att.toRotationMatrix() * extrinsic;
  // Eigen::Affine3d extrinsic_matrix = this->getHomogeneousExtrinsicMatrix();
  return extrinsic;
}

Eigen::Affine3d Camera::getHomogeneousProjectionMatrix() const {
  Eigen::Affine3d pmtx = Eigen::Affine3d::Identity();
  pmtx.matrix().topRows(3) = this->getProjectionMatrix();
  return pmtx;
}

Eigen::MatrixXd Camera::getProjectionMatrix() const {
  Eigen::MatrixXd proj_mtx(3,4);
  proj_mtx = this->getIntrinsicMatrix() * this->getExtrinsicMatrix();
  return proj_mtx;
}

std::array<double, CAMERA_INTRINSIC_PARAM> Camera::getIntrinsicParams() const {
  std::array<double, CAMERA_INTRINSIC_PARAM> params;
  this->getIntrinsicParams(params);
  return params;
}

void Camera::getIntrinsicParams(std::array<double, CAMERA_INTRINSIC_PARAM>& params) const {
  params = {dx_, dy_, up_, vp_, skew_};
}

cv::Mat Camera::getBlankCameraImage() const {
  // cv::Mat uses [rows, cols]
  cv::Mat image(this->getImageHeight(), this->getImageWidth(), 
                CV_8UC3, 
                cv::Scalar(50, 50, 50));
  return image;
}

void Camera::resetImage(cv::Mat& image) const {
  image = cv::Scalar(50, 50, 50);
}

double Camera::getImageWidth() const {
  return this->image_size_.width;
}

double Camera::getImageHeight() const {
  return this->image_size_.height;
}

void Camera::getImageMidpoint(Eigen::Vector2d& img_midpoint) const {
  img_midpoint = {this->up_, this->vp_};
}

Eigen::Vector2d Camera::getImageMidpoint() const {
  Eigen::Vector2d img_midpoint;
  this->getImageMidpoint(img_midpoint);
  return img_midpoint;
}

cv::Size2i Camera::getImageSize() const {
  return this->image_size_;
}

double Camera::getFovX() const {
  return getCameraFovX(this->getIntrinsicMatrix(), this->image_size_);
}

double Camera::getFovXDeg() const {
  return getCameraFovXDeg(this->getIntrinsicMatrix(), this->image_size_);
}

double Camera::getFovY() const {
  return getCameraFovY(this->getIntrinsicMatrix(), this->image_size_);
}

double Camera::getFovYDeg() const {
  return getCameraFovYDeg(this->getIntrinsicMatrix(), this->image_size_);
}

void Camera::world2Camera(const Eigen::Vector3d& pt, Eigen::Vector3d& pt_cam) const {
  pt_cam = this->getHomogeneousExtrinsicMatrix() * pt;
  // std::cout << "WORLD2CAMERA val: " << pt_cam.transpose() << std::endl;
  if(std::abs(pt_cam.z())<1e-9) {
    throw std::runtime_error("Point lies in the z=0 camera plane; projection not possible.");
  }
}

void Camera::world2Pixel(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const {
  // Eigen::Vector4d point = pt.homogeneous();
  Eigen::Matrix3d intrinsic = this->getIntrinsicMatrix();
  Eigen::Vector3d pt_homogeneous = intrinsic * this->world2Camera(pt);
  // std::cout << "WORLD2PIXEL Pt (homogeneous): " << pt_homogeneous.transpose() << std::endl;
  pt_pxl = pt_homogeneous.hnormalized();
}

Eigen::Vector3d Camera::world2Camera(const Eigen::Vector3d& pt) const {
  Eigen::Vector3d pt_cam;
  this->world2Camera(pt, pt_cam);
  return pt_cam;
}

void Camera::world2Camera(const std::vector<Eigen::Vector3d>& pts, std::vector<Eigen::Vector3d>& pts_cam) const {
  pts_cam.clear();
  std::vector<Eigen::Vector3d>::const_iterator it;
  for(it = pts.begin(); it != pts.end(); it++) {
    Eigen::Vector3d pt_cam;
    this->world2Camera(*it, pt_cam);
    pts_cam.push_back(pt_cam);
  }
}

Eigen::Vector2d Camera::world2Pixel(const Eigen::Vector3d& pt) const {
  Eigen::Vector2d pt_pxl;
  this->world2Pixel(pt, pt_pxl);
  return pt_pxl;
}

void Camera::world2Pixel(const std::vector<Eigen::Vector3d>& pts, std::vector<Eigen::Vector2d>& pts_pxl) const {
  pts_pxl.clear();
  std::vector<Eigen::Vector3d>::const_iterator it;
  for(it = pts.begin(); it != pts.end(); it++) {
    Eigen::Vector2d pt_pxl;
    this->world2Pixel(*it, pt_pxl);
    pts_pxl.push_back(pt_pxl);
  }
}

bool Camera::isInFrontOfCamera(const Eigen::Vector3d& pt) const {
  Eigen::Vector3d pt_wrt_cam_cam = this->world2Camera(pt);
  return pt_wrt_cam_cam.z()>0;
}
    
bool Camera::isInCameraFrame(const Eigen::Vector3d& pt_xyz, Eigen::Vector2d& pt_pxl) const  {
  if(!isInFrontOfCamera(pt_xyz)) {
    return false;
  }
  world2Pixel(pt_xyz, pt_pxl);
  return this->isInPixelArray(pt_pxl);
}
    
bool Camera::isInCameraFrame(const Eigen::Vector2d& pt_pxl) const  {
  return this->isInPixelArray(pt_pxl);
}
    
bool Camera::isInCameraFrame(const Eigen::Vector3d& pt) const {
  Eigen::Vector2d pt_pxl;
  return isInCameraFrame(pt, pt_pxl);
}

bool Camera::isInPixelArray(const Eigen::Vector2d& pt_uv) const {
  return isInImage(pt_uv, this->image_size_);
}

void Camera::setIntrinsicMatrix(const Eigen::Matrix3d& intrinsic) {
  throw std::runtime_error("Function is not yet implemented.");
}

void Camera::setExtrinsicMatrix(const Eigen::Matrix3d& extrinsic) {
  throw std::runtime_error("Function is not yet implemented.");
}

void Camera::setAttitude(const Eigen::Quaterniond& orientation) {
  this->state_.linear() = orientation.normalized().inverse().toRotationMatrix();
  // std::cout << __func__ << ": " << orientation << std::endl;
}

void Camera::setAttitude(const Eigen::Matrix3d& orientation) {
  this->setAttitude(Eigen::Quaterniond(orientation));
  // std::cout << __func__ << ":\n" << orientation << std::endl;
}

void Camera::setAttitude(const Eigen::AngleAxisd& orientation) {
  this->setAttitude(Eigen::Quaterniond(orientation));
}

void Camera::setLocation(const Eigen::Vector3d& location) {
  this->state_.translation() = location;
  // std::cout << __func__ << ": " << location.transpose() << std::endl;
}

void Camera::setPosition(const Eigen::Vector3d& location) {
  this->setLocation(location);
}

void Camera::moveCamera(const Eigen::Isometry3d& transform) {
  this->state_.rotate(transform.rotation()).pretranslate(transform.translation());
}

void Camera::moveX(const double x_offset) {
  Eigen::Vector3d movement;
  movement << x_offset, 0, 0;
  this->move(movement);
}

void Camera::moveY(const double y_offset) {
  Eigen::Vector3d movement = {0, y_offset, 0};
  this->move(movement);
}

void Camera::moveZ(const double z_offset) {
  Eigen::Vector3d movement = {0, 0, z_offset};
  this->move(movement);
}


void Camera::rotate(const Eigen::Quaterniond& rotation) {
  this->state_.rotate(rotation);
}

void Camera::rotate(const Eigen::Matrix3d& rotation) {
  this->state_.rotate(rotation);
}

void Camera::rotate(const Eigen::AngleAxisd& rotation) {
  this->rotate(rotation.toRotationMatrix());
}

// TODO: need to determine/express if the translation is in the inertial or body frame
void Camera::move(const Eigen::Vector3d& translation) {
  this->state_.pretranslate(translation);
}
void Camera::move(const Eigen::Translation3d& translation) {
  this->state_.pretranslate(translation.translation());
}

void Camera::rotateRelative(const Eigen::Quaterniond& rotation) {
  this->state_.prerotate(rotation);
}

// this is intended to make it possible to move the camera in the local frame (e.g., forward)
void Camera::moveRelative(const Eigen::Vector3d& translation) {
  this->state_.translate(translation);
}
void Camera::moveRelative(const Eigen::Translation3d& translation) {
  this->state_.translate(translation.translation());
}

void Camera::pointTo(const Eigen::Vector3d& point, const Eigen::Vector3d& up_axis) {
  if((point - this->getPosition()).norm() < 1e-3) {
    // std::cerr << "Cannot point the camera in the direction specified; too close.\n";
    throw std::runtime_error("Cannot point the camera in the direction specified; too close.");
  }
  Eigen::Matrix3d xform = lookAt(this->getPosition(), point, up_axis);
  this->state_.linear() = xform.transpose();
}

void Camera::pointTo(const double lat, const double lon, const Eigen::Vector3d& up_axis) {
  Eigen::Vector3d point = latlonalt(lat, lon, 0);
  this->pointTo(point, up_axis);
}

void Camera::moveTo(const Eigen::Vector3d& point) {
  this->setPosition(point);
}

void Camera::moveTo(const Eigen::Translation3d& translation) {
  this->setPosition(translation.translation());
}

void Camera::resetCameraState() {
  this->state_.fromPositionOrientationScale(
    Eigen::Vector3d::Zero(), 
    Eigen::Quaterniond::Identity(), 
    Eigen::Vector3d::Ones());
}

Eigen::Matrix3d Camera::projectQuadric(const Eigen::Matrix4d& quadric_envelope) const {
  // TODO: account for whether the crater rim is occluded by the lunar surface
  // Returns the envelope of the projected conic
  Eigen::MatrixXd proj = this->getProjectionMatrix();
  return proj * quadric_envelope * proj.transpose();
}

Eigen::Matrix3d Camera::projectQuadricToLocus(const Eigen::Matrix4d& quadric_locus) const {
  // TODO: account for whether the crater rim is occluded by the lunar surface
  // Returns the envelope of the projected conic
  Eigen::MatrixXd proj = this->getProjectionMatrix();
  return adjugate(proj * adjugate(quadric_locus) * proj.transpose());
}

Eigen::Matrix3d Camera::getImagePlaneLocus(const Eigen::Matrix3d& image_locus) const {
  Eigen::Matrix3d Kinv = this->getInverseIntrinsicMatrix();
  Eigen::Matrix3d envelope = adjugate(image_locus);
  return adjugate(Kinv * envelope * Kinv.transpose());
}

Eigen::Matrix3d Camera::getMoonConic(const double radius) const {
  Eigen::Matrix4d sphere = makeSphere(double(R_MOON));
  return this->projectQuadricToLocus(sphere);
}


std::ostream& operator<<(std::ostream& os, const Camera& cam) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  return os 
    << std::fixed << std::setw(8) << std::setprecision(1)
    <<   "Camera ->  XYZ = (" << cam.getLocation().transpose() << ") km"
    << std::setw(8) << std::setprecision(5)
    << "\n           Att = (" << cam.getAttitude() << ")"
    << std::setw(8) << std::setprecision(2)
    << "\n           Fov = ( " << cam.getFovXDeg() << ", " << cam.getFovYDeg() << " ) deg"
    << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
    ;
}

std::ostream& operator<<(std::ostream& os, const Camera* cam) {
  return operator<<(os, *cam);
}

/***********************Utils**************************/


void getCameraIntrinsicParams(const Eigen::Matrix3d& cam_intrinsic_mtx,
                        std::array<double, CAMERA_INTRINSIC_PARAM>& params) {
  params = {
    cam_intrinsic_mtx(0, 0),
    cam_intrinsic_mtx(0, 1),
    cam_intrinsic_mtx(0, 2),
    cam_intrinsic_mtx(1, 1),
    cam_intrinsic_mtx(1, 2)
  };
}

Eigen::Matrix3d getCameraInverseIntrinsicMatrix(const Eigen::Matrix3d& camera_intrinsic_mtx) {
  Eigen::Matrix3d inv_intrinsic;
  double dx, dy, up, vp, skew;
  std::array<double, CAMERA_INTRINSIC_PARAM> params;
  getCameraIntrinsicParams(camera_intrinsic_mtx, params);
  dx   = params.at(0);
  skew = params.at(1);
  up   = params.at(2);
  dy   = params.at(3);
  vp   = params.at(4);
  inv_intrinsic <<  1/dx, -skew/(dx*dy), (skew*vp - dy*up)/(dx*dy),
                    0, 1/dy, -vp/dy,
                    0, 0, 1;
  return inv_intrinsic;
}
double getCameraFovX(const Eigen::Matrix3d& camera_intrinsic_mtx, const cv::Size2i& image_size) {
  double width = image_size.width;
  Eigen::Matrix3d inv_intrinsic = getCameraInverseIntrinsicMatrix(camera_intrinsic_mtx);
  Eigen::Vector3d topLeft = Eigen::Vector3d::UnitZ();
  Eigen::Vector3d topRight(width, 0, 1);
  Eigen::Vector3d p1 = (inv_intrinsic * topLeft ).normalized();
  Eigen::Vector3d p2 = (inv_intrinsic * topRight).normalized();
  return std::acos(p1.dot(p2));
}

double getCameraFovXDeg(const Eigen::Matrix3d &camera_intrinsic_mtx, const cv::Size2i &image_size) {
  return rad2deg(getCameraFovX(camera_intrinsic_mtx, image_size));
}

double getCameraFovY(const Eigen::Matrix3d& camera_intrinsic_mtx, const cv::Size2i& image_size) {
  double height = image_size.height;
  Eigen::Matrix3d inv_intrinsic = getCameraInverseIntrinsicMatrix(camera_intrinsic_mtx);
  Eigen::Vector3d topLeft = Eigen::Vector3d::UnitZ();
  Eigen::Vector3d bottomLeft(0, height, 1);
  Eigen::Vector3d p1 = (inv_intrinsic * topLeft   ).normalized();
  Eigen::Vector3d p2 = (inv_intrinsic * bottomLeft).normalized();
  return std::acos(p1.dot(p2));
}

double getCameraFovYDeg(const Eigen::Matrix3d &camera_intrinsic_mtx, const cv::Size2i &image_size) {
  return rad2deg(getCameraFovY(camera_intrinsic_mtx, image_size));
}

// This produces a passive transformation
// Read Zanetti's "Rotations, Transformations, Left Quaternions, Right Quaternions?"
Eigen::Matrix3d getAttitudeTransformBetweenPoints(const Eigen::Vector3d& camera_position, const Eigen::Vector3d& desired_location) {
  Eigen::Matrix3d T_m2c = getENUFrame(desired_location - camera_position);
  Eigen::Matrix3d z_rot;
  eulerToDCM(0., 0., M_PI, z_rot);
  return z_rot * T_m2c.transpose();
}

// This produces an active rotation
Eigen::Matrix3d lookAt( const Eigen::Vector3d& camera_position, 
                        const Eigen::Vector3d& desired_location, 
                        const Eigen::Vector3d& up_vector) {
  Eigen::Vector3d zaxis((desired_location - camera_position));
  assert(!zaxis.normalized().isApprox(up_vector.normalized()));
  Eigen::Vector3d xaxis = zaxis.cross(up_vector);
  Eigen::Vector3d yaxis = zaxis.cross(xaxis);

  Eigen::Matrix3d mtx = Eigen::Matrix3d::Identity();
  // mtx << xaxis.normalized(), yaxis.normalized(), zaxis.normalized();
  mtx(0,Eigen::seq(0,Eigen::last)) = xaxis.normalized();
  mtx(1,Eigen::seq(0,Eigen::last)) = yaxis.normalized();
  mtx(2,Eigen::seq(0,Eigen::last)) = zaxis.normalized();
  // mtx(0,3) = -xaxis.dot(camera_position);
  // mtx(1,3) = -yaxis.dot(camera_position);
  // mtx(2,3) = -zaxis.dot(camera_position);
  return mtx;
}

Eigen::Matrix3d lookAt( const Eigen::Vector3d& camera_position, 
                        const Eigen::Vector3d& desired_location) {
  Eigen::Vector3d zaxis((desired_location - camera_position).normalized());
  // Form the up vector based on switching the axes of the z-axis
  Eigen::Vector3d up_vector(zaxis(1), zaxis(2), zaxis(0));
  return lookAt(camera_position, desired_location, up_vector);
}

Eigen::Matrix3d lookAt( const Eigen::Isometry3d& transform, 
                        const Eigen::Vector3d& desired_location, 
                        const Eigen::Vector3d& up_vector) {
  return lookAt(transform.translation(), desired_location, up_vector).transpose();
}

bool isInImage(const Eigen::Vector2d& pt_uv, const cv::Size2i image_size) {
  cv::Point uv(pt_uv(1), pt_uv(0));
  return isInImage(uv, image_size);
}

bool isInImage(const cv::Point& pt_uv, const cv::Size2i image_size) {
  return  pt_uv.x >= 0 && pt_uv.x <= image_size.width &&
          pt_uv.y >= 0 && pt_uv.y <= image_size.height;
}

bool isInImage(const cv::Point& pt_uv, const cv::MatSize image_size) {
  // MatSize is [col, row]
  cv::Size2i img_size;
  img_size.width = image_size[0];
  img_size.height = image_size[1];
  return isInImage(pt_uv, img_size);
}
