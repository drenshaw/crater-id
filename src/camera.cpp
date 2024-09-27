#include <iostream>
#include <iomanip>
#include <vector>
#include <opencv2/viz/types.hpp>

#include "math_utils.h"
#include "camera.h"
// TODO: Ensure that the "rotation" matrix is a transformation matrix, or transpose it

// This is the base constructor
Camera::Camera( const double dx,
                const double dy,
                const double up,
                const double vp,
                const double skew,
                const cv::Size2i& image_size,
                const Eigen::Quaterniond& attitude,
                const Eigen::Vector3d& position) :  
                  image_size_(image_size),
                  position_(position), 
                  attitude_(attitude.normalized()) {
  intrinsic_matrix_ << dx, skew, up,
                        0,   dy, vp,
                        0,    0,  1;
  Eigen::Translation3d translation(position);                        
  this->state_ = attitude * translation;                     
}

Camera::Camera() : Camera(1000, 1000, 
                          1296.5, 1024.5, 
                          0.0, 
                          cv::Size2i(2592, 2048)) {

}

// These are delegating constructors
Camera::Camera( const Eigen::Matrix3d& intrinsic_mtx, 
                const cv::Size2i& image_size,
                const Eigen::Quaterniond& attitude,
                const Eigen::Vector3d& position) : 
                  image_size_(image_size),
                  intrinsic_matrix_(intrinsic_mtx),
                  position_(position), 
                  attitude_(attitude.normalized()) {
  // Note the Eigen quaternion has the scalar value last (i.e., v0, v1, v2, s) 
  // when printing, BUUUUUUT when inputting the values, scalar part is first
  // However, you can access the elements using quat.vec() and quat.w()
  // https://tools.glowbuzzer.com/rotationconverter for visualizations
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

Eigen::Vector3d Camera::getPosition() const {
  return this->position_;
}

Eigen::Vector3d Camera::getLocation() const {
  return this->position_;
}

Eigen::Quaterniond Camera::getAttitude() const {
  return this->attitude_;
}

void Camera::getAttitude(Eigen::Quaterniond& quat) const {
  quat = this->attitude_;
}

Eigen::Quaterniond Camera::getAttitudeQuaternion() const {
  return this->getAttitude();
}

Eigen::Matrix3d Camera::getAttitudeMatrix() const {
  return this->attitude_.toRotationMatrix();
}

Eigen::Vector3d Camera::getAttitudeEuler(Eigen::Index a0, Eigen::Index a1, Eigen::Index a2) const {
  return this->getAttitudeMatrix().eulerAngles(a0, a1, a2);
}

Eigen::Matrix3d Camera::getIntrinsicMatrix() const {
  return this->intrinsic_matrix_;
}

Eigen::Matrix3d Camera::getInverseIntrinsicMatrix() const {
  Eigen::Matrix3d inv_intrinsic;
  double dx, dy, up, vp, skew;
  std::array<double, CAMERA_INTRINSIC_PARAM> params = this->getIntrinsicParams();
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

Eigen::Matrix4d Camera::getHomogeneousExtrinsicMatrix() const {
  Eigen::Matrix4d extrinsic_matrix;
  // extrinsic_matrix.topLeftCorner() == this->getAttitudeMatrix(), this->getPosition();
  extrinsic_matrix << this->getAttitudeMatrix(), this->getPosition(), Eigen::Vector4d::UnitW().transpose();
  return extrinsic_matrix;
}

Eigen::MatrixXd Camera::getExtrinsicMatrix() const {
  Eigen::Matrix4d extrinsic_matrix = this->getHomogeneousExtrinsicMatrix();
  return extrinsic_matrix.topRows(3);
}

Eigen::Matrix4d Camera::getHomogeneousProjectionMatrix() const {
  Eigen::MatrixXd nonhomogeneous_proj_mtx = this->getProjectionMatrix();
  Eigen::Vector4d unitv = Eigen::Vector4d::UnitW().transpose();
  Eigen::Matrix4d ext_mtx;
  ext_mtx << nonhomogeneous_proj_mtx, unitv.transpose();
  return ext_mtx;
}

Eigen::MatrixXd Camera::getProjectionMatrix() const {
  Eigen::Matrix3d intr = this->getIntrinsicMatrix();
  Eigen::MatrixXd extr = this->getExtrinsicMatrix();
  return intr * extr;
}

std::array<double, CAMERA_INTRINSIC_PARAM> Camera::getIntrinsicParams() const {
  std::array<double, CAMERA_INTRINSIC_PARAM> params;
  this->getIntrinsicParams(params);
  return params;
}

void Camera::getIntrinsicParams(std::array<double, CAMERA_INTRINSIC_PARAM>& params) const {
  Eigen::Matrix3d cam_intrinsic_mtx = this->getIntrinsicMatrix();
  getCameraIntrinsicParams(cam_intrinsic_mtx, params);
}

Eigen::Vector3d Camera::getPointWrtCameraFrame(const Eigen::Vector3d& pt) const {
  Eigen::Vector3d cam_position = this->getPosition();
  Eigen::Vector3d pt_cam = this->getAttitudeMatrix() * pt;
  return pt_cam - cam_position;
}

double Camera::getImageWidth() const {
  return this->image_size_.width;
}

double Camera::getImageHeight() const {
  return this->image_size_.height;
}

double Camera::getFovX() const {
  return getCameraFovX(this->intrinsic_matrix_, this->image_size_);
}

double Camera::getFovXDeg() const {
  return getCameraFovXDeg(this->intrinsic_matrix_, this->image_size_);
}

double Camera::getFovY() const {
  return getCameraFovY(this->intrinsic_matrix_, this->image_size_);
}

double Camera::getFovYDeg() const {
  return getCameraFovYDeg(this->intrinsic_matrix_, this->image_size_);
}

bool Camera::isInFrontOfCamera(const Eigen::Vector3d& pt) const {
  Eigen::Vector3d pt_wrt_cam_cam = this->getPointWrtCameraFrame(pt);
  return pt_wrt_cam_cam(2)>0;
}
    
bool Camera::isInCameraFrame(const Eigen::Vector3d& pt) const {
  Eigen::Vector2d pt_pxl;
  return isInCameraFrame(pt, pt_pxl);
}

bool Camera::pointUVinImage(const Eigen::Vector2d& pt_uv) const {
  return  pt_uv(0) >= 0 && pt_uv(0) <= this->image_size_.width &&
          pt_uv(1) >= 0 && pt_uv(1) <= this->image_size_.height;
}
    
bool Camera::isInCameraFrame(const Eigen::Vector3d& pt_xyz, Eigen::Vector2d& pt_pxl) const  {
  if(!isInFrontOfCamera(pt_xyz)) {
    return false;
  }
  projectXYZtoImage(pt_xyz, pt_pxl);
  return this->pointUVinImage(pt_pxl);
}

void Camera::projectXYZtoImage(const Eigen::Vector3d& pt, Eigen::Vector2d& pt_pxl) const {
  Eigen::Vector4d point = pt.homogeneous();
  Eigen::MatrixXd proj_mtx = this->getProjectionMatrix();
  Eigen::Vector3d proj_pt_h = proj_mtx * point;
  proj_pt_h /= proj_pt_h(2);
  pt_pxl = proj_pt_h.topRows(2);
}

Eigen::Vector2d Camera::projectXYZtoImage(const Eigen::Vector3d& pt) const {
  Eigen::Vector2d pt_pxl;
  this->projectXYZtoImage(pt, pt_pxl);
  return pt_pxl;
}


bool Camera::setAttitude(const Eigen::Quaterniond& orientation) {


  return true;
}

bool Camera::setAttitude(const Eigen::Matrix3d& orientation) {


  return true;
}

bool Camera::setAttitude(const Eigen::Vector3d& orientation) {


  return true;
}


bool Camera::setIntrinsicMatrix(const Eigen::Matrix3d& intrinsic) {


  return true;
}

bool Camera::setExtrinsicMatrix(const Eigen::Matrix3d& extrinsic) {


  return true;
}

void Camera::moveCamera(const Eigen::Transform<double, 3, Eigen::Isometry>& transform) {
  this->state_ = transform * this->state_;
}

// //Eigen::Transform<double, 3, Eigen::Isometry>
// Eigen::Matrix3d Camera::getStatePosition() const {
//   return this->state_.translation();
// }



std::ostream& operator<<(std::ostream& os, const Camera& cam) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  return os 
    << std::fixed << std::setw(8) << std::setprecision(1)
    << "Camera -> (" << cam.getLocation().transpose() << ") "
    << "\n\tAttitude: (" << cam.getAttitude() << ")"
    << "\n\tFov: " << cam.getFovXDeg() << ", " << cam.getFovYDeg() << std::endl
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

Eigen::Matrix3d pointCameraInDirection(const Eigen::Vector3d& camera_position, const Eigen::Vector3d& desired_location) {
  Eigen::Matrix3d T_m2c = getENUFrame(desired_location - camera_position);
  Eigen::Matrix3d z_rot;
  eulerToDCM(0., 0., M_PI, z_rot);
  return z_rot * T_m2c.transpose();
}
