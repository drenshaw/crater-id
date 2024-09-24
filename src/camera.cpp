// #include <iostream>
// #include <vector>
#include <opencv2/viz/types.hpp>

// #include "math_utils.h"
#include "camera.h"

cv::viz::Camera camera(1000, 1000, 1296.5, 1024.5, cv::Size(2048,1920));

Camera::Camera(const double dx,
               const double dy,
               const double up,
               const double vp,
               const double skew) {
    makeCamera(dx, dy, up, vp, skew);
}
// Camera::Camera(const double dx,
//                const double dy,
//                const double up,
//                const double vp) {

//                }
Camera::Camera(const std::vector<double>& intrinsics) {
  makeCamera(intrinsics);
}
Camera::Camera(const Eigen::Matrix3d& intrinsic) {
  makeCamera(intrinsic);
}
Camera::Camera(const Eigen::Matrix3d& intrinsic,
               const Eigen::MatrixXd& extrinsic) {
  makeCamera(intrinsic, extrinsic);
}
Camera::Camera(const Eigen::Matrix3d& intrinsic, 
               const Eigen::Matrix3d& rotation,
               const Eigen::Vector3d& translation) {
  makeCamera(intrinsic, rotation, translation);
}
Camera::Camera(const Eigen::Matrix3d& intrinsic, 
               const Eigen::Quaterniond& rotation,
               const Eigen::Vector3d& translation) {
  makeCamera(intrinsic, rotation, translation, true);
}

void Camera::makeCamera(const double dx,
                        const double dy,
                        const double up,
                        const double vp,
                        const double skew) {
  intrinsic_matrix_ << dx, skew, up,
                        0,    dy, vp,
                        0,     0,  1;
  // extrinsic_matrix_ = Eigen::MatrixXd(3,4);
  attitude_ = Eigen::Matrix3d::Identity();
  // position_ = Eigen::Vector3d;
  position_ << 0, 0, 0;
}
void Camera::makeCamera(const std::vector<double>& intrinsics) {
  // auto [dx, dy, up, vp, skew] = intrinsics; // as a tuple
  double dx   = intrinsics[0];
  double dy   = intrinsics[1];
  double up   = intrinsics[2];
  double vp   = intrinsics[3];
  double skew = intrinsics[4];
  makeCamera(dx, dy, up, vp, skew);
}
void Camera::makeCamera(const Eigen::Matrix3d& intrinsic) {
  Eigen::MatrixXd extrinsic = Eigen::MatrixXd(3,4);
  // Eigen::Vector4d origin;
  // origin << 0, 0, 0, 1;
  extrinsic << Eigen::Matrix3d::Identity(), 0, 0, 0;
  makeCamera(intrinsic, extrinsic);
}
void Camera::makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::MatrixXd& extrinsic) {

}
void Camera::makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::Matrix3d& rotation,
                        const Eigen::Vector3d& translation) {

}
void Camera::makeCamera(const Eigen::Matrix3d& intrinsic, 
                        const Eigen::Quaterniond& rotation,
                        const Eigen::Vector3d& translation, bool yes) {

}
Eigen::Vector3d Camera::getPosition() {
  return position_;
}
Eigen::Matrix3d Camera::getAttitude() {
  return attitude_;
}
void Camera::get_attitude(Eigen::Quaterniond& quaternion) {
  quaternion = attitude_;
}
// Camera(const Eigen::MatrixXd& projection_matrix);