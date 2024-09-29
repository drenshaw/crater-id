#include "camera.h"
// #include "quadrics.h"
// #include "conics.h"
#include "math_utils.h"

#include "gtest/gtest.h"
#include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "opencv2/viz/types.hpp"

class CameraTest : public testing::Test {
  protected:
    CameraTest() {
      dx = 1000, dy = 1000, skew = 0, im_height = 2048, im_width = 2592;
      up = (im_width+1)/2;
      vp = (im_height+1)/2;
      image_size = cv::Size(im_width, im_height);
      cam = new Camera(dx, dy, up, vp, skew, image_size, quat, position);
      cv_cam = new cv::viz::Camera(dx, dy, up, vp, image_size);
    }
    ~CameraTest() override {
      delete cam;
      delete cv_cam;
    }
  public:
      double dx, dy, skew, im_height, im_width, up, vp;
      cv::Size2i image_size;
      Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();
      Eigen::Vector3d position{1,2,-3e4};

      Camera* cam;
      cv::viz::Camera* cv_cam;
      // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST_F(CameraTest, InitCamera) {
  // Eigen::Matrix4d extrinsic_h = cam->getHomogeneousExtrinsicMatrix();
  // Eigen::Matrix3d intrinsic   = cam->getIntrinsicMatrix();
  // Eigen::MatrixXd extrinsic   = cam->getExtrinsicMatrix();
  // Eigen::MatrixXd proj_mtx    = cam->getProjectionMatrix();
  // std::cout << "Extrinsic Matrix:\n" << extrinsic_h << std::endl;
  // std::cout << "Intrinsic Matrix:\n" << intrinsic << std::endl;
  // std::cout << "\nProjection Matrix:\n" << cam->getHomogeneousProjectionMatrix() << std::endl;

  // Quadric quad(-89, 0, 200, "south_pole");
  // Eigen::Matrix4d q_locus = quad.getLocus();
  // Eigen::Matrix3d c_locus = proj_mtx * q_locus * proj_mtx.transpose();
  // std::cout << "Projected locus:\n" << q_locus << std::endl << c_locus << std::endl;
  // Conic con(c_locus);
  // std::cout << "Projected conic: " << con << std::endl;

  }

TEST_F(CameraTest, CameraFOVX) {
  double fovx = 84.306181166432509;
  EXPECT_DOUBLE_EQ(cam->getFovXDeg(), fovx);
  double fovd = cam->getFovX();
  EXPECT_DOUBLE_EQ(rad2deg(fovd), fovx);
}

TEST_F(CameraTest, CameraFOVY) {
  double fovy = 64.043810747703901;
  EXPECT_DOUBLE_EQ(cam->getFovYDeg(), fovy);
  double fovd = cam->getFovY();
  EXPECT_DOUBLE_EQ(rad2deg(fovd), fovy);
}

TEST_F(CameraTest, PointMovingInRightDirection) {

  for(int i = 0; i < 10; i++) {
    // Eigen::Vector3d pt = {-i*5e3,0,0};
    Eigen::Vector3d pt = 100*latlon2bearing(0.0, i*36.);
    std::cout << "Bearing: " << pt.transpose() << std::endl;
  }
  for(int i = 0; i < 10; i++) {
    // Eigen::Vector3d pt = {-i*5e3,0,0};
    Eigen::Vector3d pt = 1e5*latlon2bearing(0.0, i*36.);
    Eigen::Vector2d proj_pt;
    cam->projectXYZtoImage(pt, proj_pt);
    std::cout << "Projected point: " << proj_pt.transpose() << " | " << pt.transpose()<< (cam->isInCameraFrame(pt) ? " :: yes":" :: no") << std::endl;

  }
}

TEST_F(CameraTest, Transformation) {
  // Eigen::Quaterniond rot1 = Eigen::Quaterniond::Identity();
  // std::cout << "Quat:" << rot1 << std::endl;
  Eigen::Transform<double, 3, Eigen::Projective>  transformation;
  // Eigen::Transform<double, 3, Eigen::AffineCompact3d>  iso;
  double rx = 0, ry = 0, rz = 10, rw = 10;
  Eigen::Quaterniond           rotation(Eigen::Quaterniond(rw, rx, ry, rz).normalized());
  Eigen::Translation<double,3> translation(Eigen::Vector3d(1,2,3));
  transformation = rotation * translation;
  Eigen::Vector3d vec = {0.1, 0.2, 10};
  Eigen::Vector3d x_vec = transformation.rotation()*vec + transformation.translation();

  // std::cout << "Rotation:\n" << transformation.rotation() << "\nTranslation:\n" << transformation.translation()<< std::endl;
  // std::cout << "Transform\n" << transformation.rotation()*vec + transformation.translation() << std::endl;

  Eigen::Matrix3d t_e2m;
  t_e2m <<  0, -1, 0, 
            1,  0, 0, 
            0,  0, 1;
  EXPECT_TRUE(rotation.toRotationMatrix().isApprox(t_e2m));
  EXPECT_TRUE(transformation.translation().isApprox(Eigen::Vector3d(-2, 1, 3)));
  EXPECT_TRUE(x_vec.isApprox(Eigen::Vector3d(-2.2, 1.1, 13.0)));
}

TEST_F(CameraTest, TransformState) {
  
  double rx = 0, ry = 0, rz = 0.5, rw = 0.87;
  Eigen::Quaterniond           rotation(Eigen::Quaterniond(rw, rx, ry, rz).normalized());
  Eigen::Quaterniond           attitude = Eigen::Quaterniond::Identity();
  Eigen::Translation<double,3> translation(Eigen::Vector3d(1e3,2,0));
  Eigen::Translation<double,3> position(Eigen::Vector3d(10,20,30));
  Eigen::Transform<double, 3, Eigen::Isometry>  transformation, state;
  transformation = rotation * translation;
  state = attitude * position;
  // std::cout << "Pre -transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;
  // state = transformation * state;
  // std::cout << "Post-transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;
  // state = transformation.inverse() * state;
  // std::cout << "Reg -transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;

  std::cout << *cam << std::endl;
  std::cout << "Applying transform: translation: ( " << transformation.translation().transpose()
            << " )\n\tRotation:\n" << transformation.rotation() << std::endl;
  cam->moveCamera(transformation);
  std::cout << *cam << std::endl;
  transformation = attitude * Eigen::Translation<double,3>(transformation.translation());
  cam->setAttitude(attitude);
  std::cout << *cam << std::endl;
  cam->moveCamera(rotation);
  std::cout << *cam << std::endl;
}