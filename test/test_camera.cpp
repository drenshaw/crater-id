#include "camera.h"
#include "quadrics.h"
#include "conics.h"
// #include "math_utils.h"

#include "gtest/gtest.h"
#include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "opencv2/viz/types.hpp"

class CameraTest : public testing::Test {
  protected:
    CameraTest() {
      // std::array<double, 5> arr = {10.0, 7.0, 300.0, 50.0, 0.0};
      // Conic conic_arr(arr);
      // c_arr.SetGeometricParameters(arr);
    }
  public:
      // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST(CameraTest, InitCamera) {
  cv::viz::Camera camera(1000, 1000, 1296.5, 1024.5, cv::Size(2048,1920));
  // std::cerr << "Camera fov: \n" << rad2deg(camera.getFov()[0]) << " | " << rad2deg(camera.getFov()[1]) << std::endl;
  cv::Matx44d proj;
  camera.computeProjectionMatrix(proj);
  // std::cerr << "Projection matrix:\n" << proj;
  Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();

  Eigen::Vector3d position{1,2,-3e4};
  cv::Size2i image_size(2048, 1920);
  Camera cam(1000, 1000, 1296.5, 1024.5, 0.0, image_size, quat, position);
  Eigen::Matrix4d extrinsic_h = cam.getHomogeneousExtrinsicMatrix();
  Eigen::Matrix3d intrinsic = cam.getIntrinsicMatrix();
  Eigen::MatrixXd extrinsic = cam.getExtrinsicMatrix();
  Eigen::MatrixXd proj_mtx = cam.getProjectionMatrix();
  std::cout << "Extrinsic Matrix:\n" << extrinsic_h << std::endl;
  std::cout << "Intrinsic Matrix:\n" << intrinsic << std::endl;
  std::cout << "Intrinsic Matrix:\n" << intrinsic << std::endl;
  std::cout << "\nExtrinsic (nonhomogeneous) Matrix:\n" << extrinsic << std::endl;
  std::cout << "\nExpected Projection Matrix:\n" << intrinsic*extrinsic << std::endl;
  std::cout << "\nProjection Matrix:\n" << cam.getHomogeneousProjectionMatrix() << std::endl;

  Quadric quad(-89, 0, 200, "south_pole");
  Eigen::Matrix4d q_locus = quad.getLocus();
  Eigen::Matrix3d c_locus = proj_mtx * q_locus * proj_mtx.transpose();
  std::cout << "Projected locus:\n" << q_locus << std::endl << c_locus << std::endl;
  Conic con(c_locus);
  std::cout << "Projected conic: " << con << std::endl;

  for(int i = 0; i < 10; i++) {
    Eigen::Vector3d pt = {i*5e3,0,0};
    Eigen::Vector2d proj_pt;
    cam.projectXYZtoImage(pt, proj_pt);
    std::cout << "Projected point: " << proj_pt.transpose() << " | " << pt.transpose() << std::endl;
    // std::cout << "Point is in front of camera? " << (inFrame ? "yes":"no") << std::endl;

  }
}

// TODO LIST: flesh out the camera class and constructors, then test, then work on projections