#include "gtest/gtest.h"
#include <eigen3/Eigen/Dense>
#include <opencv2/core/core.hpp> 
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include "opencv2/viz/types.hpp"

#include "math_utils.h"

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
  std::cerr << "Camera fov: \n" << rad2deg(camera.getFov()[0]) << " | " << rad2deg(camera.getFov()[1]) << std::endl;
  cv::Matx44d proj;
  camera.computeProjectionMatrix(proj);
  std::cerr << "Projection matrix:\n" << proj;
}