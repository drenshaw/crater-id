#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

// #include "camera.h"
// #include "quadrics.h"
// #include "conics.h"

class NavigationTest : public testing::Test {
protected:
  NavigationTest() {
  }
  ~NavigationTest() override {
    // delete quadric_default;
  }

public:
};