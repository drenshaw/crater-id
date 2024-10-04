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

bool projectPointToPixel( const Eigen::Matrix3d& intrinsic, 
                          const Eigen::Isometry3d& transformation, 
                          const Eigen::Vector3d& point, Eigen::Vector2d& pixel) {
  // Transform the point from world coordinates to camera coordinates
  // std::cout << transformation.inverse().matrix() << std::endl;
  Eigen::Vector3d cameraPoint = transformation.inverse() * point;
  std::cout << "Camera frame point: " << cameraPoint.transpose() << std::endl;
  cv::Size2i image_size(2*intrinsic(0,2)-1, 2*intrinsic(1,2)-1);
  std::cout << "Image size: " << image_size << std::endl;

  // Check if the point is in front of the camera
  if (cameraPoint.z() <= 0) {
    std::cout << "The point is behind the camera." << std::endl;
    return false;
  }

  // Project the 3D point to 2D using the intrinsic matrix
  Eigen::Vector3d homogeneousPixel = intrinsic * cameraPoint;

  // Convert from homogeneous coordinates to 2D pixel coordinates
  pixel = homogeneousPixel.hnormalized();
  std::cout << "Pixel: " << pixel.transpose() << std::endl;
  return isInImage(pixel, image_size);
}

TEST_F(CameraTest, PointMovingInRightDirection) {
  Eigen::Isometry3d transformation = Eigen::Isometry3d::Identity();
  Eigen::AngleAxisd rot = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  // Eigen::AngleAxisd rot1 = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  Eigen::Vector3d translation(1, 0, 0);
  transformation.rotate(rot);
  transformation.pretranslate(translation);

  cam->moveCamera(Eigen::Quaterniond(rot));
  cam->moveCamera(translation);

  // Example 3D point
  Eigen::Vector3d point(1, 0, 5);
  Eigen::Vector2d pt_pxl;
  cam->world2Pixel(point, pt_pxl);
  Eigen::Vector2d expected_pixel(296.5, 1024.5);
  ASSERT_TRUE(pt_pxl.isApprox(expected_pixel));
  
  Eigen::Vector3d cam_pos = 1e4*latlon2bearing(0, 0);
  Eigen::Vector3d look_at_me = 1e4*Eigen::Vector3d::Zero();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->setPosition(cam_pos);
  cam->pointTo(look_at_me, up_vector);
  // Eigen::Matrix3d mtx = lookAt(cam_pos, look_at_me, up_vector);
  // std::cout << Eigen::Quaterniond(mtx) << std::endl;
  // cam->setAttitude(mtx);
  // std::cout << *cam << std::endl;
  for(int i = 0; i < 10; i++) {
    Eigen::Vector3d pt = {i*2.5e3,0,0};
    Eigen::Vector2d pixel;
    if(i == 4) {
      ASSERT_THROW(cam->world2Pixel(pt, pixel), std::runtime_error);
      continue;
    }
    cam->world2Pixel(pt, pixel);
    if(i < 4){
      EXPECT_TRUE(cam->isInCameraFrame(pt));
    }
    else {
      EXPECT_FALSE(cam->isInCameraFrame(pt));
    }
    // std::cout 
    //   << "Projected point: " << pixel.transpose() << " | " 
    //   << pt.transpose() << " | (cam) "
    //   << cam->world2Camera(pt).transpose()
    //   << (cam->isInCameraFrame(pt) ? " :: yes":" :: no") 
    //   << std::endl;

  }
}

TEST_F(CameraTest, Transformation) {
  Eigen::Transform<double, 3, Eigen::Projective>  transformation;
  // Eigen::Transform<double, 3, Eigen::AffineCompact3d>  iso;
  double rx = 0, ry = 0, rz = 10, rw = 10;
  Eigen::Quaterniond           rotation(Eigen::Quaterniond(rw, rx, ry, rz).normalized());
  Eigen::Translation3d translation(Eigen::Vector3d(1,2,3));
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
  
  // // double rx = 0, ry = 0, rz = 0.5, rw = 0.87;
  // double roll = 0, pitch = 0, yaw = 60;
  // Eigen::Quaterniond            euler(Eigen::Matrix3d(Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) * 
  //                                                             Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) * 
  //                                                             Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX())));
  // // Eigen::Quaterniond            rotation(Eigen::Quaterniond(rw, rx, ry, rz).normalized());
  // Eigen::Translation3d  translation(Eigen::Vector3d(50,2,0));
  // Eigen::Vector3d position(10,20,30);
  // Eigen::Isometry3d  transformation, state;
  // transformation = euler * translation;
  // state.setIdentity();
  // state.translate(position);

  // // TODO: Figure out the prerotate and pretranslate stuff; weird things happening
  // // TODO: Also ensure that transformation info is output, not active rotation
  // std::cout << "Pre -transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;
  // // state.rotate(euler);
  // // std::cout << "Post-transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;
  // // state.rotate(euler.inverse());
  // // std::cout << "Reg -transform:\n" << state.rotation() << "\nTrans:\n" << state.translation() << std::endl;

  // // std::cout << *cam << std::endl;
  // // std::cout << "Applying transform: translation: ( " << transformation.translation().transpose()
  // //           << " )\n\tRotation:\n" << Eigen::Quaterniond(transformation.rotation()) << std::endl;
  // std::cout << *cam << std::endl;
  // cam->moveCamera(transformation);
  // std::cout << *cam << std::endl;
  // // std::cout << *cam << std::endl;
  // // transformation = attitude * Eigen::Translation3d(transformation.translation());
  // // cam->setAttitude(attitude);
  // // std::cout << *cam << std::endl;
  // cam->moveCamera(euler);
  // std::cout << *cam << std::endl;
  // cam->moveCamera(translation);
  // std::cout << *cam << std::endl;
  // cam->moveCamera(translation);
}

TEST_F(CameraTest, EnsurePassiveXform) {
  // TODO: wait for Paul to send me unit tests for this one
}

void printTransformation(const Eigen::Isometry3d& xform, const std::string& str) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  Eigen::Quaterniond mtx = Eigen::Quaterniond(xform.rotation());
  std::cout << str << std::fixed 
            // << std::setw(8) << std::setprecision(1)
            // << xform.translation().transpose()
            << std::setw(8) << std::setprecision(5)
            // << std::endl 
            << mtx.x() << ", "
            << mtx.y() << ", "
            << mtx.z() << ", "
            << mtx.w()
            << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
            << std::endl;
}

TEST_F(CameraTest, MoveRelative) {
  // // state.fromPositionOrientationScale(Eigen::Vector3d::Zero(), Eigen::Quaterniond::Identity(), Eigen::Vector3d::Ones());
  // double roll = deg2rad(5), pitch = deg2rad(15), yaw = deg2rad(60);
  // Eigen::Quaterniond            euler(Eigen::Matrix3d(Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) * 
  //                                                             Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) * 
  //                                                             Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX())));
  // Eigen::Quaterniond            stata(Eigen::Matrix3d(Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitZ()) * 
  //                                                             Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitY()) * 
  //                                                             Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitX())));
  // Eigen::Isometry3d state = Eigen::Isometry3d::Identity();
  // state.fromPositionOrientationScale(Eigen::Vector3d::Zero(), stata, Eigen::Vector3d::Ones());
  // Eigen::Vector3d translation(-10, 20, 30);
  // printTransformation(state, "Transformation before:   ");

  // state.rotate(euler);
  // state.translate(translation);
  // printTransformation(state, "Rotate->   Translate:    ");

  // state.fromPositionOrientationScale(Eigen::Vector3d::Zero(), stata, Eigen::Vector3d::Ones());
  // // std::cout << "Transformation before xform:\n" << state.rotation() << "\n---" << state.translation().transpose() << std::endl;

  // state.rotate(euler);
  // state.pretranslate(translation);
  // printTransformation(state, "Rotate->   Pretranslate: ");

  // state.fromPositionOrientationScale(Eigen::Vector3d::Zero(), stata, Eigen::Vector3d::Ones());
  // // std::cout << "Transformation before xform:\n" << state.rotation() << "\n---" << state.translation().transpose() << std::endl;

  // state.prerotate(euler);
  // state.pretranslate(translation);
  // printTransformation(state, "Prerotate->Pretranslate: ");

  // Eigen::Transform<double, 3, Eigen::Projective>  projection;
  // projection.setIdentity();

  // state.fromPositionOrientationScale(Eigen::Vector3d::Zero(), stata, Eigen::Vector3d::Ones());

  // state.prerotate(euler);
  // state.translate(translation);
  // printTransformation(state, "Prerotate->   Translate: ");

}

TEST_F(CameraTest, pointInDirection) {
  // Eigen::Vector3d cam_pos(2e3,0,0);
  // Eigen::Vector3d location(0,0,1e3);
  // // TODO: is this backwards?
  // Eigen::Matrix3d xform = getAttitudeTransformBetweenPoints(cam_pos, location);
  // std::cout << "Xform:\n" << xform << std::endl;
  // Eigen::Affine3d proj_mtx = Eigen::Affine3d::Identity();
  // proj_mtx.translate(Eigen::Vector3d(50, 100, 200));
  // std::cout << "Projection matrix:\n" << proj_mtx.translation() << "\n" << proj_mtx.matrix() << std::endl;
  // std::cout << "Looking at " << location.transpose() << "\n" << lookAt(cam_pos, location, Eigen::Vector3d::UnitZ()) << std::endl;
  // std::cout << "Looking at " << location.transpose() << "\n" << lookAt(cam_pos, location) << std::endl;
  // std::cout << "getAttMtx:\n" << getAttitudeTransformBetweenPoints(cam_pos, location) << std::endl;
  // std::cout << "getAttMtx:\n" << getENUFrame(cam_pos - location) << std::endl;
}

TEST_F(CameraTest, Multiplication) {
  // Eigen::Isometry3d state = Eigen::Isometry3d::Identity(), stata = Eigen::Isometry3d::Identity(), statb = Eigen::Isometry3d::Identity();
  // double droll = deg2rad(15), dpitch = deg2rad(8), dyaw = deg2rad(20);
  // double aroll = deg2rad(3), apitch = deg2rad(9), ayaw = deg2rad(4);
  // Eigen::Quaterniond att( Eigen::Matrix3d(Eigen::AngleAxisd(ayaw, Eigen::Vector3d::UnitZ()) * 
  //                           Eigen::AngleAxisd(apitch, Eigen::Vector3d::UnitY()) * 
  //                           Eigen::AngleAxisd(aroll, Eigen::Vector3d::UnitX())));
  // Eigen::Quaterniond euler( Eigen::Matrix3d(Eigen::AngleAxisd(dyaw, Eigen::Vector3d::UnitZ()) * 
  //                           Eigen::AngleAxisd(dpitch, Eigen::Vector3d::UnitY()) * 
  //                           Eigen::AngleAxisd(droll, Eigen::Vector3d::UnitX())));
  // Eigen::Vector3d vec(1,0,0);


  // state.linear() = att.toRotationMatrix();     
  // stata.linear() = att.toRotationMatrix();
  // statb.linear() = att.toRotationMatrix();
  // printTransformation(state, "Original: ");
  // assert(state.isApprox(stata));
  // assert(state.isApprox(statb));
  // std::cout << std::endl;

  // for(uint i = 0; i < 10; i++) {
  //   state.rotate(euler);
  //   stata.prerotate(euler);
  //   statb.linear() = euler.toRotationMatrix().inverse() * statb.linear();
  //   printTransformation(state, "Post : ");
  //   printTransformation(stata, "Prea : ");
  //   printTransformation(statb, "Preb : ");
  //   std::cout << "State: " << (state * vec).transpose() << std::endl;
  //   std::cout << "Stata: " << (stata * vec).transpose() << std::endl;
  //   std::cout << "Statb: " << (statb * vec).transpose() << std::endl;
  //   std::cout << std::endl;
  


  // Test case 1: 90 degrees rotation around Z-axis
  Eigen::Vector3d vec1(1, 0, 0);
  Eigen::Matrix3d rotMatrix1;
  rotMatrix1 = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d expected1(0, 1, 0);
  Eigen::Vector3d result1 = rotMatrix1 * vec1;
  ASSERT_TRUE(result1.isApprox(expected1));

  // Test case 2: 180 degrees rotation around Z-axis
  Eigen::Vector3d vec2(1, 0, 0);
  Eigen::Matrix3d rotMatrix2;
  rotMatrix2 = Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ());
  Eigen::Vector3d expected2(-1, 0, 0);
  Eigen::Vector3d result2 = rotMatrix2 * vec2;
  ASSERT_TRUE(result2.isApprox(expected2));

  // Test case 3: 90 degrees rotation around Y-axis
  Eigen::Vector3d vec3(1, 0, 0);
  Eigen::Matrix3d rotMatrix3;
  rotMatrix3 = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY());
  Eigen::Vector3d expected3(0, 0, -1);
  Eigen::Vector3d result3 = rotMatrix3 * vec3;
  ASSERT_TRUE(result3.isApprox(expected3));

  // Test case 4: No rotation (identity matrix)
  Eigen::Vector3d vec4(1, 2, 3);
  Eigen::Matrix3d rotMatrix4 = Eigen::Matrix3d::Identity();
  Eigen::Vector3d expected4(1, 2, 3);
  Eigen::Vector3d result4 = rotMatrix4 * vec4;
  ASSERT_TRUE(result4.isApprox(expected4));
}

TEST_F(CameraTest, IsInFrame) {
  Eigen::Matrix3d intrinsic;
  // TODO: am I mixing things up by having pretranslate in the camera class? things might be wrong there
  // intrinsic << 800, 0, 320.5,
  //               0, 800, 240.5,
  //               0, 0, 1;
  intrinsic = cam->getIntrinsicMatrix();

  // Example transformation matrix (arbitrary rotation and translation)
  Eigen::Isometry3d transformation = Eigen::Isometry3d::Identity();
  Eigen::AngleAxisd rot = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  // Eigen::AngleAxisd rot1 = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  transformation.rotate(rot);
  transformation.pretranslate(Eigen::Vector3d(1, 0, 0));

  // Example 3D point
  Eigen::Vector3d point(1, 0, 5);

  // Variable to hold the resulting pixel coordinate
  Eigen::Vector2d pixel;

  // Project the point to pixel coordinates
  Eigen::Vector2d img_midpoint(intrinsic(0,2), intrinsic(1,2));
  bool success;

  ASSERT_TRUE((success=projectPointToPixel(intrinsic, transformation, point, pixel)));
  ASSERT_FALSE(img_midpoint.isApprox(pixel));
  std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
            << (success ? " ":"NOT ") << "in the image" << std::endl;
            
  transformation.pretranslate(Eigen::Vector3d(0, 0, -5));
  // Project the point to pixel coordinates
  ASSERT_TRUE((success=projectPointToPixel(intrinsic, transformation, point, pixel)));
  ASSERT_FALSE(img_midpoint.isApprox(pixel));
  std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
            << (success ? "":"NOT ") << "in the image" << std::endl;

  transformation.rotate(rot);            
  // Project the point to pixel coordinates
  ASSERT_FALSE((success=projectPointToPixel(intrinsic, transformation, point, pixel)));
  ASSERT_FALSE(img_midpoint.isApprox(pixel));
  std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
            << (success ? "":"NOT ") << "in the image" << std::endl;


  transformation.pretranslate(Eigen::Vector3d(16, 4, 0));
  Eigen::Matrix3d look_at = lookAt(transformation, point, Eigen::Vector3d::UnitZ());
  
  Eigen::Matrix3d mtx = getAttitudeTransformBetweenPoints(transformation.translation(), point);
  // transformation.linear() = mtx.transpose();
  transformation.linear() = look_at;
  ASSERT_TRUE((success=projectPointToPixel(intrinsic, transformation, point, pixel)));
  ASSERT_TRUE(img_midpoint.isApprox(pixel));
  std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
            << (success ? " ":"NOT ") << "in the image" << std::endl;

  transformation.linear() = Eigen::Matrix3d::Identity();
  transformation.rotate(mtx.transpose());
  ASSERT_TRUE((success=projectPointToPixel(intrinsic, transformation, point, pixel)));
  ASSERT_TRUE(img_midpoint.isApprox(pixel));
  std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
            << (success ? " ":"NOT ") << "in the image" << std::endl;
}
