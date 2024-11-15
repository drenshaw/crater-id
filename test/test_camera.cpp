#include "camera.h"
#include "quadrics.h"
#include "conics.h"
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
  
  Eigen::Vector3d cam_pos = 1e4*latlon2bearing(-15, 0);
  Eigen::Vector3d look_at_me = 1e4*Eigen::Vector3d::Zero();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->setPosition(cam_pos);
  cam->pointTo(look_at_me, up_vector);
  
  Eigen::MatrixXd proj_mtx    = cam->getProjectionMatrix();

  Quadric quad(-5, 0, 200, "south_pole");
  Eigen::Matrix4d q_envelope = quad.getEnvelope();
  Eigen::Matrix3d c_envelope = proj_mtx * q_envelope * proj_mtx.transpose();
  Conic con(adjugate(c_envelope));
  Conic conic = quad.projectToImage(cam->getProjectionMatrix());

  ASSERT_EQ(con, conic);
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

TEST_F(CameraTest, GetImageMidpoint) {
  Eigen::Matrix3d intrinsic = cam->getIntrinsicMatrix();
  Eigen::Vector2d img_midpoint1(intrinsic(0,2), intrinsic(1,2));
  Eigen::Vector2d img_midpoint = cam->getImageMidpoint();
  ASSERT_TRUE(img_midpoint1.isApprox(img_midpoint));
}

TEST_F(CameraTest, InverseIntrinsicMatrix) {
  Eigen::Matrix3d K = cam->getIntrinsicMatrix();
  Eigen::Matrix3d Kinv = cam->getInverseIntrinsicMatrix();

  ASSERT_TRUE(K.inverse().isApprox(Kinv));
}

TEST_F(CameraTest, PointMovingInRightDirection) {
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

TEST_F(CameraTest, PointJumpsToWrongPixel) {
  
  Eigen::Vector3d cam_pos = 1e4*latlon2bearing(-15, 0);
  Eigen::Vector3d look_at_me = 1e4*Eigen::Vector3d::Zero();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->setPosition(cam_pos);
  cam->pointTo(look_at_me, up_vector);
  for(int i = 0; i < 10; i++) {
    Eigen::Vector3d pt = {i*2.5e3,0,0};
    Eigen::Vector2d pixel;
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
  Eigen::Isometry3d transform;
  double rx = 0, ry = 0, rz = 10, rw = 10;
  Eigen::Quaterniond rotation(Eigen::Quaterniond(rw, rx, ry, rz).normalized());
  Eigen::Translation3d translation(Eigen::Vector3d(1,2,3));
  transform = rotation * translation;
  Eigen::Vector3d vec(0.1, 0.2, 100);
  transform.translation() = vec;
  ASSERT_TRUE(transform.translation().isApprox(vec));

  Eigen::Quaterniond qrand = Eigen::Quaterniond::UnitRandom();
  Eigen::Matrix3d mrand = qrand.toRotationMatrix();
  Eigen::Vector3d v_rand = Eigen::Vector3d::Random();
  ASSERT_TRUE((qrand*v_rand).isApprox(mrand*v_rand));
}

TEST_F(CameraTest, EnsurePassiveXform) {
  // TODO: wait for Paul to send me unit tests for this one
}

void printTransformation(const Eigen::Isometry3d& xform, const std::string& str) {
  std::streamsize ss = std::cout.precision();
  std::streamsize sw = std::cout.width();
  Eigen::Quaterniond quat = Eigen::Quaterniond(xform.rotation());
  std::cout << str << std::fixed 
            << std::setw(8) << std::setprecision(1)
            << xform.translation().transpose()
            << std::setw(8) << std::setprecision(5)
            // << std::endl 
            << quat.x() << ", "
            << quat.y() << ", "
            << quat.z() << ", "
            << quat.w()
            << std::setprecision(ss) << std::setw(sw) << std::defaultfloat
            << std::endl;
}

TEST_F(CameraTest, MoveRelative) {

}

TEST_F(CameraTest, pointInDirection) {
  Eigen::Vector3d cam_pos(2e3,0,0);
  Eigen::Vector3d location(0,0,1e3);
  Eigen::Vector2d pixel;
  Eigen::Vector2d cam_mid = Eigen::Vector2d(1296.5, 1024.5);
  double prec = 1e-3;
  // // TODO: is this backwards?
  cam->move(cam_pos);
  cam->pointTo(location, Eigen::Vector3d::UnitZ());
  cam->world2Pixel(location, pixel);
  EXPECT_TRUE(pixel.isApprox(cam_mid, prec));
  // std::cout << "Pointing at ( " << location.transpose() << " ) which becomes: ( " << pixel.transpose() << " ) in the camera view\n";

  location = Eigen::Vector3d(0,0,0.5e3);
  cam->world2Pixel(location, pixel);
  EXPECT_TRUE(pixel.isApprox(Eigen::Vector2d(1296.5, 1064.15), prec));
  // std::cout << "Pointing at ( " << location.transpose() << " ) which becomes: ( " << pixel.transpose() << " ) in the camera view\n";

  location = Eigen::Vector3d::Zero();
  cam->setPosition(Eigen::Vector3d(0,0,-1e4));
  cam->pointTo(location, -Eigen::Vector3d::UnitY());
  cam->world2Pixel(location, pixel);
  EXPECT_TRUE(pixel.isApprox(cam_mid, prec));
  // std::cout << "Pointing at ( " << location.transpose() << " ) which becomes: ( " << pixel.transpose() << " ) in the camera view\n";
  
  location = -1e3*Eigen::Vector3d::UnitZ();
  cam->world2Pixel(location, pixel);
  EXPECT_TRUE(pixel.isApprox(cam_mid, prec));
  // std::cout << "Pointing at ( " << location.transpose() << " ) which becomes: ( " << pixel.transpose() << " ) in the camera view\n";

  location = -1e3*Eigen::Vector3d::UnitX();
  cam->world2Pixel(location, pixel);
  EXPECT_TRUE(pixel.isApprox(Eigen::Vector2d(1196.5, 1024.5), prec));
  // std::cout << "Pointing at ( " << location.transpose() << " ) which becomes: ( " << pixel.transpose() << " ) in the camera view\n";
}

TEST_F(CameraTest, ProjectQuadricAfterRotation) {

  Eigen::Matrix3d rotMatrix1 = Eigen::AngleAxisd(-M_PI / 6, Eigen::Vector3d::UnitZ()).toRotationMatrix();
  cam->rotate(rotMatrix1);
  // std::cout << *cam << std::endl;
  Eigen::Matrix3d intrinsic = cam->getIntrinsicMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();
  Quadric quad(-70, 0, 100, "rotation@SouthPole");
  Eigen::MatrixXd proj_mtx = intrinsic * extrinsic;
  Eigen::MatrixXd proj1 = cam->getProjectionMatrix();
  EXPECT_TRUE(proj_mtx.isApprox(proj1));
  Eigen::Matrix3d conic_envelope = proj_mtx * quad.getEnvelope() * proj_mtx.transpose();
  Eigen::Matrix3d conic_locus = adjugate(conic_envelope);
  Conic conic(conic_locus);
  // std::cout << conic << std::endl;

  Eigen::Matrix3d c_envelope = cam->projectQuadric(quad.getEnvelope());
  Eigen::Matrix3d c_locus = adjugate(c_envelope);
  Conic c(c_locus);
  // std::cout << c << std::endl;

  EXPECT_EQ(conic, c);

  // TODO: ensure that this active rotation results in the proper frame transformation
  // Then figure out the projection matrix issue
}

TEST_F(CameraTest, Rotations) {
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
  // TODO: am I mixing things up by having pretranslate in the camera class? things might be wrong there
  Eigen::Vector2d img_midpoint = cam->getImageMidpoint();

  // Example transformation matrix (arbitrary rotation and translation)
  Eigen::AngleAxisd rot = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  // Eigen::AngleAxisd rot1 = Eigen::AngleAxisd(M_PI / 4, Eigen::Vector3d::UnitY());
  Eigen::Vector3d offset(1, 0, 0);
  cam->rotate(rot);
  cam->setPosition(offset);

  // Example 3D point
  Eigen::Vector3d point(1, 0, 5);

  // Variable to hold the resulting pixel coordinate
  Eigen::Vector2d pixel, other_pxl(296.5, 1024.5);

  // Project the point to pixel coordinates
  bool success;

  success = cam->isInCameraFrame(point, pixel);
  ASSERT_TRUE(success);
  ASSERT_TRUE(other_pxl.isApprox(pixel));
  other_pxl(296.5,1024.5);

  Eigen::Vector3d translation1(0, 0, -5);
  cam->move(translation1);
  // Project the point to pixel coordinates
  ASSERT_TRUE(cam->isInCameraFrame(point, pixel));
  ASSERT_TRUE(other_pxl.isApprox(pixel));
  // std::cout << "Pixel coordinates: " << pixel.transpose() << " are " 
  //           << (success ? "":"NOT ") << "in the image" << std::endl;
  
  cam->rotate(rot);            
  // Project the point to pixel coordinates
  EXPECT_THROW(cam->isInCameraFrame(point, pixel), std::runtime_error);

  Eigen::Vector3d translation2(16, 4, 0);
  cam->move(translation2);
  
  cam->pointTo(point, Eigen::Vector3d::UnitZ());
  EXPECT_TRUE(cam->isInCameraFrame(point, pixel));
  EXPECT_TRUE(img_midpoint.isApprox(pixel));

  Eigen::AngleAxisd rot1 = Eigen::AngleAxisd(M_PI / 12, Eigen::Vector3d::UnitY());
  cam->rotate(rot1);
  Eigen::Vector2d rot1_exp_pixel(1028.55, 1024.5);
  EXPECT_TRUE(cam->isInCameraFrame(point, pixel));
  EXPECT_TRUE(rot1_exp_pixel.isApprox(pixel, 1e-5));

  cam->resetCameraState();
  cam->pointTo(point, Eigen::Vector3d::UnitZ());
  ASSERT_TRUE(cam->isInCameraFrame(point, pixel));
  ASSERT_TRUE(img_midpoint.isApprox(pixel));
  
}

TEST_F(CameraTest, ProjectQuadric) {
  Quadric quad1( 10,  15, 50, "Crater 1");
  Quadric quad2(-10,  15, 100, "Crater 2");
  Quadric quad3(-30, -10, 75, "Crater 3");
  Quadric quad4(  0, -20, 150, "Crater 4");

  Eigen::Vector3d up_axis(0,0,1);
  cam->move(Eigen::Vector3d(1e4, 0, 0));
  cam->pointTo(quad1.getLocation(), up_axis);
  

  // Eigen::Matrix3d c_envelope = cam->projectQuadric(quad1.getEnvelope());
  // Eigen::Matrix3d c_locus = adjugate(c_envelope);
  // std::cout << "Image Ellipse: " << Conic(c_locus) << std::endl;
  // std::cout << "Image Plane Ellipse: " << Conic(cam->getImagePlaneLocus(c_locus)) << std::endl;
}

TEST_F(CameraTest, MoveCamera) {
  cam->resetCameraState();
  Eigen::AngleAxisd rot1 = Eigen::AngleAxisd(M_PI / 12, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd rot2 = Eigen::AngleAxisd(M_PI / 5, Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd rot3 = Eigen::AngleAxisd(15*M_PI / 7, Eigen::Vector3d::UnitX());
  Eigen::Vector3d move1(10,20,30);
  Eigen::Vector3d move2(37,50,300);
  Eigen::Vector3d move3(-110,-20,0);
  Eigen::Vector3d move4(-80,-30,-230);
  cam->move(move1);
  ASSERT_TRUE(cam->getPosition().isApprox(move1));
  cam->rotate(rot1);

  cam->move(move2);
  ASSERT_TRUE(cam->getPosition().isApprox(move1+move2));
  cam->rotate(rot2);

  cam->move(move3);
  ASSERT_TRUE(cam->getPosition().isApprox(move1+move2+move3));
  cam->rotate(rot3);

  cam->move(move4);
  ASSERT_TRUE(cam->getPosition().isApprox(move1+move2+move3+move4));
  
}

TEST_F(CameraTest, RotateCamera) {
  cam->resetCameraState();
  Eigen::AngleAxisd rot1 = Eigen::AngleAxisd( M_PI / 12, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd att1 = Eigen::AngleAxisd(-M_PI / 12, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd rot2 = Eigen::AngleAxisd( M_PI / 5, Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd att2 = Eigen::AngleAxisd(-M_PI / 5, Eigen::Vector3d::UnitZ());
  Eigen::AngleAxisd rot3 = Eigen::AngleAxisd( 15*M_PI / 7, Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd att3 = Eigen::AngleAxisd(-15*M_PI / 7, Eigen::Vector3d::UnitX());

  // One rotation
  cam->rotate(rot1);
  Eigen::Quaterniond q_att = cam->getAttitude();
  ASSERT_TRUE(q_att.isApprox(Eigen::Quaterniond(att1)));

  // Two rotations
  cam->resetCameraState();
  cam->rotate(rot1);
  cam->rotate(rot2);
  
  Eigen::Quaterniond att12 = Eigen::Quaterniond(att2*att1);
  Eigen::Quaterniond q_att2 = cam->getAttitude();
  cam->resetCameraState();
  cam->rotate(rot1*rot2);
  Eigen::Quaterniond q_att2_mult = cam->getAttitude();
  ASSERT_TRUE(att12.isApprox(q_att2));
  // TODO: Figure out why this is incorrect in the first vector element
  ASSERT_TRUE(q_att2_mult.isApprox(q_att2) || q_att2_mult.coeffs().isApprox(-q_att2.coeffs()));

  // Three rotations
  cam->resetCameraState();
  cam->rotate(rot1);
  cam->rotate(rot2);
  cam->rotate(rot3);
  
  Eigen::Quaterniond att123 = Eigen::Quaterniond(att3*att2*att1);
  Eigen::Quaterniond q_att3 = cam->getAttitude();
  cam->resetCameraState();
  cam->rotate(rot1*rot2*rot3);
  Eigen::Quaterniond q_att3_mult = cam->getAttitude();
  // Eigen does not handle equality of negative quaternions
  ASSERT_TRUE(att123.isApprox(q_att3) || att123.coeffs().isApprox(-q_att3.coeffs()));
  // TODO: Figure out why this is incorrect in the first vector element
  ASSERT_TRUE(q_att3_mult.isApprox(q_att3));
}

TEST_F(CameraTest, SetAttitude) {
  cam->resetCameraState();
  Eigen::AngleAxisd rot1 = Eigen::AngleAxisd( M_PI / 12, Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd att1 = Eigen::AngleAxisd(-M_PI / 12, Eigen::Vector3d::UnitY());

  cam->rotate(rot1);
  Eigen::Quaterniond res_rot = cam->getAttitude();

  cam->resetCameraState();
  cam->setAttitude(att1);
  Eigen::Quaterniond res_att = cam->getAttitude();

  ASSERT_TRUE(res_att.isApprox(res_rot));
}

TEST_F(CameraTest, ProjectionMatrix) {
  cam->resetCameraState();
  Eigen::Vector3d move1(10,20,30);
  cam->move(move1);
  cam->pointTo(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
  Eigen::Quaterniond att = cam->getAttitude();
  Eigen::Vector3d pos = cam->getPosition();
  Eigen::MatrixXd ext(3,4);

  ext.topLeftCorner(3,3) = Eigen::Matrix3d::Identity();
  ext.topRightCorner(3,1) = -pos;
  ext = att.toRotationMatrix() * ext;
  // Eigen::AffineCompact3d proj;
  // std::cout << "Transformation matrix:\n" << att.toRotationMatrix() << std::endl;
  // std::cout << "Extrinsic matrix:\n" << ext << std::endl;
}
