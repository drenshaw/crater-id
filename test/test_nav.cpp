#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp>
#include <random>

#include "navigation.h"
#include "camera.h"
#include "quadrics.h"
#include "conics.h"
#include "visuals.h"

class NavigationTest : public testing::Test {
protected:
  NavigationTest() {
    latitude = -30, longitude = 0, radius = 200;
    id = "defaultQuadric";
    quadric_default = new Quadric(latitude, longitude, radius, id);
    dx = 1000, dy = 1000, skew = 0, im_height = 2048, im_width = 2592;
    up = (im_width+1)/2;
    vp = (im_height+1)/2;
    image_size = cv::Size(im_width, im_height);
    cam = new Camera(dx, dy, up, vp, skew, image_size, quat, position);
    cv_cam = new cv::viz::Camera(dx, dy, up, vp, image_size);
  }
  ~NavigationTest() override {
    // delete quadric_default;
      delete cam;
      delete cv_cam;
  }
  public:
  double latitude, longitude, radius;
  std::string id;
  Quadric* quadric_default;
    double dx, dy, skew, im_height, im_width, up, vp;
    cv::Size2i image_size;
    Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();
    Eigen::Vector3d position{1e4, 0, 0};

    Camera* cam;
    cv::viz::Camera* cv_cam;

};

TEST_F(NavigationTest, ShiuConicBackprojection) {
  // // Camera cam;
  // // Quadric quad(-30,0,100,"backprojection");
  // // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  // Eigen::Vector3d location = R_MOON*Eigen::Vector3d::UnitZ();
  // Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // // cam->resetCameraState();
  // cam->moveX(1e4);
  // cam->pointTo(location, up_vector);
  
  // Eigen::MatrixXd proj = cam->getProjectionMatrix();
  // Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();
  // // Eigen::Matrix3d Kinv = cam->getInverseIntrinsicMatrix();
  
  // std::array<Eigen::Vector3d, 2> centers, normals;
  // Eigen::Vector3d center1, center2, normal1, normal2;
  // // std::cout << *quadric_default << std::endl;
  // // std::cout << cam << std::endl;

  // // Eigen::Matrix3d c_locus = cam->projectQuadricToLocus(quadric_default->getLocus());
  // // Eigen::Matrix3d plane_locus = Kinv * c_locus * Kinv.transpose();
  // // double dist;
  // // double dist2center = (quad.getLocation() - cam->getPosition()).norm();
  // // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  // Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(extrinsic);
  // // std::cout << "Plane:\n" << quadric_default->getLocus() << std::endl;
  // // std::cout << "Plane1:\n" << conic_locus << std::endl;
  // Conic cc(conic_locus);
  // std::cout << cc << std::endl;
  // Shiu::conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
}

TEST_F(NavigationTest, ChristianBackprojection) {
  // Camera cam;
  // Quadric quad(-30,0,100,"backprojection");
  // Eigen::Vector3d location = Eigen::Vector3d::Zero();
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->resetCameraState();
  cam->moveX(1e4);
  cam->moveY(-1e3);
  cam->moveZ(-2e3);
  cam->pointTo(look_here, up_vector);
  
  // Eigen::MatrixXd proj = cam->getProjectionMatrix();
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Eigen::Matrix3d att = cam->getAttitudeMatrix();
  Eigen::Vector3d centerWrtCamWorld = quadric_default->getLocation() - cam->getPosition();
  Eigen::Vector3d centerWrtCam = att * centerWrtCamWorld;
  // double dist2center = (centerWrtCamWorld).norm();
  Eigen::Vector3d normal_test = att * quadric_default->getNormal();
  // std::cout << "---\n---Truth Data---\n";
  // std::cout << "---DISTANCE TO QUADRIC CENTER: " << dist2center << std::endl;
  // std::cout << "---CAMERA TO QUADRIC CENTER: " << centerWrtCam.transpose() << std::endl;
  // std::cout << "---QUADRIC NORMAL WRT CAMERA: " << normal_test.transpose() << std::endl;
  
  std::array<Eigen::Vector3d, 2> centers, normals;
  Eigen::Vector3d center1, center2, normal1, normal2;
  // double dist;
  // Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus = quadric_default->projectToPlaneLocus(extrinsic);
  Christian::conicBackprojection(conic_locus, quadric_default->getRadius(), centers, normals);
  center1 = centers.at(0);
  center2 = centers.at(1);
  normal1 = normals.at(0);
  normal2 = normals.at(1);
  // std::cout << "+++\n+++Calculations+++\n";
  // std::cout << "Center (cam) 1: " << center1.transpose() << std::endl;
  // std::cout << "Center (cam) 2: " << center2.transpose() << std::endl;
  // std::cout << "Normal (cam) 1: " << normal1.transpose() << std::endl;
  // std::cout << "Normal (cam) 2: " << normal2.transpose() << std::endl;
  // Eigen::Vector3d q_normal = Eigen::Vector3d::UnitZ();
  // std::cout << "Angle 1: " << rad2deg(std::acos(q_normal.dot(-normal1))) << std::endl;
  // std::cout << "Angle 2: " << rad2deg(std::acos(q_normal.dot(-normal2))) << std::endl;
  ASSERT_TRUE(centerWrtCam.isApprox(center1, 1e-3) || centerWrtCam.isApprox(center2, 1e-3));
  ASSERT_TRUE( normal_test.isApprox(normal1, 1e-3) ||  normal_test.isApprox(normal2, 1e-3));
  // conicBackprojection(conic_locus, quad.getRadius(), normal, dist);
}

TEST_F(NavigationTest, ChooseSupportingPlanes) {
  cam->moveX(1e4);
  // cam->moveY(-1e3);
  cam->moveZ(-2e3);
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->pointTo(look_here, up_vector);
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Quadric q1(30, -30, 100, "crater1");
  Quadric q2(-10, 30, 200, "crater2");
  Eigen::Matrix3d conic_locus1 = q1.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus2 = q2.projectToPlaneLocus(extrinsic);
  std::array<Eigen::Vector3d, 2> centers1, normals1, centers2, normals2;
  Christian::conicBackprojection(conic_locus1, q1.getRadius(), centers1, normals1);
  Christian::conicBackprojection(conic_locus2, q2.getRadius(), centers2, normals2);
  double angle = q1.getAngleBetweenQuadrics(q2);
  // std::cout << "Angle between q1 and q2: " << angle << " (" << rad2deg(angle) << ")\n";
  uint planes_idx = chooseSupportingPlanes(angle, normals1, normals2);

  uint idx1 = planes_idx / 2;
  uint idx2 = planes_idx % 2;
  Eigen::Vector3d normal1 = normals1.at(idx1);
  Eigen::Vector3d normal2 = normals2.at(idx2);
  double angle_calc = getAngleBetweenVectors(normal1, normal2);
  ASSERT_NEAR(angle, angle_calc, 1e-8);
}

TEST_F(NavigationTest, ChooseVectorOfCraters) {
  cam->moveX(1e4);
  // cam->moveY(-1e3);
  cam->moveZ(-2e3);
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->pointTo(look_here, up_vector);
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Quadric q1( 30, -30, 100, "crater1");
  Quadric q2(-10,  30, 200, "crater2");
  Quadric q3( 10,   0,  50, "crater3");
  Quadric q4(-10,   0,  30, "crater3");
  std::vector<Quadric> quadrics = {q1, q2, q3, q4};
  std::vector<Eigen::Matrix3d> locii;
  for(std::vector<Quadric>::iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    locii.push_back((*it).projectToPlaneLocus(extrinsic));
  }

  std::vector<Eigen::Vector3d> centers, normals;
  reprojectLociiToQuadrics(quadrics, locii, centers, normals);

  Eigen::Matrix3d att = cam->getAttitudeMatrix();
  // std::vector<Eigen::Vector3d>::iterator the_chosen;
  std::vector<Eigen::Vector3d>::iterator the_chosen;
  for(the_chosen = normals.begin(); the_chosen != normals.end(); the_chosen++) {
    int index = getIndex(normals.begin(), the_chosen);
    Quadric quad = quadrics.at(index);
    Eigen::Vector3d centerWrtCamWorld = quad.getLocation() - cam->getPosition();
    Eigen::Vector3d centerWrtCam = att * centerWrtCamWorld;
    // double dist2center = (centerWrtCamWorld).norm();
    Eigen::Vector3d normalWrtCam = att * quad.getNormal();
    // std::cout << "---\n---Truth Data---\n";
    // // std::cout << "---DISTANCE TO QUADRIC CENTER: " << dist2center << std::endl;
    // std::cout << "---\tCENTER: " << centerWrtCam.transpose() << std::endl;
    // std::cout << "---\tNORMAL: " << normalWrtCam.transpose() << std::endl;
    // // std::cout << "Crater " << index << std::endl;
    Eigen::Vector3d center = centers.at(index);
    // std::cout << "---Calculated Data---\n";
    // std::cout << "\tCenter: " << center.transpose() << std::endl;
    // std::cout << "\tNormal: " << (*the_chosen).transpose() << std::endl;
    EXPECT_TRUE(centerWrtCam.isApprox(center, EPS));
    EXPECT_TRUE(normalWrtCam.isApprox(*the_chosen, EPS));
  }
}

TEST_F(NavigationTest, CalculateHomography) {
  cam->resetCameraState();
  cam->moveX(2.5e3);
  cam->moveY(-1e2);
  cam->moveZ(-2e2);
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->pointTo(look_here, up_vector);
  // std::cout << *cam << std::endl;
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Quadric q1( 10, -10, 100, "crater1");
  Quadric q2(-10,  10, 200, "crater2");
  Quadric q3( 10,   0,  50, "crater3");
  Quadric q4(-10,   0,  30, "crater4");
  Quadric q5(  0,  -3, 150, "crater5");
  std::vector<Quadric> quadrics = {q1, q2, q3, q4, q5};

  // cv::Mat image = cam->getBlankCameraImage();
  // cv::Mat outImg;
  // double scaling = 0.4;
  // // Showing image inside a window 
  // viz::drawEllipses(image, *cam, quadrics, viz::CV_colors);
  // // viz::drawEllipses(image, conics, viz::CV_colors);
  // cv::resize(image, outImg, cv::Size(), scaling, scaling);
  // cv::imshow("Projecting Moon to Camera", outImg); 
  // cv::waitKey(0); 

  std::vector<Eigen::Matrix3d> locii;
  for(std::vector<Quadric>::iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    locii.push_back((*it).projectToPlaneLocus(extrinsic));
    std::cout << Conic((*it).projectToPlaneLocus(extrinsic)) << std::endl;
  }

  Eigen::Quaterniond attitude;
  Eigen::Vector3d position;
  solve_navigation_problem(quadrics, locii, attitude, position);
  ASSERT_TRUE(cam->getAttitude().isApprox(attitude, 1e-3));
  ASSERT_TRUE(cam->getPosition().isApprox(position, 1e-3));
}

TEST_F(NavigationTest, KanataniEllipseCheck) {
  // double semimajor = 100, semiminor = 40, xc = 40, yc = 120, phi = deg2rad(10);
  // Conic conic(semimajor, semiminor, xc, yc, phi);
}

TEST_F(NavigationTest, ImageConic2PlaneConic) {
  // Quadric quad0( 0,  0, 100, "Crater 1");
  // // Quadric quad1( 10,  15, 50, "Crater 1");
  // // Quadric quad2(-10,  15, 100, "Crater 2");
  // // Quadric quad3(-30, -10, 75, "Crater 3");
  // // Quadric quad4(  0, -20, 150, "Crater 4");

  // Eigen::Vector3d up_axis(0,0,1);
  // cam->move(Eigen::Vector3d(1e4, 0, 0));
  // cam->pointTo(quad0.getLocation(), up_axis);
  

  // Eigen::Matrix3d c_envelope = cam->projectQuadric(quad0.getEnvelope());
  // Eigen::Matrix3d c_locus = adjugate(c_envelope);
  // Conic image_conic(c_locus);
  // Conic plane_conic(cam->getImagePlaneLocus(c_locus));
  // Eigen::Vector2d ellipse_center = image_conic.getCenter();
  // Eigen::Vector2d plane_center = plane_conic.getCenter();
  // Eigen::Vector2d image_center = cam->getImageMidpoint();
  // Eigen::Vector2d zeros = Eigen::Vector2d::Zero();
  // ASSERT_TRUE(image_center.isApprox(ellipse_center, 1e-8));
  // // TODO: For some reason the next line fails, but the subsequent line passes
  // // ASSERT_TRUE(zeros.isApprox(plane_center, 1e-8));
  // ASSERT_LE((plane_center-zeros).norm(), 1e-9);
}

void addNoise(const double mean, const double st_dev, std::vector<Eigen::Vector2d>& points) {
  if(mean == 0 && st_dev == 0) {
    return;
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> dist(mean, st_dev); // Mean 0, standard deviation 1

  // Create an Eigen matrix and fill it with noise
  Eigen::MatrixXd A(3, 3);
  // for (int i = 0; i < A.rows(); ++i) {
  for(std::vector<Eigen::Vector2d>::iterator it = points.begin(); it != points.end(); it++) {
    // std::cout << "Points: " << (*it).transpose();
    (*it)(0) += dist(gen);
    (*it)(1) += dist(gen);
    // std::cout << "\tAdded noise: " << (*it).transpose() << std::endl;
  }
}

TEST_F(NavigationTest, QuadricPointsWNoise) {
  int n_pts = 10;

  cam->resetCameraState();
  cam->moveX(2.5e3);
  cam->moveY(-1e2);
  cam->moveZ(-2e2);
  Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  cam->pointTo(look_here, up_vector);
  // std::cout << *cam << std::endl;
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Quadric q1( 10, -10, 100, "crater1");
  Quadric q2(-10,  10, 200, "crater2");
  Quadric q3( 10,   0,  50, "crater3");
  Quadric q4(-10,   0,  30, "crater4");
  Quadric q5(  0,  -3, 150, "crater5");
  std::vector<Quadric> quadrics = {q1, q2, q3, q4, q5};

  double mean = 0, st_dev = 0.0;
  std::vector<Eigen::Vector2d> all_pts_pxl;
  all_pts_pxl.reserve(n_pts * quadrics.size());
  std::vector<Eigen::Matrix3d> locii;
  for(std::vector<Quadric>::const_iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    std::vector<Eigen::Vector3d> pts_world;
    pts_world.reserve(n_pts);
    (*it).getRimPoints(n_pts, pts_world);
    std::vector<Eigen::Vector2d> pts_pxl;
    cam->world2Pixel(pts_world, pts_pxl);

    // Conic no_noise(pts_pxl);
    // // std::cout << "====>No noise: " << no_noise << std::endl;
    // addNoise(mean, st_dev, pts_pxl);
    // Conic from_proj = (*it).projectToImage(cam->getExtrinsicMatrix());
    // std::cout << "+++++From projection: " << from_proj << std::endl;
    Conic noisy(pts_pxl);
    // conics.push_back(noisy);
    Eigen::Matrix3d plane_locus = cam->getImagePlaneLocus(noisy.getLocus());
    Conic cc(plane_locus);
    std::cout << "<====Noisy: " << cc << std::endl << std::endl;
    std::cout << "Angle: " << cc.getAngle() << " (" << cc.getAngleDeg() << ") " << M_PI_2 << std::endl;
    locii.push_back(plane_locus);
    for(std::vector<Eigen::Vector2d>::iterator it2 = pts_pxl.begin(); it2 != pts_pxl.end(); it2++){
      all_pts_pxl.push_back(*it2);
    }
  }

  std::vector<Eigen::Matrix3d> locii2;
  for(std::vector<Quadric>::iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    locii2.push_back((*it).projectToPlaneLocus(extrinsic));
    std::cout << Conic((*it).projectToPlaneLocus(extrinsic)) << std::endl;
  }

  // cv::Mat image = cam->getBlankCameraImage();
  // cv::Mat outImg;
  // double scaling = 0.4;
  // // Showing image inside a window 
  // viz::drawEllipses(image, *cam, quadrics, viz::CV_colors);
  // // viz::drawEllipses(image, conics, viz::CV_colors);
  // cv::resize(image, outImg, cv::Size(), scaling, scaling);
  // cv::imshow("Projecting Moon to Camera", outImg); 
  // cv::waitKey(0); 

  Eigen::Quaterniond attitude;
  Eigen::Vector3d position;
  solve_navigation_problem(quadrics, locii2, attitude, position);
  ASSERT_TRUE(cam->getAttitude().isApprox(attitude, 1e-3));
  ASSERT_TRUE(cam->getPosition().isApprox(position, 1e-3));
  // int n_pts = 10;
  // cam->resetCameraState();
  // cam->moveX(2.5e3);
  // cam->moveY(-1e2);
  // cam->moveZ(-2e2);
  // Eigen::Vector3d look_here = -1e2*Eigen::Vector3d::UnitZ();
  // Eigen::Vector3d up_vector = Eigen::Vector3d::UnitZ();
  // cam->pointTo(look_here, up_vector);
  // Quadric q1( 10, -10, 100, "crater1");
  // Quadric q2(-10,  10, 200, "crater2");
  // Quadric q3( 10,   0,  50, "crater3");
  // Quadric q4(-10,   0,  30, "crater4");
  // Quadric q5(  0,  -3, 150, "crater5");
  // const std::vector<Quadric> quadrics = {q1, q2, q3, q4, q5};
  // std::vector<Conic> conics;

  // double mean = 0, st_dev = 0.0;
  // std::vector<Eigen::Vector2d> all_pts_pxl;
  // all_pts_pxl.reserve(n_pts * quadrics.size());
  // std::vector<Eigen::Matrix3d> locii;
  // for(std::vector<Quadric>::const_iterator it = quadrics.begin(); it != quadrics.end(); it++) {
  //   std::vector<Eigen::Vector3d> pts_world;
  //   pts_world.reserve(n_pts);
  //   (*it).getRimPoints(n_pts, pts_world);
  //   std::vector<Eigen::Vector2d> pts_pxl;
  //   cam->world2Pixel(pts_world, pts_pxl);

  //   Conic no_noise(pts_pxl);
  //   // std::cout << "====>No noise: " << no_noise << std::endl;
  //   addNoise(mean, st_dev, pts_pxl);
  //   // Conic from_proj = (*it).projectToImage(cam->getExtrinsicMatrix());
  //   // std::cout << "+++++From projection: " << from_proj << std::endl;
  //   Conic noisy(pts_pxl);
  //   conics.push_back(noisy);
  //   Eigen::Matrix3d plane_locus = cam->getImagePlaneLocus(noisy.getLocus());
  //   std::cout << "<====Noisy: " << Conic(plane_locus) << std::endl << std::endl;
  //   locii.push_back(plane_locus);
  //   for(std::vector<Eigen::Vector2d>::iterator it2 = pts_pxl.begin(); it2 != pts_pxl.end(); it2++){
  //     all_pts_pxl.push_back(*it2);
  //   }
  // }

  // // Eigen::Quaterniond attitude;
  // // Eigen::Vector3d position;
  // // solve_navigation_problem(quadrics, locii, attitude, position);
  // // EXPECT_TRUE(cam->getAttitude().isApprox(attitude, 1e-3));
  // // EXPECT_TRUE(cam->getPosition().isApprox(position, 1e-3));
  
  // cv::Mat image = cam->getBlankCameraImage();
  // // Add the moon
  // Eigen::Matrix4d sphere = makeSphere(double(R_MOON));
  // Eigen::Matrix3d moon_locus = cam->projectQuadricToLocus(sphere);
  // Conic moon(moon_locus);
  // // std::cout << "Moon ellipse: " << moon << std::endl;
  // if(cam->isInCameraFrame(moon.getCenter())) {
  //   viz::drawEllipse(image, moon, cv::viz::Color::gray());
  // }
  // else {
  //   std::cout << "The Moon center is not in the image: " << moon.getCenter().transpose() << std::endl;
  // }
  // viz::draw3dAxes(image, *cam);
  // viz::drawPoints(image, all_pts_pxl, viz::CV_colors);
  // // Showing image inside a window 
  // double scaling = 0.5;
  // cv::Mat outImg;
  // cv::resize(image, outImg, cv::Size(), scaling, scaling);
  // viz::interactiveZoom(outImg);
}

