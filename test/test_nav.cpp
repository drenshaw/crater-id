#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

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

void reprojectLociiToQuadrics(const Camera& camera, 
                              const Quadric& q1, 
                              const Quadric& q2, 
                              const Quadric& q3) {
  Eigen::MatrixXd extrinsic = camera.getExtrinsicMatrix();
  Eigen::Matrix3d conic_locus1 = q1.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus2 = q2.projectToPlaneLocus(extrinsic);
  Eigen::Matrix3d conic_locus3 = q3.projectToPlaneLocus(extrinsic);
  std::array<Eigen::Vector3d, 2> centers1, normals1, centers2, normals2, centers3, normals3;
  double radius1 = q1.getRadius();
  double radius2 = q2.getRadius();
  double radius3 = q3.getRadius();
  double angle12 = q1.getAngleBetweenQuadrics(q2);
  double angle23 = q2.getAngleBetweenQuadrics(q3);
  double angle31 = q3.getAngleBetweenQuadrics(q1);

  Christian::conicBackprojection(conic_locus1, radius1, centers1, normals1);
  Christian::conicBackprojection(conic_locus2, radius2, centers2, normals2);
  Christian::conicBackprojection(conic_locus3, radius3, centers3, normals3);

  uint planes12_idx = chooseSupportingPlanes(angle12, normals1, normals2);
  uint planes23_idx = chooseSupportingPlanes(angle23, normals2, normals3);
  uint planes31_idx = chooseSupportingPlanes(angle31, normals3, normals1);

  uint idx1_1 = planes12_idx / 2;
  uint idx2_1 = planes12_idx % 2;
  uint idx2_2 = planes23_idx / 2;
  uint idx3_1 = planes23_idx % 2;
  uint idx3_2 = planes31_idx / 2;
  uint idx1_2 = planes31_idx % 2;
  assert(idx1_1 == idx1_2);
  assert(idx2_1 == idx2_2);
  assert(idx3_1 == idx3_2);
  Eigen::Vector3d normal1 = normals1.at(idx1_1);
  Eigen::Vector3d normal2 = normals2.at(idx2_1);
  Eigen::Vector3d normal3 = normals3.at(idx3_1);

  double angle_calc12 = getAngleBetweenVectors(normal1, normal2);
  double angle_calc23 = getAngleBetweenVectors(normal2, normal3);
  double angle_calc31 = getAngleBetweenVectors(normal3, normal1);
  assert(std::abs(angle12 - angle_calc12) < 1e-8);
  assert(std::abs(angle23 - angle_calc23) < 1e-8);
  assert(std::abs(angle31 - angle_calc31) < 1e-8);
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

void reprojectionsToPlanes( const std::vector<Eigen::Vector3d>& centers,
                            const std::vector<Eigen::Vector3d>& normals,
                            std::vector<Eigen::Hyperplane<double, 3> >& planes) {
  std::vector<Eigen::Vector3d>::const_iterator the_chosen;
  for(the_chosen = normals.begin(); the_chosen != normals.end(); the_chosen++) {
    int index = getIndex(normals.begin(), the_chosen);
    Eigen::Vector3d center = centers.at(index);
    Eigen::Hyperplane<double, 3> plane(*the_chosen, center);
    planes.push_back(plane);
  }
}

void calculateHomography( const std::vector<Eigen::Hyperplane<double, 3> >& planes_world,
                          const std::vector<Eigen::Hyperplane<double, 3> >& planes_cam,
                          Eigen::Quaterniond& attitude, Eigen::Vector3d& position) {
  uint PLANE_DIM = 4;
  uint n_planes = planes_world.size();
  assert(n_planes == planes_cam.size());
  Eigen::VectorXd lhs(PLANE_DIM*n_planes);
  Eigen::MatrixXd rhs(PLANE_DIM*n_planes, PLANE_DIM*(PLANE_DIM-1));
  for(size_t i = 0; i < n_planes; i++) {
    lhs.block<4,1>(4*i,0) = planes_cam.at(i).coeffs();
    Eigen::Hyperplane<double, 3> world_plane = planes_world.at(i);
    lhs(4*i+3) -= world_plane.offset();
    
    double a = world_plane.normal()(0);
    Eigen::Matrix4d diag_a = a*Eigen::Matrix4d::Identity();
    double b = world_plane.normal()(1);
    Eigen::Matrix4d diag_b = b*Eigen::Matrix4d::Identity();
    double c = world_plane.normal()(2);
    Eigen::Matrix4d diag_c = c*Eigen::Matrix4d::Identity();
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*0) = diag_a;
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*1) = diag_b;
    rhs.block<4,4>(PLANE_DIM*i, PLANE_DIM*2) = diag_c;
  }

  Eigen::BDCSVD<Eigen::MatrixXd> SVD(rhs, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd solu = SVD.solve(lhs);
  Eigen::MatrixXd P_est = solu.reshaped(4,3).transpose();
  Eigen::Matrix3d dcm = P_est.topLeftCorner(3,3);
  attitude = Eigen::Quaterniond(dcm).inverse();
  position = P_est.topRightCorner(3,1);
}

void solve_navigation_problem(const std::vector<Quadric>& quadrics,
                              const std::vector<Eigen::Matrix3d>& locii,
                              Eigen::Quaterniond& attitude, Eigen::Vector3d& position) {
  std::vector<Eigen::Hyperplane<double, 3> > planes_world;
  for(std::vector<Quadric>::const_iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    planes_world.push_back((*it).getPlane());
  }
  std::vector<Eigen::Vector3d> centers;
  std::vector<Eigen::Vector3d> normals;
  reprojectLociiToQuadrics(quadrics, locii, centers, normals);
  std::vector<Eigen::Hyperplane<double, 3> > planes_cam;
  reprojectionsToPlanes(centers, normals, planes_cam);
  calculateHomography(planes_world, planes_cam, attitude, position);
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

  cv::Mat image = cam->getBlankCameraImage();
  cv::Mat outImg;
  double scaling = 0.4;
  // Showing image inside a window 
  viz::drawEllipses(image, *cam, quadrics, viz::CV_colors);
  // viz::drawEllipses(image, conics, viz::CV_colors);
  cv::resize(image, outImg, cv::Size(), scaling, scaling);
  cv::imshow("Projecting Moon to Camera", outImg); 
  cv::waitKey(0); 

  std::vector<Eigen::Matrix3d> locii;
  for(std::vector<Quadric>::iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    locii.push_back((*it).projectToPlaneLocus(extrinsic));
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
