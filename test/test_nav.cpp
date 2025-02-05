#include "navigation.h"
#include "camera.h"
#include "quadrics.h"
#include "conics.h"
#include "visuals.h"

#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Core>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/calib3d/calib3d.hpp"
#include <opencv2/core/eigen.hpp>
#include <random>
#include <fstream>
#include <string>

class NavigationTest : public testing::Test {
protected:
  NavigationTest() {
    latitude = -30, longitude = 0, radius = 200;
    id = "defaultQuadric"; // .as<std::string>()???
    quadric_default = new Quadric(latitude, longitude, radius, id);
    dx = 1000, dy = 1000, skew = 0, im_height = 2048, im_width = 2592;
    up = (im_width+1)/2.;
    vp = (im_height+1)/2.;
    image_size.width  = im_width;
    image_size.height = im_height;
    cam = new Camera(dx, dy, up, vp, skew, image_size, quat, position);

    double cassini_f = 165714.2857142857;
    Eigen::Vector3d cassini_pos;
    cassini_pos << 360276.0,  90684.3, -44379.2;
    Eigen::Quaterniond cassini_att(-0.58826, 0.524272, -0.407365, -0.461674);
    c_dx = cassini_f, c_dy = cassini_f, c_skew = 0, c_im_height = 1024, c_im_width = 1024;
    c_up = (c_im_width+1)/2.;
    c_vp = (c_im_height+1)/2.;
    c_image_size.width  = c_im_width;
    c_image_size.height = c_im_height;
    cassini_cam = new Camera(c_dx, c_dy, c_up, c_vp, c_skew, c_image_size, cassini_att, cassini_pos);
    // cassini_cam->setAttitude(cassini_att);
    // cassini_cam->setPosition(cassini_pos);
    cv_cam = new cv::viz::Camera(dx, dy, up, vp, image_size);
  }
  ~NavigationTest() override {
    delete quadric_default;
    delete cam;
    delete cv_cam;
    delete cassini_cam;
  }
  public:
  double latitude, longitude, radius;
  std::string id;
  Quadric* quadric_default;
    double dx, dy, skew, up, vp;
    uint im_height, im_width;
    double c_dx, c_dy, c_skew, c_up, c_vp;
    uint c_im_height, c_im_width;
    cv::Size2i image_size;
    cv::Size2i c_image_size;
    Eigen::Quaterniond quat = Eigen::Quaterniond::Identity();
    Eigen::Vector3d position{1e4, 0, 0};

    Camera* cam;
    cv::viz::Camera* cv_cam;
    Camera* cassini_cam;

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
    // std::cout << Conic((*it).projectToPlaneLocus(extrinsic)) << std::endl;
  }

  Eigen::Quaterniond attitude;
  Eigen::Vector3d position;
  solve_navigation_problem(quadrics, locii, attitude, position);
  ASSERT_TRUE(cam->getAttitude().isApprox(attitude, 1e-3));
  ASSERT_TRUE(cam->getPosition().isApprox(position, 1e-3));
}

// TEST_F(NavigationTest, KanataniEllipseCheck) {
//   // double semimajor = 100, semiminor = 40, xc = 40, yc = 120, phi = deg2rad(10);
//   // Conic conic(semimajor, semiminor, xc, yc, phi);
// }

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

TEST_F(NavigationTest, QuadricPointsWNoise) {
  const int n_pts = 15;

  cam->resetCameraState();
  cam->moveX(2.5e3);
  cam->moveY(-1e2);
  cam->moveZ(3e2);
  Eigen::Vector3d look_here = 1e2*Eigen::Vector3d::Random();
  Eigen::Vector3d up_vector = Eigen::Vector3d::Random();
  cam->pointTo(look_here, up_vector);
  // std::cout << *cam << std::endl;
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();

  Quadric q1( 15, -20, 100, "crater1");
  Quadric q2(-20,  15, 200, "crater2");
  Quadric q3( 10,  05, 150, "crater3");
  Quadric q4(-10,  10, 130, "crater4");
  Quadric q5(-30, -20, 250, "crater5");
  std::vector<Quadric> quadrics = {q1, q2, q3, q4, q5};

  double mean = 0, st_dev = 1.5;
  std::vector<std::vector<Eigen::Vector2d> > noisy_pts;
  noisy_pts.reserve(quadrics.size());
  std::vector<Eigen::Matrix3d> locii_noisy, locii_points, locii_truth;
  std::vector<Conic> conics, noisy_conics;
  Eigen::Matrix3d T_e2m = cam->getAttitudeMatrix();

  for(std::vector<Quadric>::const_iterator it = quadrics.begin(); it != quadrics.end(); it++) {
    Eigen::Vector3d normal_cam = T_e2m * (*it).getNormal();
    // std::cout << "Normal dot: " << normal_cam.normalized() << std::endl;
    if(Eigen::Vector3d::UnitZ().dot(normal_cam.normalized()) > 0) {
      std::cout << "Normal is pointing in the wrong direction: " << normal_cam.normalized().transpose() << std::endl;
      continue;
    }
    std::vector<Eigen::Vector3d> pts_world;
    pts_world.reserve(n_pts);
    (*it).getRimPoints(n_pts, pts_world);
    std::vector<Eigen::Vector2d> pts_pxl;
    cam->world2Pixel(pts_world, pts_pxl);
    Conic fromQuadric = (*it).projectToImage(cam->getProjectionMatrix());
    Eigen::Matrix3d plane_points = cam->getImagePlaneLocus(fromQuadric.getLocus());
    locii_points.push_back(plane_points);

    std::random_device rd;
    std::mt19937 gen(rd());
    const double offset_std = 0;
    std::normal_distribution<double> dist(0, offset_std);
    addNoise(mean + dist(gen) - offset_std/2.0, st_dev, pts_pxl);
    Conic noisy(pts_pxl);
    noisy_conics.push_back(noisy);
    Eigen::Matrix3d locus_noisy = cam->getImagePlaneLocus(noisy.getLocus());
    locii_noisy.push_back(locus_noisy);
    noisy_pts.push_back(pts_pxl);
    conics.push_back(fromQuadric);

    Eigen::Matrix3d plane_truth = (*it).projectToPlaneLocus(extrinsic);
    locii_truth.push_back(plane_truth);
  }

  Eigen::Quaterniond attitude;
  Eigen::Vector3d position;
  solve_navigation_problem(quadrics, locii_noisy, attitude, position);
  // std::cout << "ESTIMATION:"
  //           << "\n\tLocation: " << position.transpose()
  //           << "\n\tAttitude: " << attitude << std::endl;
  // std::cout << "ACTUAL CAMERA POSE:\n" << *cam << std::endl;
  // std::cout << "Attitude Error: " << rad2deg(attitudeError(attitude, cam->getAttitude())) << " deg" << std::endl;
  // std::cout << "Position Error: " << (position - cam->getPosition()).norm() << std::endl;
  EXPECT_TRUE(cam->getAttitude().isApprox(attitude, 3));
  EXPECT_TRUE(cam->getPosition().isApprox(position, 1e2));

  // // Visuals
  // viz::drawNoisyPoints(*cam, conics, noisy_pts);
}

void drawLabels(cv::Mat& image, const Camera& camera, const std::vector<Quadric>& craters) {
  Eigen::MatrixXd proj_mtx = camera.getProjectionMatrix();

  int font = cv::FONT_HERSHEY_SIMPLEX;
  float font_scale = 0.2; // scale factor from base size
  int thickness = 1; //in pixels

  for(const Quadric& crater : craters) {
    Eigen::Matrix3d conic_envelope = proj_mtx * crater.getEnvelope() * proj_mtx.transpose();
    Eigen::Matrix3d locus = adjugate(conic_envelope);
    Conic conic(locus);
    cv::Size2i img_size(image.cols, image.rows);
    cv::Point2d ellipse_center;
    conic.getCenter(ellipse_center);
    if(!isInImage(ellipse_center, img_size)) {
      std::cerr << __func__ << ": Ellipse is not in the image: " << ellipse_center << std::endl;
      return;
    }
    cv::Scalar rect_color = cv::viz::Color::orange();
    cv::Scalar text_color = cv::viz::Color::azure();
    cv::Point2d ll_offset(-5, -10);
    cv::Point2d ur_offset( 5*crater.getID().length(),  30);
    cv::rectangle (image, ellipse_center+ll_offset, ellipse_center+ur_offset, rect_color, cv::FILLED, cv::LINE_8, 0);
    cv::putText(image, crater.getID(), ellipse_center, font, font_scale,
                text_color, thickness, cv::LINE_AA, false);
  }
}

// Using https://planetarynames.wr.usgs.gov/Feature/___ to get locations of craters
// E.g., for Beaumont, use https://planetarynames.wr.usgs.gov/Feature/653
std::vector<std::array<double, 2> > theophilus = {{627, 473}, {633, 468}, {650, 465}, {661, 472}, {665, 484}, {663, 493}, {656, 504}, {644, 508}, {632, 505}, {625, 499}};
std::vector<std::array<double, 2> > catharina = {{665, 388}, {672, 383}, {684, 378}, {697, 382}, {704, 392}, {705, 404}, {697, 414}, {668, 413}, {663, 407}};
std::vector<std::array<double, 2> > beaumont = {{609, 393}, {613, 388}, {622, 388}, {628, 394}, {627, 401}, {622, 406}, {616, 407}, {608, 404}};
std::vector<std::array<double, 2> > madler = {{594, 489}, {596, 486}, {599, 485}, {604, 488}, {604, 492}, {601, 495}, {596, 495}, {594, 493}};
std::vector<std::array<double, 2> > tacitus = {{739, 417}, {746, 417}, {750, 423}, {747, 429}, {742, 433}, {737, 432}, {733, 428}};
std::vector<std::array<double, 2> > alfraganus = {{735, 567}, {739, 566}, {743, 568}, {742, 572}, {739, 574}, {736, 574}, {734, 572}, {733, 570}};
std::vector<std::array<double, 2> > polybius = {{655, 336}, {658, 333}, {665, 334}, {670, 338}, {669, 344}, {666, 348}, {659, 349}, {655, 347}, {653, 344}};
std::vector<std::array<double, 2> > kant = {{719, 494}, {723, 493}, {727, 495}, {730, 499}, {728, 504}, {724, 507}, {717, 505}, {715, 502}};
std::vector<std::array<double, 2> > bohnemberger = {{475, 412}, {479, 409}, {483, 408}, {486, 411}, {486, 417}, {483, 420}, {477, 421}, {474, 418}};
std::vector<std::array<double, 2> > plinius = {{671, 836}, {676, 833}, {683, 833}, {688, 839}, {687, 844}, {684, 848}, {676, 850}, {672, 848}, {670, 844}};
std::vector<std::array<double, 2> > macrobius = {{412, 891}, {418, 888}, {427, 891}, {433, 899}, {434, 907}, {429, 912}, {422, 911}, {413, 906}, {409, 901}};
std::vector<std::array<double, 2> > taruntius = {{382, 691}, {391, 688}, {398, 692}, {401, 700}, {399, 709}, {391, 712}, {384, 710}, {379, 703}};
std::vector<std::array<double, 2> > cook = {{378, 385}, {384, 381}, {389, 381}, {393, 389}, {391, 395}, {386, 400}, {378, 398}, {375, 394}};

Quadric Theophilus(-11.45, 26.28, 98.59/2., "Theophilus");
Quadric Catharina(-17.98, 23.55, 98.77/2., "Catharina");
Quadric Beaumont(-18.08, 28.82, 50.69/2., "Beaumont");
Quadric Madler(-11.04, 29.76, 27.58/2., "Madler");
Quadric Tacitus(-16.20, 18.95, 39.81/2., "Tacitus");
Quadric Alfraganus(-5.42, 18.97, 20.52/2., "Alfraganus");
Quadric Polybius(-22.46, 25.63, 40.81/2., "Polybius");
Quadric Kant(-10.62, 20.20, 30.85/2., "Kant");
Quadric Bohnemberger(-16.24, 40.06, 31.74/2., "Plinius");
Quadric Plinius(15.36, 23.61, 41.31/2., "Macrobius");
Quadric Macrobius(21.26, 45.97, 62.79/2., "Polybius");
Quadric Taruntius(5.50, 46.54, 57.32/2., "Taruntius");
Quadric Cook(-17.50, 48.81, 45.16/2., "Cook");

const std::vector<std::vector<std::array<double, 2> > > craters = {
  theophilus,
  catharina,
  beaumont,
  madler,
  tacitus,
  alfraganus,
  polybius,
  kant,
  bohnemberger,
  plinius,
  macrobius,
  taruntius,
  cook
};
const std::vector<Quadric> quadrics = {
  Theophilus,
  Catharina,
  Beaumont,
  Madler,
  Tacitus,
  Alfraganus,
  Polybius,
  Kant,
  Bohnemberger,
  Plinius,
  Macrobius,
  Taruntius,
  Cook
};

std::vector<cv::Point2d> Generate2DPoints() {
  std::vector<cv::Point2d> points;
  for(const std::vector<std::array<double, 2> >& crater : craters) {

    Eigen::MatrixXd crater_pts = toEigenArray(crater);
    std::array<double, IMPLICIT_PARAM> impl = ellipseFitLstSq(crater_pts);
    Conic conic(impl);
    cv::Point2d point;
    conic.getCenter(point);
    points.push_back(point);
    // std::cout << "Conic center: " << point << std::endl;
  }
  return points;
}

std::vector<cv::Point3d> Generate3DPoints() {
  std::vector<cv::Point3d> points;
  for(const Quadric& crater : quadrics) {
    Eigen::Vector3d location = crater.getLocation();
    points.push_back(cv::Point3d(location(0), location(1), location(2)));
  }
  return points;
}

TEST_F(NavigationTest, Triangulation) {
  double f = 165714.2857142857;
  double dx = f, dy = f, skew = 0;
  uint im_height, im_width;
  im_height = 1024, im_width = 1024;
  cv::Size2i image_size(im_height, im_width);
  up = (double(im_width )+1)/2;
  vp = (double(im_height)+1)/2;
  image_size.width  = im_width;
  image_size.height = im_height;
  Camera camera(dx, dy, up, vp, skew, image_size);
  // std::cout << camera << std::endl;
  // Read points
  std::vector<cv::Point2d> imagePoints = Generate2DPoints();
  std::vector<cv::Point3d> objectPoints = Generate3DPoints();

  // std::cout << "There are " << imagePoints.size() << " imagePoints and " << objectPoints.size() << " objectPoints." << std::endl;
  cv::Mat cameraMatrix(3,3,cv::DataType<double>::type);
  // cv::setIdentity(cameraMatrix);
  cameraMatrix = camera.getCvIntrinsicMatrix();

  // std::cout << "Initial cameraMatrix:\n" << cameraMatrix << std::endl;

  cv::Mat distCoeffs(4,1,cv::DataType<double>::type);
  distCoeffs.at<double>(0) = 0;
  distCoeffs.at<double>(1) = 0;
  distCoeffs.at<double>(2) = 0;
  distCoeffs.at<double>(3) = 0;

  cv::Mat rvec(3,1,cv::DataType<double>::type);
  cv::Mat tvec(3,1,cv::DataType<double>::type);

  cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec, tvec);

  // std::cout << "rvec: " << rvec << std::endl;
  cv::Mat rotMtx(3,3,cv::DataType<double>::type);
  cv::Rodrigues(rvec, rotMtx);
  // std::cout << "Rot matrix:\n" << rotMtx << std::endl;
  Eigen::Matrix3d rMtx;
  cv::cv2eigen(rotMtx, rMtx);
  Eigen::Vector3d tVec;
  cv::cv2eigen(tvec, tVec);
  Eigen::Quaterniond q(rMtx.transpose());
  std::cout << "Camera wrt Moon in Camera Frame:\n " << tvec << std::endl;
  std::cout << "Camera wrt Moon:\n" << q * tVec << std::endl;
  std::cout << "Distance: " << (q * tVec).norm() << std::endl;
  std::cout << "Quat: " << q << std::endl;
  std::cout << "Matrix: " << q.toRotationMatrix() << std::endl;

  std::vector<cv::Point2d> projectedPoints;
  cv::projectPoints(objectPoints, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);

  // for(unsigned int i = 0; i < projectedPoints.size(); ++i)
  //   {
  //   std::cout << "Image point: " << imagePoints[i] << " Projected to " << projectedPoints[i] << std::endl;
  //   }
}


TEST_F(NavigationTest, RealImage) {

  std::vector<Conic> conics;
  std::vector<Eigen::Matrix3d> locii_noisy;
  std::vector<std::vector<std::array<double, 2> > >::const_iterator crater;
  for(crater = craters.begin(); crater != craters.end(); crater++) {
    Eigen::MatrixXd crater_pts = toEigenArray(*crater);
    std::array<double, IMPLICIT_PARAM> impl = ellipseFitLstSq(crater_pts);
    conics.push_back(Conic(impl));
    locii_noisy.push_back(Conic(impl).getLocus());
  }

  Eigen::Quaterniond attitude;
  Eigen::Vector3d position;
  solve_navigation_problem(quadrics, locii_noisy, attitude, position);
  std::cout << "ESTIMATION:"
            << "\n\tLocation: " << position.transpose()
            << "\n\tAttitude: " << attitude << std::endl;

  double f = 165714.2857142857;
  double dx = f, dy = f, skew = 0;
  // double dx = 7291.6666, dy = 7291.6666, skew = 0;
  uint im_height, im_width;
  im_height = 1024, im_width = 1024;
  cv::Size2i image_size(im_height, im_width);
  up = (double(im_width )+1)/2;
  vp = (double(im_height)+1)/2;
  image_size.width  = im_width;
  image_size.height = im_height;
  // cassini_cam->resetCameraState();
    Eigen::Vector3d cassini_pos;
    cassini_pos << 360276.0,  90684.3, -44379.2;
    Eigen::Quaterniond cassini_att(0.58826, 0.524272, -0.407365, -0.461674);
  Camera camera(dx, dy, up, vp, skew, image_size);
  camera.moveX(360276.0);
  camera.moveY( 90684.3);
  camera.moveZ(-44379.2);
  Eigen::Matrix3d att;
  att << 0.241822,    0.116029,   -0.963358,
        -0.970309,   0.0239921,   -0.240677,
        -0.00481262,    0.992956,    0.118386;

  // att <<  0.2342060407733891, -0.9721684019135322, -0.006010722597369456,
  //         0.09936144043346457, 0.01778604106378229, 0.9948924368484588,
  //        -0.9670960834478861, -0.2336070526849872, 0.1007616510166092;
  camera.setAttitude(att.transpose());
  std::cout << camera << std::endl;
  std::cout << "Cassini: " << *cassini_cam << std::endl;

  std::string filename = "/home/ndvr/data/pds/cassini/image.png";
  cv::Mat image_gray = cv::imread(filename, 0);
  cv::Mat image;
  cv::cvtColor(image_gray, image, 1);
  EXPECT_EQ(image.rows, c_im_height);
  EXPECT_EQ(image.cols, c_im_width);
  // Add the moon
  // viz::drawMoon(image, camera);
  viz::draw3dAxes(image, *cassini_cam);
  cv::Point center(774, 554);
  cv::drawMarker(image, center, cv::viz::Color::yellow());

  // viz::drawEllipses(image, conics, viz::CV_colors);

  viz::drawCraters(image, *cassini_cam, quadrics);
  Eigen::MatrixXd proj = cassini_cam->getProjectionMatrix();
  for(uint i = 0; i < conics.size(); i++) {
    Quadric crater = quadrics.at(i);
    Conic conic = conics.at(i);
    Conic projected = crater.projectToImage(proj);

    cv::Point2d ellipse_center, projected_center;
    conic.getCenter(ellipse_center);
    projected.getCenter(projected_center);
    cv::arrowedLine(image, projected_center, ellipse_center, cv::viz::Color::red());
  }
  // drawLabels(image, camera, quadrics);
  cv::imshow("Cassini Image", image);
  cv::waitKey(0);
}

TEST_F(NavigationTest, Eccentricity) {
  double smajor = 6378137, sminor = 6356752, xc = 30, yc = 20, phi = 0;
  Conic ellipse(smajor, sminor, xc, yc, phi);
  ASSERT_NEAR(ellipse.getEccentricity(), 0.0818198, 1e-6);
}

TEST_F(NavigationTest, AngleOfEllipse) {
  // Create quadric
  Quadric quad(0,0,100,"Angle");
  // Place camera
  cam->setLocation(Eigen::Vector3d(4e3, 0e3, 3e3));
  // Point camera to quadric
  cam->pointTo(quad.getLocation(), Eigen::Vector3d::UnitZ());
  cam->moveY(-2e3);

  // Get distance to quadric
  Eigen::Matrix3d conic = quad.projectToImageLocus(cam->getProjectionMatrix());
  // std::cout << Conic(conic) << std::endl;
  std::array<Eigen::Vector3d, 2> centers, normals;

  Eigen::Matrix3d att = cam->getAttitudeMatrix();
  Eigen::Vector3d centerWrtCamWorld = quad.getLocation() - cam->getLocation();
  Eigen::Vector3d centerWrtCam = att * centerWrtCamWorld;
  Eigen::Vector3d normal_test = att * quad.getNormal();
  double dist2center = (centerWrtCamWorld).norm();
  std::cout << "---\n---Truth Data---\n";
  std::cout << "---DISTANCE TO QUADRIC CENTER: " << dist2center << std::endl;
  std::cout << "---CAMERA TO QUADRIC CENTER: " << centerWrtCam.transpose() << std::endl;
  std::cout << "---QUADRIC NORMAL WRT CAMERA: " << normal_test.transpose() << std::endl;

  Eigen::Vector3d center, normal;
  // double dist;
  Eigen::MatrixXd extrinsic = cam->getExtrinsicMatrix();
  Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extrinsic);
  Christian::conicBackprojection(conic_locus, quad.getRadius(), centers, normals);

  uint correct_index = 0;
  center = centers.at(correct_index);
  normal = normals.at(correct_index);
  std::cout << "+++\n+++Calculations+++\n";
  std::cout << "Center (cam): " << center.transpose() << std::endl;
  std::cout << "Normal (cam): " << normal.transpose() << std::endl;
  Eigen::Vector3d q_normal = Eigen::Vector3d::UnitZ();
  std::cout << "Angle 1: " << rad2deg(std::acos(q_normal.dot(-normal))) << std::endl;
  std::cout << "Correct angle: " << rad2deg(std::acos(q_normal.dot(-normal_test))) << std::endl;

  EXPECT_TRUE(centerWrtCam.isApprox(center, 1e-3));
  EXPECT_TRUE( normal_test.isApprox(normal, 1e-3));

  // Determine angle of supporting plane wrt camera plane
  double eccentricity = Conic(conic).getEccentricity();
  Conic plane_conic(conic_locus);
  double theta = std::asin(eccentricity / std::sqrt(1+std::pow(quad.getRadius(),2)));
  std::cout << "eccentricity: " << eccentricity << std::endl;
  std::cout << "Angle calc: " << theta << " (" << rad2deg(theta) << ")" << std::endl;

  double target_theta_deg = 52.9414;
  double target_theta = deg2rad(target_theta_deg);
  double r = std::sqrt(std::pow(eccentricity, 2)/std::pow(std::sin(target_theta), 2) - 1);
  std::cout << "Expected radius: " << r << std::endl;

  // Eigen::Vector3d k_normal;
  // double dist;
  // Kanatani::conicBackprojection(conic_locus, quad.getRadius(), k_normal, dist);
  // std::cout << "Distance: " << dist << std::endl;
  // std::cout << "Normal: " << k_normal.transpose() << std::endl;
  // EXPECT_TRUE(k_normal.isApprox(normal,  1e-3) || k_normal.isApprox(-normal,  1e-3));

  // double lambda1, lambda2, lambda3;
  // Eigen::Vector3d u1, u2, u3;
  // // We need to switch the 1 & 2 indices for Kanatani's method
  // Preprocessing::getEigenParts(conic_locus, lambda2, lambda1, lambda3, u2, u1, u3);

  // double scalar1 = std::sqrt((lambda2 - lambda1)/(lambda2 - lambda3));
  // double scalar2 = std::sqrt((lambda1 - lambda3)/(lambda2 - lambda3));
  // Eigen::Vector3d normal2 = -(scalar1*u2 + scalar2*u3);
  // std::cout << "Alternative normal: " << normal2.transpose() << std::endl;

}

TEST_F(NavigationTest, AngleToTheta) {
  // Create quadric
  Quadric quad(10,15,100,"Angle");
  Eigen::Matrix4d quadric_envelope = quad.getEnvelope();
  // Place camera
  double offset = 0.5e3;
  uint n_steps = 240;
  double d_angle = M_PI_2 / (double)n_steps;
  uint scalar = 300;
  cv::Mat image(scalar, 360,
                CV_8UC3,
                cv::Scalar(50, 50, 50));
  const double center_dist = quad.getCenter().norm();
  cv::Mat locations(offset + center_dist, offset + center_dist,
                    CV_8UC3,
                    cv::Scalar(50, 50, 50));
  for(uint i = 0; i <= n_steps; i++) {
    double angle = i * d_angle;
    cam->setLocation(Eigen::Vector3d(center_dist + offset*std::sin(angle), offset*std::cos(angle), 0));
    cam->pointTo(quad.getLocation(), Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d conic_locus;
    Eigen::MatrixXd proj = cam->getProjectionMatrix();
    try {
      conic_locus = cam->projectQuadricToLocus(quadric_envelope);
    }
    catch (const std::exception& e) {
      std::cerr << "Error in conic locus\n";
      continue;
    }
    Conic conic(conic_locus);
    // std::cout << conic << std::endl;
    // std::cout << 1296.5 - conic.getCenterX() << std::endl;
    double smajor = conic.getSemiMajorAxis();
    double sminor = conic.getSemiMinorAxis();

    double center_off = 1296.5 - conic.getCenterX();
    double axis_ratio = sminor / smajor;
    cv::Point center(4*rad2deg(angle), 0.5*scalar+10*(center_off));
    cv::Point sine(4*rad2deg(angle), 0.5*scalar+10*(-1+std::abs(std::sin(angle))));
    cv::Point axes(4*rad2deg(angle), scalar*axis_ratio);
    // cv::drawMarker(image, center, cv::viz::Color::red(), cv::MARKER_STAR, 5);
    cv::drawMarker(image, axes, cv::viz::Color::green(), cv::MARKER_CROSS, 5);
    cv::drawMarker(image, sine, cv::viz::Color::azure(), cv::MARKER_CROSS, 5);
    cv::Point location(center_dist + offset*std::sin(angle), offset*std::cos(angle));
    cv::drawMarker(locations, location, cv::viz::Color::azure(), cv::MARKER_CROSS, 5);

    // std::cout << "Angle: " << rad2deg(angle);// << std::endl;//<< " -> " << conic << std::endl;
    // // std::cout << "\tRelative Position: " << (quad.getLocation() - cam->getLocation()).transpose() << std::endl;
    // // std::cout << "\tCenter offset: " << center_off << std::endl;
    // std::cout << "\tAxis ratio: " << axis_ratio << std::endl;// << std::endl;
  }
  // // Show the image with detected circles
  // cv::imshow("Axes and Offsets", image);
  // cv::imshow("Locations in x/y", locations);
  // cv::waitKey(0);

}

double getAngleFromBoresight(const Camera& camera, const Eigen::Vector2d& point) {
  // const Eigen::Matrix3d Kinv = camera.getInverseIntrinsicMatrix();
  return 1.0;
}

// double getPhi(const Eigen::Matrix3d& conic_locus) {
//   double lambda1, lambda2, lambda3;
//   Eigen::Vector3d u1, u2, u3;
//   Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);
//   const double abs_l1 = std::abs(lambda1);
//   const double abs_l2 = std::abs(lambda2);
//   const double abs_l3 = std::abs(lambda3);
//   const double radical = (abs_l1 - abs_l2)/(abs_l2 + abs_l3);
//   const double phi = std::atan(std::sqrt(radical));
//   // const double phi = std::atan(abs_l3/abs_l1*std::sqrt(radical));
//   return phi;
// }

// double getPhi2(const Eigen::Matrix3d& conic_locus) {
//   double lambda1, lambda2, lambda3;
//   Eigen::Vector3d u1, u2, u3;
//   Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);

//   const double kx2 = 1/std::abs(lambda1);
//   const double ky2 = 1/std::abs(lambda2);
//   const double kz2 = 1/std::abs(lambda3);
//   const double alpha = (ky2 - kx2) / (ky2 + kz2);
//   const double aangle = std::atan(-std::sqrt(alpha)/(-kx2/kz2));
//   std::cout << "Lambdas: " << lambda1 << ", " << lambda2 << ", " << lambda3 << "\n\n";
//   return aangle;
// }

double getRho(const Eigen::Matrix3d& conic_locus) {
  double lambda1, lambda2, lambda3;
  Eigen::Vector3d u1, u2, u3;
  Preprocessing::getEigenParts(conic_locus, lambda1, lambda2, lambda3, u1, u2, u3);

  const double kx2 = 1/std::abs(lambda1);
  const double ky2 = 1/std::abs(lambda2);
  const double kz2 = 1/std::abs(lambda3);
  const double alpha = (ky2 - kx2) / (ky2 + kz2);
  const double kxz2 = kx2/kz2;
  const double rho = std::sqrt(alpha+kxz2)/(1-alpha);
  return rho;
}

// void appendToYaml(const std::string& filePath, const std::string& key, const YAML::Node& value) {
//     YAML::Node doc = YAML::LoadFile(filePath); // Load the existing YAML file
//     doc[key] = value; // Add the new entry
//     std::ofstream fout(filePath); // Open the file for writing
//     fout << doc; // Write the modified YAML node back to the file
//     fout.close(); // Close the file
// }

void putYamlElement(const std::string filename,
                    const uint indent,
                    const std::string& label,
                    const std::string& value) {
  std::string indentation(indent, ' ');
  std::string text = indentation + label + ": " + value + "\n";
  std::ofstream file(filename, std::ios_base::app);
  if(file.is_open()) {
    file << text;
    file.close();
  }
  else {
    std::cerr << "File: `" << filename << "` not open!\n";
  }
}

void putVectorIntoYaml( const std::string filename, 
                        const uint indent,
                        const std::string& label,
                        const Eigen::Vector3d& vec) {
  putYamlElement(filename, indent, label, "");
  putYamlElement(filename, indent+2, "x", std::to_string(vec(0)));
  putYamlElement(filename, indent+2, "y", std::to_string(vec(1)));
  putYamlElement(filename, indent+2, "z", std::to_string(vec(2)));
}

void putCraterIntoYaml( const std::string filename, 
                        const Eigen::Vector3d& location, 
                        const Eigen::Vector3d& normal, 
                        const double radius,
                        const Eigen::Vector3d& ext_circle, 
                        const double ext_radius,
                        const std::string& label,
                        const uint indent) {
  std::string str_name = "- name";
  std::string str_loc = "location";
  std::string str_norm = "normal";
  std::string str_rad = "radius";
  std::string str_circ = "circle";
  std::string str_rcirc = "ext_radius";
  putYamlElement(filename, indent, str_name, label);
  putVectorIntoYaml(filename, indent+2, str_loc, location);
  putVectorIntoYaml(filename, indent+2, str_norm, normal);
  putYamlElement(filename, indent+2, str_rad, std::to_string(radius));
  putVectorIntoYaml(filename, indent+2, str_circ, ext_circle);
  putYamlElement(filename, indent+2, str_rcirc, std::to_string(ext_radius));
}

void putCameraIntoYaml( const std::string filename, 
                        const Eigen::Vector3d& location) {
  putYamlElement(filename, 0, "camera", "");
  putVectorIntoYaml(filename, 2, "location", location);
}

TEST_F(NavigationTest, CalcAngle) {
  // Clear contents of output.yaml
  std::string filename = "../data/circles.yaml";
  std::ofstream file(filename, std::ios::trunc);
  file.close();
  // Create quadric
  Quadric quad1(30,-5,100,"Crater1");
  Quadric quad2(5,40,200,"Crater2");
  Quadric quad3(-15,5,150,"Crater3");
  Quadric quad4(15,15,230,"Crater4");
  Quadric quad5(-35,-35,250,"Crater5");

  const double offset = 0.5e3;
  Eigen::Vector3d crater_center = quad1.getLocation();
  Eigen::Vector3d look_at(0,0,0);
  // Eigen::Vector3d look_at = (quad1.getLocation() + quad2.getLocation())/2;
  Eigen::Vector3d location(crater_center(0)+2*offset, 0, -1*offset);
  cam->setPosition(location);
  cam->pointTo(look_at, Eigen::Vector3d::UnitZ());

  std::vector<Quadric> quads = {quad1, quad2, quad3, quad4, quad5};
  // const double eff_foc = 1e0;
  const Eigen::MatrixXd extr_mtx = cam->getExtrinsicMatrix();
  const Eigen::MatrixXd proj_mtx = cam->getProjectionMatrix();
  // const Eigen::Quaterniond att = cam->getAttitude();
  const Eigen::Vector3d cam_pos = cam->getLocation();
  putCameraIntoYaml(filename, cam_pos);
  putYamlElement(filename, 0, "craters", "");

  uint count = 0;
  for(auto& quad : quads) {
    count++;
    // Eigen::Matrix3d conic_locus = cam->projectQuadricToLocus(quad.getEnvelope());
    Eigen::Matrix3d conic_locus = quad.projectToPlaneLocus(extr_mtx);
    
    // std::cout << quad << std::endl;
    Conic conic(conic_locus);
    const double q_radius = quad.getRadius();
    const Eigen::Vector3d q_loc = quad.getLocation();
    const Eigen::Vector3d q_nor = quad.getNormal();

    // const double Rmax = conic.getSemiMajorAxis();
    // const double Rmin = conic.getSemiMinorAxis();
    // const double alpha = calcAngleFromEllipseAxes(Rmax, Rmin, eff_foc);
    // // const double h = q_radius * (eff_foc/Rmin*std::cos(alpha));
    const double angle = calcAngleFromEllipseAxes(conic_locus);

    // Eigen::Vector3d normal_test = att * q_nor;
    Eigen::Vector3d r_q_cam_world = cam->getPointWrtCameraWorld(q_loc);
    Eigen::Vector3d u_q_cam_world = r_q_cam_world.normalized();
    const double expected_angle = std::acos(q_nor.dot(-u_q_cam_world));
    // const double expected_angle = std::acos(normal_test.dot(-Eigen::Vector3d::UnitZ()));

    const double distance = Shiu::getDistanceFromConicLocus(conic_locus, q_radius);

    // distance up from the quadric to the locus of possible points for the cam
    const double rho = distance * std::cos(angle);
    Eigen::Vector3d circ_center = q_loc + rho * q_nor;
    const double ell = std::sqrt(std::pow(distance, 2) - std::pow(rho, 2));

    EXPECT_NEAR(angle, expected_angle, 1e-1);
    // EXPECT_NEAR(distance, r_q_cam_world.norm(), 1e1);
    // const double rho1 = getRho(quad.projectToImageLocus(proj_mtx));
    // EXPECT_NEAR(h, r_q_cam_world.norm(), 1e1);

    putCraterIntoYaml(filename, q_loc, q_nor, q_radius, circ_center, ell, quad.getID(), 2);
  }
}

// TEST_F(NavigationTest, FindCircles) {
//     // Load image
//     cv::Mat src = imread("/home/ndvr/data/pds/cassini/image.png", cv::IMREAD_COLOR);
//     if (src.empty()) {
//         std::cout << "Could not open or find the image!" << std::endl;
//         return;
//     }

//     // Convert to grayscale
//     cv::Mat gray;
//     cvtColor(src, gray, cv::COLOR_BGR2GRAY);

//     // Blur image to reduce noise
//     cv::GaussianBlur(gray, gray, cv::Size(9, 9), 2, 2);

//     // Detect circles using HoughCircles
//     std::vector<cv::Vec3f> circles;
//     HoughCircles(gray, circles, cv::HOUGH_GRADIENT, 1,
//                  gray.rows / 8,  // Minimum distance between detected circles
//                  100, 30, 1, 30); // Canny edge threshold, accumulator threshold, min/max radius

//     // Draw detected circles on the original image
//     for (size_t i = 0; i < circles.size(); i++) {
//         cv::Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
//         int radius = cvRound(circles[i][2]);

//         // Draw the circle center
//         circle(src, center, 3, cv::Scalar(0, 0, 255), -1, 8, 0);

//         // Draw the circle outline
//         circle(src, center, radius, cv::Scalar(0, 255, 0), 2, 8, 0);
//     }

//     // Show the image with detected circles
//     cv::imshow("Detected Circles", src);
//     cv::waitKey(0);
// }

// TEST_F(NavigationTest, Circles) {

//     cv::Mat color = cv::imread("/home/ndvr/data/pds/cassini/image.png");
//     cv::namedWindow("input"); cv::imshow("input", color);

//     cv::Mat canny;

//     cv::Mat gray;
//     /// Convert it to gray
//     cv::cvtColor( color, gray, cv::COLOR_BGR2GRAY );

//     // compute canny (don't blur with that image quality!!)
//     cv::Canny(gray, canny, 200,20);
//     cv::namedWindow("canny2"); cv::imshow("canny2", canny>0);

//     std::vector<cv::Vec3f> circles;

//     /// Apply the Hough Transform to find the circles
//     double min_radius = 700, max_radius = 750;
//     cv::HoughCircles( gray, circles, cv::HOUGH_GRADIENT, 1, 60, 200, 20, min_radius, max_radius );

//     // Eigen::Matrix3d moon_locus = camera.getMoonConic(R_MOON);
//     // Conic moon(moon_locus);
//     /// Draw the circles detected
//     for( size_t i = 0; i < circles.size(); i++ )
//     {
//         cv::Point center(cvRound(circles[i][0]), cvRound(circles[i][1]));
//         std::cout << "Center: " << center << std::endl;
//         int radius = cvRound(circles[i][2]);
//         // if(radius < min_radius || radius > max_radius) {continue;}
//         cv::circle( color, center, 3, cv::Scalar(0,255,255), -1);
//         cv::circle( color, center, radius, cv::Scalar(0,0,255), 1 );
//     }

//     //compute distance transform:
//     cv::Mat dt;
//     cv::distanceTransform(255-(canny>0), dt, cv::DIST_L2 ,3);
//     cv::namedWindow("distance transform"); cv::imshow("distance transform", dt/255.0f);

//     // test for semi-circles:
//     float minInlierDist = 2.0f;
//     for( size_t i = 0; i < circles.size(); i++ )
//     {
//       // test inlier percentage:
//       // sample the circle and check for distance to the next edge
//       unsigned int counter = 0;
//       unsigned int inlier = 0;

//       cv::Point2f center((circles[i][0]), (circles[i][1]));
//       float radius = (circles[i][2]);

//       if(radius < min_radius || radius > max_radius) {continue;}

//       // maximal distance of inlier might depend on the size of the circle
//       float maxInlierDist = radius/25.0f;
//       if(maxInlierDist<minInlierDist) maxInlierDist = minInlierDist;

//       //TODO: maybe paramter incrementation might depend on circle size!
//       for(float t =0; t<2*3.14159265359f; t+= 0.1f)
//       {
//         counter++;
//         float cX = radius*cos(t) + circles[i][0];
//         float cY = radius*sin(t) + circles[i][1];

//         if(dt.at<float>(cY,cX) < maxInlierDist)
//         {
//           inlier++;
//           cv::circle(color, cv::Point2i(cX,cY),3, cv::Scalar(0,255,0));
//         }
//         else
//           cv::circle(color, cv::Point2i(cX,cY),3, cv::Scalar(255,0,0));
//       }
//       std::cout << 100.0f*(float)inlier/(float)counter << " % of a circle with radius " << radius << " detected" << std::endl;
//     }

//     cv::namedWindow("output"); cv::imshow("output", color);
//     cv::imwrite("houghLinesComputed.png", color);

//     cv::waitKey(-1);
// }
