#include "gtest/gtest.h"
#include <eigen3/Eigen/Geometry>
#include <opencv2/imgproc.hpp> 
#include <opencv2/highgui/highgui.hpp> 
// #include <random>

#include "quadrics.h"
#include "camera.h"
#include "visuals.h"
#include "math_utils.h"
#include "conics.h"

#define R_MOON 1737.4

class QuadricTest : public testing::Test {
protected:
  QuadricTest() {
    latitude = 0, longitude = 0, radius = 200;
    id = "defaultQuadric";
    quadric_default = new Quadric(latitude, longitude, radius, id);
    // This is a real crater with lat/lon and radius from a lunar crater
    quadric_real    = new Quadric(61.6941, 172.365, 83.2894, "01-1-000026");
    dx = 1000, dy = 1000, skew = 0, im_height = 2048, im_width = 2592;
    up = (double(im_width)+1)/2;
    vp = (double(im_height)+1)/2;
    image_size = cv::Size(im_width, im_height);
    origin = Eigen::Vector3d::Zero();
    quat = Eigen::Quaterniond::Identity();
    position << 1e2,2e2,-3e4;
    cam = new Camera(dx, dy, up, vp, skew, image_size, quat, position);
  }
  ~QuadricTest() override {
    delete quadric_default;
  }

public:
  // ~QueueTest() override = default;
  double latitude, longitude, radius;
  std::string id;
  Quadric* quadric_default;
  Quadric* quadric_real;
  Eigen::Vector3d origin;
  double dx, dy, up, vp, skew;
  int im_height, im_width;
  cv::Size2i image_size;
  Camera* cam;
  Eigen::Quaterniond quat;
  Eigen::Vector3d position;
};

TEST_F(QuadricTest, MakingQuadric) {
  Eigen::Matrix3d t_e2m;
  // For clarity: we use a passive transformation, not an active rotation
  //   local x col vector is component in the Y-direction
  //   local y col vector is component in the Z-direction
  //   local z col vector is component in the X-direction
  t_e2m <<  0, 0, 1, 
            1, 0, 0, 
            0, 1, 0;
  ASSERT_TRUE(quadric_default->getQuadricTransformationMatrix().isApprox(t_e2m));
  ASSERT_EQ(quadric_default->getID(), id);
}

TEST_F(QuadricTest, InvalidCraterRadius) {
  double invalid_radius = 17000;
  double zero_radius = 0;
  std::string id_excess = "excessRadius";
  std::string id_zero = "zeroRadius";
  EXPECT_THROW(Quadric quad(latitude, longitude, invalid_radius, id_excess), std::runtime_error);
  EXPECT_THROW(Quadric quad(latitude, longitude, zero_radius, id_zero), std::runtime_error);
}

TEST_F(QuadricTest, InvalidPositionWithDependentNormal) {
  std::string id = "dependentNormal";
  EXPECT_THROW(Quadric quad(origin, radius, id), std::runtime_error);
}

TEST_F(QuadricTest, InvalidSurfaceNormal) {
  Eigen::Vector3d position = 1e5*latlon2bearing(latitude, longitude);
  Eigen::Vector3d surface_normal(0, 0, 0);
  std::string id = "zeroNormal";
  EXPECT_THROW(Quadric quad(position, radius, surface_normal, id), std::runtime_error);
}

TEST_F(QuadricTest, AngleBetweenQuadrics) {
  double lat1 =  0, lon1 = 0, radius1 = 50;
  double lat2 = 30, lon2 = 0, radius2 = 300;
  Quadric q1(lat1, lon1, radius1);
  Quadric q2(lat2, lon2, radius2);
  Quadric q1_copy(lat1, lon1, radius1);
  double angle12 = q1.getAngleBetweenQuadrics(q2);
  ASSERT_DOUBLE_EQ(rad2deg(angle12), std::abs(lat2-lat1));
  double angle_zero = q1.getAngleBetweenQuadrics(q1_copy);
  ASSERT_DOUBLE_EQ(angle_zero, 0.);
}

TEST_F(QuadricTest, AxisOfRotationQuadrics) {
  double lat1 =  0, lon1 = 0, radius1 = 50;
  double lat2 = 30, lon2 = 0, radius2 = 50;
  Quadric q1(lat1, lon1, radius1);
  Quadric q2(lat2, lon2, radius2);
  Quadric q1_copy(lat1, lon1, radius1);
  std::optional<Eigen::Vector3d> axis_normal = q1.getAxisNormalToQuadrics(q2);
  Eigen::Vector3d axis_check(3);
  // Moving from 0 degrees latitude up to 30 degrees latitude makes axis in -y direction
  axis_check << 0, -1, 0;
  ASSERT_TRUE(axis_normal.has_value());
  ASSERT_TRUE(axis_normal.value().isApprox(axis_check));
  ASSERT_FALSE(q1.getAxisNormalToQuadrics(q1_copy).has_value());
}

TEST_F(QuadricTest, QuadricFrom3Points) {

  Eigen::Vector3d p1 = {10,  10, 0};
  Eigen::Vector3d p2 = {10,  00, 10};
  Eigen::Vector3d p3 = {10, -10, 0};

  Quadric quad(p1, p2, p3, "threePoints");
  Eigen::Vector3d center = quad.getLocation();
  double radius = quad.getRadius();
  ASSERT_DOUBLE_EQ(radius, 10);
  Eigen::Vector3d center_check = {10, 0, 0};
  ASSERT_TRUE(center.isApprox(center_check));
  Eigen::Vector3d normal_check = Eigen::Vector3d::UnitX();
  EXPECT_TRUE(quad.getPlane().normal().isApprox(normal_check) || 
              quad.getPlane().normal().isApprox(-normal_check));
}

TEST_F(QuadricTest, SamePlane) {
  Eigen::Vector3d p1 = {10,  10,   0};
  Eigen::Vector3d p2 = {10,   0,  10};
  Eigen::Vector3d p3 = {10, -10,   0};
  Eigen::Vector3d p4 = {10,   0, -10.00000001};

  Quadric quad1(p1, p2, p3, "threePoints1");
  Quadric quad2(p1, p2, p4, "threePoints2");
  ASSERT_TRUE(isSamePlane(quad1.getPlane(), quad2.getPlane()));
}

TEST_F(QuadricTest, RadiusCheckFromPoints) {
  double r = 100;
  double rho = calculateCraterRimFromRadius(r);
  double phi0 = deg2rad(10.), phi1 = deg2rad(115.), phi2 = deg2rad(280.);
  Eigen::Vector3d p1 = {rho, r*std::cos(phi0), r*std::sin(phi0)};
  Eigen::Vector3d p2 = {rho, r*std::cos(phi1), r*std::sin(phi1)};
  Eigen::Vector3d p3 = {rho, r*std::cos(phi2), r*std::sin(phi2)};

  Quadric quad(p1, p2, p3, "fromRadius");
  Quadric quad_check(0, 0, r, "pts_check");

  EXPECT_DOUBLE_EQ(r, quad.getRadius());

  EXPECT_TRUE(isSamePlane(quad, quad_check));
  EXPECT_TRUE(quad == quad_check);

}

TEST_F(QuadricTest, InPlaneXYQuadric) {
  double r = 1e4;
  double z = 0;
  double phi0 = deg2rad(10.), phi1 = deg2rad(115.), phi2 = deg2rad(280.);
  Eigen::Vector3d p1 = {r*std::sin(phi0), r*std::cos(phi0), z};
  Eigen::Vector3d p2 = {r*std::sin(phi1), r*std::cos(phi1), z};
  Eigen::Vector3d p3 = {r*std::sin(phi2), r*std::cos(phi2), z};
  Quadric quad_inplane(p1, p2, p3, "inPlane");
  Eigen::Hyperplane<double, 3> plane(Eigen::Vector3d(0,0,-1),0);
  EXPECT_TRUE(isSamePlane(plane, quad_inplane.getPlane()));
}

TEST_F(QuadricTest, InPlaneXZQuadric) {
  double r = 1e4;
  double y = 0;
  // TODO: there are some instances where points like (1,0,0),(0,1,0), and (-1,0,0)
  // fail to form a valid plane, which is strange
  Eigen::Vector3d p1 = {0, y,  r};
  Eigen::Vector3d p2 = {r, y,  0};
  Eigen::Vector3d p3 = {0, y, -r};
  Quadric quad_inplane(p1, p2, p3, "inStrangePlane");
  Eigen::Hyperplane<double, 3> plane(Eigen::Vector3d(0,-1,0),0);
  EXPECT_TRUE(isSamePlane(plane, quad_inplane.getPlane()));
}

bool checkAdjugacy(const Eigen::Matrix4d& mtx) {
  Eigen::Matrix4d adj_mtx = adjugate(mtx);
  Eigen::Matrix4d adj_adj = adjugate(adj_mtx);

  double det = mtx.determinant();
  int n_dim = mtx.rows();
  double scalar = std::pow(det,n_dim-2);
  // std::cout << "\nMatrix:\n" << mtx << std::endl;
  // std::cout << "\nMatrix adjugate:\n" << adj_mtx << std::endl;
  // std::cout << "\nAA*:\n" << mtx * adj_mtx/det << std::endl;
  // std::cout << "det^(n-2)*mtx (rank deficient)-> det: " << det << std::endl
  //           << scalar * mtx << std::endl;
  return (scalar * mtx).isApprox(adj_adj);
}

TEST_F(QuadricTest, RealCraterData) {

  Eigen::Matrix4d envelope_check;
  envelope_check << 659792,     -89375.2, -1.24902e+06,     -815.596,
                  -89375.2,      5043.65,       167431,      109.331,
              -1.24902e+06,       167431,  2.33291e+06,       1527.9,
                  -815.596,      109.331,       1527.9,            1;
  Eigen::Matrix4d envelope = quadric_real->getEnvelope();
  envelope /= envelope(Eigen::last, Eigen::last);  
  // normalizeDeterminant(envelope);                
  EXPECT_TRUE(envelope.isApprox(envelope_check, 1e-3));
  // std::cout << *quadric_real << std::endl;
  // Since this is a rank-deficient matrix, expect this to fail check
  EXPECT_FALSE(checkAdjugacy(envelope));
}

TEST_F(QuadricTest, SampleCrater) {
  Quadric lat0lon0(0, 0, 100, "SittingPretty@0,0");
  EXPECT_TRUE(lat0lon0.getNormal().isApprox(Eigen::Vector3d(1,0,0)));
  Eigen::Matrix4d envelope = lat0lon0.getEnvelope();
  envelope /= envelope(Eigen::last, Eigen::last);  
  // normalizeDeterminant(envelope);                
  // EXPECT_TRUE(envelope.isApprox(envelope_check, 1e-3));
  // std::cout << lat0lon0 << std::endl;
  // Since this is a rank-deficient matrix, expect this to fail check
  EXPECT_FALSE(checkAdjugacy(envelope));
}

TEST_F(QuadricTest, EnvelopeToLocus) {
  // TODO: Since the envelope will always be rank 3, we will never abide by
  // the typical adjugate matrix rules; for now, we never use the quadric
  // locus, so we will ignore this for now.
}

void extractEllipseParameters(const Eigen::Matrix3d& A, double& semiMajor, double& semiMinor, Eigen::Vector2d& center, double& angle) {
    // Extract the top-left 2x2 submatrix
    Eigen::Matrix2d subMatrix = A.topLeftCorner<2, 2>();
    
    // Eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(subMatrix);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue decomposition failed!" << std::endl;
        return;
    }
    
    // Eigenvalues and eigenvectors
    Eigen::Vector2d eigenvalues = eigensolver.eigenvalues();
    Eigen::Matrix2d eigenvectors = eigensolver.eigenvectors();
    
    // Semi-major and semi-minor axes
    // std::cout << "Eigenvalues: " << eigenvalues.transpose() << std::endl;
    semiMajor = 1.0 / std::sqrt(eigenvalues.minCoeff());
    semiMinor = 1.0 / std::sqrt(eigenvalues.maxCoeff());
    
    // Center of the ellipse
    center = -subMatrix.inverse() * A.topRightCorner<2, 1>();
    
    // Angle with respect to the x-axis
    angle = std::atan2(eigenvectors(1, 0), eigenvectors(0, 0));
}

TEST_F(QuadricTest, ProjectQuadric) {
  double radius = 50;
  double lat = -89;
  std::vector<Quadric> quadrics;
  for(int i = 0; i < 24; i++) {
    Quadric quad(lat+i*1.5, i*13, radius, "");
    quadrics.push_back(quad);
  }

  Eigen::Vector3d pos(0, 0,-2.5e3);
  cam->moveTo(pos);
  cam->pointTo(origin, -Eigen::Vector3d::UnitY());
  // Eigen::MatrixXd proj_mtx = cam->getProjectionMatrix();
  // viz::drawCraters(*cam, quadrics);
}

TEST_F(QuadricTest, MakeSphere) {
  Eigen::Vector3d pos(2.2e3, 0 ,0);
  cam->setPosition(pos);
  Eigen::Vector3d point_to(double(R_MOON)*Eigen::Vector3d::UnitZ());
  cam->pointTo(point_to, Eigen::Vector3d::UnitZ());

  double radius = 50;
  double dist = 0.3;
  std::vector<Quadric> craters;
  std::vector<Conic> conics;
  for(int i = 0; i < 24; i++) {
    double theta = deg2rad(i*13);
    double lat = dist*rad2deg(std::cos(theta));
    double lon = dist*rad2deg(std::sin(theta));
    Quadric quad(lat , lon, radius, "");
    craters.push_back(quad);
  }
  // viz::drawCraters(*cam, craters);
}

TEST_F(QuadricTest, ProjectMoonCenter) {
  Eigen::Vector3d pos(2.5e3, 0, 0);
  // cam->moveCamera(pos);
  cam->moveTo(pos);
  cam->pointTo(origin, Eigen::Vector3d::UnitZ());
  // std::cout << "Origin: " << origin.transpose() << std::endl;
  
  // Add the moon
  Eigen::Matrix4d sphere = makeSphere(double(R_MOON));
  Eigen::Vector2d pixel;
  cv::Mat image = cam->getBlankCameraImage();
  cv::Mat outImg;
  // double scaling = 0.4;
  double semiMajor, semiMinor, angle;
  Eigen::Vector2d center;
  Eigen::AngleAxisd rot = Eigen::AngleAxisd(-M_PI / 48, Eigen::Vector3d::UnitX());
  for(int i = 0; i < 10; i++) {
    cam->resetImage(image);
    Eigen::Matrix3d moon_locus = cam->projectQuadricToLocus(sphere);
    extractEllipseParameters(moon_locus, semiMajor, semiMinor, center, angle);
    Conic moon(moon_locus);
    // std::cout << "New params: " << semiMajor << ", " << semiMinor << ", " << center.transpose() << ", " << angle << std::endl;
    cam->world2Pixel(Eigen::Vector3d::Zero(), pixel);
    // std::cout << "Moon center: " << pixel.transpose() << std::endl;

    // viz::drawEllipse(image, moon, cv::viz::Color::gray());
    // // Showing image inside a window 
    // cv::resize(image, outImg, cv::Size(), scaling, scaling);
    // cv::imshow("Projecting Moon to Camera", outImg); 
    // cv::waitKey(0); 

    cam->rotate(rot);
  }
}

TEST_F(QuadricTest, QuadricPoints) {
  int n_pts = 10;
  std::vector<Eigen::Vector3d> pts_cam;
  std::vector<Eigen::Vector2d> pxl_cam;
  pts_cam.reserve(n_pts);
  Quadric quad(35, 25, 300, "getpoints");
  quad.getRimPoints(n_pts, pts_cam);
  cam->setPosition(Eigen::Vector3d(3e3, 2.0e3, 1e3));
  cam->pointTo(Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ());
  cv::Mat image = cam->getBlankCameraImage();
  std::vector<Eigen::Vector2d> pts_pxl;
  pts_pxl.reserve(n_pts);
  Conic cc_from_proj = quad.projectToImage(cam->getProjectionMatrix());
  // std::cout << "toImage: " << cc_from_proj << std::endl;
  // Eigen::Matrix3d conic_locus = cam->projectQuadricToLocus(quad.getLocus());
  // Conic con;
  // con.fromLocus(conic_locus);
  // std::cout << con << std::endl;

  for (auto pt_it = pts_cam.begin(); pt_it != pts_cam.end(); pt_it++) {
    // int index = std::distance(pts_cam.begin(), pt_it);
    Eigen::Vector2d pt_pxl;
    cam->world2Pixel(*pt_it, pt_pxl);
    pts_pxl.push_back(pt_pxl);
  }
  Conic cc_from_pts(pts_pxl);
  EXPECT_EQ(cc_from_proj, cc_from_pts);
  // Add the moon
  Eigen::Matrix4d sphere = makeSphere(double(R_MOON));
  Eigen::Matrix3d moon_locus = cam->projectQuadricToLocus(sphere);
  Conic moon(moon_locus);
  // std::cout << "Moon ellipse: " << moon << std::endl;
  if(cam->isInCameraFrame(moon.getCenter())) {
    viz::drawEllipse(image, moon, cv::viz::Color::gray());
  }
  else {
    std::cout << "The Moon center is not in the image: " << moon.getCenter().transpose() << std::endl;
  }
  viz::draw3dAxes(image, *cam);
  viz::drawPoints(image, pts_pxl, viz::CV_colors);
  // Showing image inside a window 
  double scaling = 0.5;
  cv::Mat outImg;
  cv::resize(image, outImg, cv::Size(), scaling, scaling);
  // viz::interactiveZoom(outImg);
  // // cv::imshow("Points", outImg); 
  // COpenCVWindowExt window ("Moon stuff"); 
  // window.ImShow (outImg);
  // // window.ImRead("/home/ndvr/Pictures/IMG_20230904_121221.jpg");
  // cv::waitKey(0); 
}

TEST_F(QuadricTest, ProjectCrater) {

  Eigen::Vector3d pos(2.4e3, 0, 0);
  // cam->moveCamera(pos);
  cam->moveTo(pos);
  cam->pointTo(origin, Eigen::Vector3d::UnitZ());
  
  // Add the moon
  Quadric zero(0,0,300,"@zero");
  Eigen::Vector2d pixel;
  cv::Mat image = cam->getBlankCameraImage();
  cv::Mat outImg;
  double scaling = 0.5;
  Eigen::AngleAxisd rot = Eigen::AngleAxisd(-M_PI / 36, Eigen::Vector3d::UnitX());
  const uint n_pts = 4;
  for(int i = 0; i < 10; i++) {
    cam->resetImage(image);
    Eigen::MatrixXd proj = cam->getProjectionMatrix();
    Eigen::Matrix3d locus = zero.projectToImageLocus(proj);
    // Eigen::Matrix3d locus = cam->projectQuadricToLocus(zero.getLocus());
    // std::cout << "Count:\n" << zero.getLocus() << std::endl;
    // std::cout << "Locus:\n" << locus << std::endl;
    Conic crater(locus);
    Eigen::Vector2d pixel_n, pixel_s, pixel_e, pixel_w;

    std::vector<Eigen::Vector3d> pts_world;
    std::vector<Eigen::Vector2d> pts_pxl;
    zero.getRimPoints(n_pts, pts_world);
    cam->world2Pixel(pts_world, pts_pxl);

    viz::drawEllipse(image, crater, cv::viz::Color::gray());
    std::vector<cv::viz::Color> colors = {cv::viz::Color::red(), cv::viz::Color::cyan(), cv::viz::Color::orange(), cv::viz::Color::blue()};
    viz::drawPoints(image, pts_pxl, viz::CV_colors);
    // Showing image inside a window 
    cv::resize(image, outImg, cv::Size(), scaling, scaling);
    // COpenCVWindowExt window ("Src"); 
    // // window.ImShow (outImg);
    // window.ImRead("/home/ndvr/Pictures/IMG_20230904_121221.jpg");
    cv::imshow("Projecting Moon to Camera", outImg); 
    cv::waitKey(0); 

    cam->rotate(rot);
  }
}

TEST_F(QuadricTest, Reprojection) {

}

TEST_F(QuadricTest, ProjectionMatrix) {
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
}

TEST_F(QuadricTest, ParameterExtraction) {
  // // Example 3x3 matrix representation of an ellipse
  //   Eigen::Matrix3d A;
  //   A << 5, 1, -2,
  //        1, 4, -3,
  //        -2, -3, 1;
  //   Eigen::Matrix3d A_adj = A.adjoint();
  //   std::cout << "Adj(A):\n" << A_adj << std::endl;
  //   if (!isEllipse(A)) {
  //       std::cerr << "The matrix does not represent an ellipse." << std::endl;
  //       ASSERT_TRUE(false);
  //   }
    
  //   double semiMajor, semiMinor, angle;
  //   Eigen::Vector2d center;
    
  //   extractEllipseParameters(A, semiMajor, semiMinor, center, angle);
    
  //   std::cout << "Semi-major axis: " << semiMajor << std::endl;
  //   std::cout << "Semi-minor axis: " << semiMinor << std::endl;
  //   std::cout << "Center: (" << center.x() << ", " << center.y() << ")" << std::endl;
  //   std::cout << "Angle with respect to x-axis: " << angle << " radians" << std::endl;
    
  //   ASSERT_TRUE(true);
}
