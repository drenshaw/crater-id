#include <iostream>
#include <tuple>
#include <eigen3/Eigen/Dense>
#include <spatialindex/capi/sidx_api.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/program_options.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
// #include <boost/log/core.hpp>
// #include <boost/log/expressions.hpp>
// #include <boost/log/trivial.hpp>
// #include <Eigen/src/LU/InverseImpl.h>

#define BOOST_LOG_DYN_LINK 1
// #include <boost/log/utility/setup/console.hpp>

#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
#include "quadrics.h"
#include "visuals.h"
#include "camera.h"

#define RUN_ADJUGATE 0
#define RUN_VTK 0
#define RUN_RTREE 0
#define RUN_TRIADS 0
#define RUN_LOGGING 1
#define RUN_QUADRICS 1
#define RUN_CONICS 0
#define RUN_INVARIANTS 0

void plot_ellipse(cv::Mat& image, Conic& ellipse, const cv::Scalar color=cv::Scalar(0, 255, 0)) {
  cv::Point center;
  Eigen::Vector2d semiaxes;
  ellipse.GetCenter(center);
  ellipse.GetSemiAxes(semiaxes);
  cv::Size axes(semiaxes[0], semiaxes[1]);
  double angle = ellipse.GetAngle();
  cv::ellipse(image, center, axes, angle, 0, 360, color, -1, cv::LINE_AA);
  int font = cv::FONT_HERSHEY_SIMPLEX;
  float font_scale = 1; // scale factor from base size
  int thickness = 1; //in pixels
  cv::Scalar text_color(color[1], color[2], color[1]);
  cv::putText(image, std::to_string(ellipse.GetID()), center, font, font_scale,
              text_color, thickness, true);
}

void print_triads(const std::vector<std::tuple<uint, uint, uint>> triads, 
                  const std::vector<lunar_crater> craters,
                  const uint max_iter=10) {

  uint t_count = 0;
  for(const auto& [i, j, k] : triads) {
    std::cout << "-" << t_count++ << "-\t"
              << " "  << i
              << ", " << j
              << ", " << k
              << "\n\t" << craters[i]
              << "\n\t" << craters[j]
              << "\n\t" << craters[k]
              << std::endl;
    if(t_count >= max_iter) {
      break;
    }
  }
}

#if RUN_LOGGING
// namespace logging = boost::log;
// namespace keywords = boost::log::keywords;
// void init_logging() {     
//   // logging::add_console_log(std::clog, keywords::format = "%TimeStamp%: %Message%");
//   // logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::warning);
// }
#endif

#if RUN_ADJUGATE
template <typename Derived, int size>
Derived getCofactorTemplate(const Eigen::Matrix<Derived, size, size>& matrix, size_t cf_row, size_t cf_col) {
  size_t nrow = matrix.rows();
  size_t ncol = matrix.cols();
  size_t i = 0, j = 0;
  Eigen::MatrixX<Derived> cf_temp(nrow-1, ncol-1);

  // Looping for each element of the matrix
  for (size_t irow = 0; irow < nrow; irow++) {
    for (size_t icol = 0; icol < ncol; icol++) {
      //  Copying into temporary matrix only those
      //  element which are not in given row and
      //  column
      if (irow != cf_row && icol != cf_col) {
        cf_temp(i, j++) = matrix(irow, icol);

        // Row is filled, so increase row index and reset col index
        if (j == nrow - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
  return cf_temp.determinant();
}

template <typename Derived, int size>
Eigen::Matrix<Derived, size, size> getCofactorMatrixTemplate(const Eigen::Matrix<Derived, size, size>& matrix) {
  size_t nrow = matrix.rows();
  size_t ncol = matrix.cols();
  Eigen::MatrixX<Derived> cofactor_matrix(nrow, ncol);
  int cofactor_sign;
  for(size_t row=0;row<nrow; row++) {
    for(size_t col=0; col<ncol; col++) {
      cofactor_sign = (row+col)%2==0 ? 1 : -1;
      cofactor_matrix(row, col) = cofactor_sign*getCofactorTemplate(matrix, row, col);
    }
  }
  return cofactor_matrix;
}
#endif


int main(int argc, char** argv) {
#if RUN_ADJUGATE
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(4,4);
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(5,5);
    Eigen::MatrixXd madj = getAdjugateMatrix(m1);
    Eigen::MatrixXd madj2 = getCofactorMatrixTemplate(m2).transpose();
    double det = m1.determinant();
    Eigen::MatrixXd invv = m1.inverse();
    Eigen::MatrixXd m4 = det*invv;
    Eigen::MatrixXd m5 = det*m1.adjoint();
    // std::cout << "Matrix adjugate: \n" << std::endl << madj << std::endl;
    Eigen::MatrixXd res = madj - m4;
    
    double prec = sqrt(Eigen::NumTraits<double>::epsilon());
    if(!madj.isApprox(m4, prec)) {
    // if(res.isApproxToConstant(0.0)) {
    // if(res.norm() > prec) {
    // if(madj.isApprox(m4)) {
      std::cout << "Limits exceeded: "
      << std::endl << res << std::endl 
      << "\tNorm: " << res.norm() << std::endl
      << "\tEpsilon: " << prec << std::endl;
    }
#endif

#if RUN_LOGGING
  // init_logging();
    // BOOST_LOG_TRIVIAL(trace) << "A trace severity message";
    // BOOST_LOG_TRIVIAL(debug) << "A debug severity message";
    // BOOST_LOG_TRIVIAL(info) << "An informational severity message";
    // BOOST_LOG_TRIVIAL(warning) << "A warning severity message";
    // BOOST_LOG_TRIVIAL(error) << "An error severity message";
    // BOOST_LOG_TRIVIAL(fatal) << "A fatal severity message";
  // BOOST_LOG_TRIVIAL(error) << "Testing the log error";
#endif
#if RUN_VTK
  // VIS::Other();
  VIS::Sphere();
  // VIS::Cone();
  // VIS::Ellipsoid();
  // VIS::Cylinder();
  // VIS::HyperboloidOneSheet();
  // VIS::HyperboloidTwoSheets();
  // VIS::HyperbolicParaboloid();
  // VIS::EllipticParaboloid();
#endif

#if RUN_RTREE
  const std::string degrees = "°";
  // RTree
  std::vector<Rect> objects;
  objects.push_back(Rect(0, 0, 2, 2));
  objects.push_back(Rect(5, 5, 7, 7));
  objects.push_back(Rect(8, 5, 9, 6));
  objects.push_back(Rect(7, 1, 9, 2));
  objects.push_back(Rect(1, 8));
  objects.push_back(Rect(4, 7));
  objects.push_back(Rect(8, 2));
  objects.push_back(Rect(7, 5));
  objects.push_back(Rect(6, 6));
  objects.push_back(Rect(3, 7));
  objects.push_back(Rect(7, 7));
  objects.push_back(Rect(7, 9));

  // // Fill objects 

  // create the R* variant of the rtree
  bgi::rtree< value, bgi::rstar<16> > rtree;

  // insert some values to the rtree
  for ( size_t i = 0 ; i < objects.size() ; ++i )
  {
    // create a box
    box b = objects[i].calculate_bounding_box();
    // insert new value
    rtree.insert(std::make_pair(b, i));
  }

  // find values intersecting some area defined by a box
  box query_box(point(0, 0), point(5, 5));
  std::vector<value> result_s;
  rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
  std::cout << "Number of intersection results: " << result_s.size() << std::endl;
  for(auto& res : result_s) {
    std::cout << std::get<1>(res) << std::get<0>(res) << std::endl;
  }

  // find 5 nearest values to a point
  std::vector<value> result_n;
  rtree.query(bgi::nearest(point(0, 0), 5), std::back_inserter(result_n));
  std::cout << "Number of nearest neighbor results: " << result_n.size() << std::endl;
  for(auto& res : result_n) {
    std::cout << std::get<1>(res) << std::get<0>(res) << std::endl;
  }
#endif
    
#if RUN_TRIADS
  std::string fname;
  // Crater Reading
  // cout<<"Enter crater file name: ";
  // cin>>fname;
  // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
  fname = "/home/dqr0509/data/craters/lunar_craters.csv";

  std::vector<lunar_crater> craters;
  runCraterReader(fname, craters);
  // uint c_idx = 0;
  // for(auto& crater : craters) {
  //   std::cout << "Crater " << c_idx++ << " " << crater << std::endl;
  // }
  
  std::vector<std::tuple<uint, uint>> valids;
  float max_angle = 10.;
  specialCombination(craters, valids, max_angle);
  
  std::vector<std::tuple<uint, uint, uint>> triads;
  formTriads(valids, triads);

  std::cout << triads.size() << " valid triads found." << std::endl;
  // uint idx = 0;
  // for(const auto& [i, j, k]: triads) {
  //   std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
  // }

  int max_iter = 10;
  print_triads(triads, craters, max_iter);
#endif
  
#if RUN_QUADRICS
  /* QUADRICS */
  // lunar_crater crater = craters[0];
  // double radius = crater.diam/2;
  double radius = 30.0;
  // double r_radius_rim = calculateCraterRimFromRadius(radius);
  // Eigen::Vector3d position = {2.3, 1.25, -1.75};
  Eigen::Vector3d position = {0, 1, 2};
  // Eigen::Vector3d position = {0, 0, 0}; // should fail assert in quadric constructor
  // Eigen::Vector3d orientation = {4., 2.5, -1.5};
  Eigen::Vector3d orientation = {1, 0, 0};
  // orientation.normalize();
  // double lat = crater.lat;
  // double lon = crater.lon;
  // Quadric quad("TestQuadric", position, radius);
  // Quadric quad("LatLonQuadric", lat, lon, radius);
  Quadric quad("SmallQuadric", position, radius, orientation);
  // std::cout << "Location: (" << r_radius_rim*latlon2bearing(lat, lon).transpose() << ")" << std::endl;
  std::cout << quad << std::endl;

  Eigen::Matrix4d locus = quad.GetLocus();
  // locus(3,3) = 0;
  double maxVal = locus.cwiseAbs().maxCoeff();
  // maxVal = locus(0, 0);
  std::cout << "QLocus:\n" << locus/maxVal << std::endl;
#endif

#if RUN_CONICS
  std::vector<double> intrin = {10, 10, 511.5, 241.5, 0};
  Camera cam(intrin);

  cv::Mat image(500, 500, CV_8UC3,
                cv::Scalar(50, 50, 50));
  if (!image.data) {
    std::cout << "Could not open or find the image\n";
    return 0;
  }

  // const cv::Vec3b green(0, 255, 0);
  const cv::Scalar green(0, 255, 0);
  const cv::Scalar light_blue(255, 190, 150);
  const cv::Scalar purple(180, 15, 160);
  const cv::Scalar red(0, 0, 255);
  Conic conicA(10, 7, 300, 50, 0);
  Conic conicB(15, 12, 100, 200, 0);
  Conic conicC(12, 8, 50, 200, 0);
  Conic conicD(12, 8, 400, 20, 0);
  plot_ellipse(image, conicA, green);
  plot_ellipse(image, conicB, light_blue);
  plot_ellipse(image, conicC, purple);
  plot_ellipse(image, conicD, red);

  // Showing image inside a window
  cv::imshow("Output", image);
  cv::waitKey(0);
  assert(conicA == conicA);
  assert(!(conicA == conicB));
  assert(conicA != conicB);
  std::cout << "Does this match? It should: " << (conicA == conicA) << std::endl;
  std::cout << "Does this match? It shouldn't: " << (conicA == conicB) << std::endl;
  std::cout << "IDs? " << conicA.get_id() << " | " << conicB.get_id() << std::endl;
#endif


  int sz[3] = {2,2,2};
  cv::Mat L(3,sz, CV_8UC(1), cv::Scalar::all(0));
  std::cout << L.size() << std::endl;
  // std::cout << L.row(0).col(0).row(0) << std::endl;
  // // Eigen::MatrixXd locus = conicA.GetLocus();
  // // Eigen::MatrixXd envelope = getAdjugateMatrix(locus);
  // Eigen::Matrix4d test;
  // test << 5,  -2,  2,  7,
  //         1,   0,  0,  3,
  //         -3,  1,  5,  0,
  //         3,  -1, -9,  4;
  // // test << 1,2,3, 4,5,6, 7,8,9;
  // Eigen::MatrixXd adjugate = getAdjugateMatrix(test);
  // // std::cout << "Adjugate:\n" << adjugate << std::endl;

#if RUN_INVARIANTS
  // Invariants
  std::vector<double> invariantsABC, invariantsBCD, invariantsCDA, invariantsDAB;
  if(!computeCraterTriadInvariants(conicA, conicB, conicC, invariantsABC)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicB, conicC, conicD, invariantsBCD)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicC, conicD, conicA, invariantsCDA)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  if(!computeCraterTriadInvariants(conicD, conicA, conicB, invariantsDAB)) {
    std::cerr << "Error in `computeCraterTriadInvariants`" << std::endl;
  }
  printVector(invariantsABC, "Invariants ABC: ");
  printVector(invariantsBCD, "Invariants BCD: ");
  printVector(invariantsCDA, "Invariants CDA: ");
  printVector(invariantsDAB, "Invariants DAB: ");
#endif

    return 0;
}

