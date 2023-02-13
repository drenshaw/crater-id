#include <iostream>
#include <Eigen/Dense>
#include <spatialindex/capi/sidx_api.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/program_options.hpp>
#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
#include "quadrics.h"
#include "visuals.h"
#include "camera.h"

// #include <vtk3DS.h>
// #include <vtkActor.h>
// #include "gnuplot-iostream.h"
namespace po = boost::program_options;
using namespace boost::program_options;

void on_age(int age)
{
  std::cout << "On age: " << age << '\n';
}

int main(int argc, char** argv) {
  try
  {
    boost::program_options::options_description desc
        ("\nMandatory arguments marked with '*'.\n"
           "Invocation : <program> --host <hostname> --port <port> <web_app_name> <web_app_schema_file> \nAgruments");
    // po::options_description tester_options("Tester options");
    // boost::program_options::options_description desc{"Options"};
    // desc.add_options()
    //   ("help,h", "Help screen");
    //   ("pi", value<float>()->default_value(3.14f), "Pi")
    //   ("age", value<int>()->notifier(on_age), "Age");

    // variables_map vm;
    // store(parse_command_line(argc, argv, desc), vm);
    // notify(vm);

    // if (vm.count("help"))
    //   std::cout << desc << '\n';
    // else if (vm.count("age"))
    //   std::cout << "Age: " << vm["age"].as<int>() << '\n';
    // else if (vm.count("pi"))
    //   std::cout << "Pi: " << vm["pi"].as<float>() << '\n';
  }
  catch (const error &ex)
  {
    std::cerr << ex.what() << '\n';
  }
    // VIS::Other();
    // VIS::Sphere();
    // VIS::Cone();
    // VIS::Ellipsoid();
    // VIS::Cylinder();
    // VIS::HyperboloidOneSheet();
    // VIS::HyperboloidTwoSheets();
    // VIS::HyperbolicParaboloid();
    // VIS::EllipticParaboloid();
  

    std::string fname;
    const std::string degrees = "°";
/*    // RTree
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
*/
    

    // Crater Reading
    // cout<<"Enter crater file name: ";
    // cin>>fname;
    // fname = "/home/ndvr/Documents/Research/CraterDatabase/crater_reduced.csv";
    fname = "/home/dqr0509/data/craters/lunar_craters.csv";

    std::vector<lunar_crater> craters;
    runCraterReader(fname, craters);
    // uint c_idx = 0;
    // for(auto& crater : craters) {
    //     std::cout << "Crater " << c_idx++ << " " << crater << std::endl;
    // }
    
    std::vector<std::tuple<uint, uint>> valids;
    specialCombination(craters, valids, 10.);
    
    std::vector<std::tuple<uint, uint, uint>> triads;
    formTriads(valids, triads);
    std::cout << triads.size() << " valid triads found." << std::endl;
    // uint idx = 0;
    // for(const auto& [i, j, k]: triads) {
    //     std::cout << "IDX: " << idx++ << " | " << i << ", " << j << ", " << k << std::endl;
    // }

    uint t_count = 0;
    for(const auto& [i, j, k] : triads) {
        std::cout << "-" << t_count++ << "-\t"
                  << " "  << i
                  << ", " << j
                  << ", " << k
                  << "\t" << craters[i]
                  << "\t" << craters[j]
                  << "\t" << craters[k]
                  << std::endl;
        if(t_count > 10) {
            break;
        }
    }
    
    // lunar_crater crater = craters[0];
    // double radius = crater.diam/2;
    // double r_radius_rim = calculateCraterRimFromRadius(radius);
    Eigen::Vector3d position = {2., 1.25, -1.25};
    // Eigen::Vector3d orientation = {4., 2.5, -1.5};
    // orientation.normalize();
    // double lat = crater.lat;
    // double lon = crater.lon;
    // Quadric quad("TestQuadric", position, radius);
    // Quadric quad("LatLonQuadric", lat, lon, radius);
    Quadric quad("SmallQuadric", position, 3.0);
    // std::cout << "Location: (" << r_radius_rim*latlon2unitVector(lat, lon).transpose() << ")" << std::endl;
    std::cout << quad << std::endl;

    Eigen::Matrix4d locus = quad.getLocus();
    // locus(3,3) = 0;
    double maxVal = locus.cwiseAbs().maxCoeff();
    maxVal = locus(0, 0);
    std::cout << "QLocus:\n" << locus/maxVal/2.5 << std::endl;

    std::vector<double> intrin = {10, 10, 511.5, 241.5, 0};
    Camera cam(intrin);


    // // Conic conicA(10, 7, -100, -50, 0);
    // // Eigen::MatrixXd locus = conicA.getLocus();
    // // Eigen::MatrixXd envelope = getMatrixAdjugate(locus);
    // Eigen::Matrix4d test;
    // test << 5,  -2,  2,  7,
    //         1,   0,  0,  3,
    //         -3,  1,  5,  0,
    //         3,  -1, -9,  4;
    // // test << 1,2,3, 4,5,6, 7,8,9;
    // Eigen::MatrixXd adjugate = getMatrixAdjugate(test);
    // // std::cout << "Adjugate:\n" << adjugate << std::endl;
/*    // Invariants 
    // Invariants
    Conic conicB(15, 12, 100, 200, 0);
    Conic conicC(12, 8, 50, -200, 0);
    Conic conicD(12, 8, -500, -20, 0);
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
    printVector(invariantsABC, "Invariants: ");
    printVector(invariantsBCD, "Invariants: ");
    printVector(invariantsCDA, "Invariants: ");
    printVector(invariantsDAB, "Invariants: ");
*/
    return 0;
}

