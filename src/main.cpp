#include <iostream>
#include <Eigen/Dense>
#include "structs.h"
#include "crater-id.h"
#include "combinatorics.h"
#include "io.h"
#include "conics.h"
#include "quadrics.h"
#include <vtk-9.2/vtk3DS.h>
#include <spatialindex/capi/sidx_api.h>
// #include "RTree.h"
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
// #include "gnuplot-iostream.h"

#define NDIM 2

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef int ValueType;

typedef bg::model::point<float, NDIM, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, size_t> value;

struct Rect
{
  Rect()  {}

  Rect(int a_minX, int a_minY, int a_maxX, int a_maxY)
  {
    min[0] = a_minX;
    min[1] = a_minY;

    max[0] = a_maxX;
    max[1] = a_maxY;
  }
  Rect(int X, int Y)
  {
    min[0] = X;
    min[1] = Y;

    max[0] = X;
    max[1] = Y;
  }

  int min[2];
  int max[2];

  box calculate_bounding_box() {
    point point1(min[0], min[1]);
    point point2(max[0], min[1]);
    box bbox(point1, point2);
    return bbox;
  }

};
 
std::ostream& operator<<(std::ostream& os, const Rect& rectangle) {
    return os 
        // << std::fixed << std::setw(10) 
        // << std::setprecision(3) 
        // << std::setfill('_')
        << "->Rect: [" 
        << rectangle.min[0] << "," 
        << rectangle.min[1] << "] ["
        << rectangle.max[0] << "," 
        << rectangle.max[1] << "] ";
}
 
std::ostream& operator<<(std::ostream& os, const box& bbox) {
    return os 
        // << std::fixed << std::setw(10) 
        // << std::setprecision(3) 
        // << std::setfill('_')
        << "->Box: [" 
        << bg::get<bg::min_corner, 0>(bbox) << "," 
        << bg::get<bg::min_corner, 1>(bbox) << "] ["
        << bg::get<bg::max_corner, 0>(bbox) << "," 
        << bg::get<bg::max_corner, 1>(bbox) << "] ";
}

int main() {
    std::string fname;
    const std::string degrees = "°";
/*
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
    
/*
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
                  << std::endl;
        if(t_count > 20) {
            break;
        }
    }
    */
    Eigen::Vector3d surface_point1 = {30, 40, 50};
    Eigen::Vector3d surface_point2 = {0, 0, 10};
    Eigen::Vector3d surface_point3 = {1, 0, 0};
    getENUFrame(surface_point1);
    getENUFrame(surface_point2);
    getENUFrame(surface_point3);
    
    double radius = 100;
    Eigen::Vector3d position = {1000, 0., 0.};
    Quadric quad("TestQuadric", position, radius);
    std::cout << quad << std::endl;
    Conic conicA(10, 7, -100, -50, 0);
    Eigen::MatrixXd locus = conicA.getLocus();
    // Eigen::MatrixXd envelope = getMatrixAdjugate(locus);
    Eigen::Matrix4d test;
    test << 5,  -2,  2,  7,
            1,   0,  0,  3,
            -3,  1,  5,  0,
            3,  -1, -9,  4;
    // test << 1,2,3, 4,5,6, 7,8,9;
    Eigen::MatrixXd adjugate = getMatrixAdjugate(test);
    std::cout << "Adjugate:\n" << adjugate << std::endl;
/*
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

