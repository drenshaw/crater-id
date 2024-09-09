#include "gtest/gtest.h"
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

#include "conics.h"
// #include "io.h"

template <typename T>
void printVector(const std::vector<T> vec, const std::string prepend) {
  std::cout << prepend;
  for(auto& idx : vec) {
    std::cout << idx << ", ";
  }
  std::cout << std::endl;
}

template <typename T, size_t SIZE>
void printVector(const std::array<T, SIZE> arr, const std::string prepend) {
  std::cout << prepend;
  for(auto& elem : arr) {
    std::cout << elem << ", ";
  }
  std::cout << std::endl;
}

class InvariantTest : public testing::Test {
  protected:
    InvariantTest() {
      // std::array<double, 5> arr = {10.0, 7.0, 300.0, 50.0, 0.0};
      // Conic conic_arr(arr);
      // c_arr.SetGeometricParameters(arr);
    }
  public:
      // // ~QueueTest() override = default;
  // Queue<int> q0_;
  // Queue<int> q1_;
  // Queue<int> q2_;
};

TEST(InvariantTest, InvariantTriad) {
  Conic conicA(10, 7, 300, 50, 0);
  Conic conicB(15, 12, 100, 200, 0);
  Conic conicC(12, 8, 50, 200, 0);
  Conic conicD(12, 8, 400, 20, 0);
  
  // Invariants
  std::array<double, NONCOPLANAR_INVARIANTS> invariantsABC, invariantsBCD, invariantsCDA, invariantsDAB;
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
}

