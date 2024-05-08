#include "combinatorics.h"

template <typename T>
std::vector<std::vector<T>> Combination(const std::vector<T> choices, const uint r) {
  // https://stackoverflow.com/questions/9430568/generating-combinations-in-c
  const uint n = choices.size();
  std::vector<bool> v(n);
  std::vector<std::vector<T>> output;
  std::vector<T> subvector;
  std::fill(v.begin(), v.begin() + r, true);

  do {
    subvector.clear();
    uint idx = 0;
    for(auto& choice : choices) {
    // for (int i = 0; i < n; ++i) {
      if (v[idx]) {
        subvector.push_back(choice);
      }
      idx++;
    }
    output.push_back(subvector);
  } while (std::prev_permutation(v.begin(), v.end()));
  return output;
}

template <typename T>
void Permutation(std::vector<T> v)
{
  std::sort(v.begin(), v.end());
  do {
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
  } while (std::next_permutation(v.begin(), v.end()));
}

void formTriads(const std::vector<std::tuple<uint, uint>> pairs,
                std::vector<std::tuple<uint, uint, uint>>& triads) {
  uint ki, kj; 
  uint idx;
  std::vector<std::tuple<uint, uint>>::const_iterator iti, itj;

  idx = 0;
  for(const auto& [i, j] : pairs) {
    iti = std::find_if(pairs.begin()+idx, pairs.end(), [&i](const std::tuple<uint,uint>& e) {return std::get<0>(e) == i;});
    do {
      itj = std::find_if(pairs.begin()+idx, pairs.end(), [&j](const std::tuple<uint,uint>& e) {return std::get<0>(e) == j;});
      do {
        ki = std::get<1>(*iti);
        kj = std::get<1>(*itj);
        if(    ki == kj 
            && iti != pairs.end() 
            && std::get<0>(*iti) == i 
            && itj != pairs.end() 
            && std::get<0>(*itj) == j) {
          triads.push_back({i, j, ki});
        }
        itj++;
      } while(itj != pairs.end() && std::get<0>(*itj) == j);
      iti++;
    } while(iti != pairs.end() && std::get<0>(*iti) == i);
    idx++;
  }
}