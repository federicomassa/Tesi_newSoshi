#ifndef COMBINATIONS_H
#define COMBINATIONS_H

#include "TMath.h"
#include <vector>
#include <iostream>
#include <ParticleAnalysis/combinations.h>

template<class T> void go(const std::vector<T>& initVector, std::vector<T>& combVector, std::vector<std::vector<T> >& resultVector, int offset, int k) {

if (k == 0) {
    resultVector.push_back(combVector);
    return;
  }
  for (int i = offset; i <= initVector.size() - k; ++i) {
    combVector.push_back(initVector[i]);
    go<T>(initVector, combVector, resultVector, i+1, k-1);
    combVector.pop_back();
  }
}

template<class T> std::vector< std::vector<T> > combinations (const std::vector<T>& initVector, int k) {

  std::vector<T> combVector;
  std::vector<std::vector<T> > resultVector;

  if (initVector.size() < k) {
    std::cout << "ERROR\t in combinations(...): class is smaller than vector size. Exiting" << std::endl;
    exit(1);
  }
 
  go<T>(initVector, combVector, resultVector, 0, k);
  
  return resultVector;
}

// void prova() {

//   std::vector<int> test;
  
//   test.push_back(1);
//   test.push_back(2);
//   test.push_back(3);
//   test.push_back(4);
//   test.push_back(5);

//   std::vector<std::vector<int> > result = combinations<int>(test, 3);

//   for (std::vector<std::vector<int> >::iterator itr = result.begin(); itr != result.end(); itr++) {
//     for (std::vector<int>::iterator itrEl = (*itr).begin(); itrEl != (*itr).end(); itrEl++)
//       std::cout << (*itrEl) << '\t';
//     std::cout << std::endl;
//   }

// }

#endif
