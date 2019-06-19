
#ifndef SH_WRAP_HPP
#define SH_WRAP_HPP

#include <string>
#include <vector>
#include <iostream>
//#include <numpy/arrayobject.h>
#include "sh/spherical_harmonics.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <memory>
#include <algorithm>

#include <omp.h>

class SH_Wrapper {
public:
  SH_Wrapper() : SH_Wrapper(0) {};

  SH_Wrapper(int _verbose);

  int run(
    int K,
    std::vector< double > vecs,
    std::vector< double > scalars,
    std::vector< int > vecs_shape,
    std::vector< int > scalars_shape,
    std::vector< double >& output,
    bool only_even=false
  );

  int run_all_same_vecs(
    int K,
    std::vector< double > vecs,
    std::vector< double > scalars,
    std::vector< int > vecs_shape,
    std::vector< int > scalars_shape,
    std::vector< double >& output,
    bool only_even=false
  );

  int run_even(
    int K,
    std::vector< double > vecs,
    std::vector< double > scalars,
    std::vector< int > vecs_shape,
    std::vector< int > scalars_shape,
    std::vector< double >& output,
    bool all_same_vecs=false
  );

  int project_from_SH(
    int K,
    std::vector<double> sh_coefs,
    std::vector<double> vecs,
    std::vector<int> vecs_shape,
    std::vector<double>& output
  );

  int project_from_even_SH(
    int K,
    std::vector<double> sh_coefs,
    std::vector<double> vecs,
    std::vector<int> vecs_shape,
    std::vector<double>& output
  );

  int projection_matrix(
    int K,
    std::vector<double> vecs,
    std::vector<int> vecs_shape,
    std::vector<double>& output,
    std::vector<int>& out_shape
  );

  int projection_matrix_even(
    int K,
    std::vector<double> vecs,
    std::vector<int> vecs_shape,
    std::vector<double>& output,
    std::vector<int>& out_shape
  );

  int check_num_coef(
    int K
  );

  int check_num_even_coef(
    int K
  );

private:

  int n_vox;
  int n_points;

  int verbose;

};







#endif



