
#include "sh_wrap.hpp"

SH_Wrapper::SH_Wrapper(int _verbose){

  verbose = _verbose;

  return;
}

//
//
//
int SH_Wrapper::run(
  int K,
  std::vector< double > vecs,
  std::vector< double > scalars,
  std::vector< int > vecs_shape,
  std::vector< int > scalars_shape,
  std::vector< double >& output,
  bool only_even
  ){

  //TODO: generalize to "same vecs for all"
  if(vecs_shape.size() != 3 && scalars_shape.size() != 2){
    std::cout << "[shmpy/sh_wrap] Bad numpy dim pass. Aborting."
      << std::endl
      << "\tvecs dim is " << vecs_shape.size() << " should be 3"
      << std::endl
      << "\tscalars dim is " << scalars_shape.size() << " should be 2"
      << std::endl;
    return 1;
  }

  int num_vox = vecs_shape[0];
  int num_points = vecs_shape[1];
  int sphere_dim = vecs_shape[2];
  int mat_size = vecs_shape[1]*vecs_shape[2];

  if(sphere_dim != 3){
    std::cout << "[shmpy/sh_wrap] vecs_shape passed without 3rd dim = 3."
      " Aborting" << std::endl;
    return 1;
  }

  if(verbose > 0){
    std::cout << "Dimensions\t"
      << num_vox << "\t"
      << num_points << "\t"
      << sphere_dim << "\t"
      << std::endl;
  }

  //
  // set up outputs
  //

  output.clear();
  int num_coef = -1;
  if(only_even){
    num_coef = check_num_even_coef(K);
  } else {
    num_coef = check_num_coef(K);
  }
  output.reserve( num_coef * num_vox );
  for(int i = 0; i < num_coef * num_vox; ++i)
    output.push_back(0);

  //
  // iterate through voxels and fit
  //
#pragma omp parallel for schedule(dynamic,10) shared(output)
  for(int i = 0; i < num_vox; ++i){
    std::unique_ptr< std::vector<double> > output_local;
    std::vector<Eigen::Vector3d> samp;
    std::vector<double> sub_vec;

    ////THE WAY IT SHOULD HAVE BEEN
    //Eigen::Map<
    //  Eigen::Matrix<,Eigen::Dynamic,3, Eigen::RowMajor>
    //> M(
    //  vecs.data() + i*mat_size,
    //  num_points,
    //  3
    //);

    for( int pt_idx = 0; pt_idx < num_points; ++pt_idx )
      samp.push_back(Eigen::Vector3d(vecs.data() + i*mat_size + pt_idx*3));
    sub_vec = std::vector<double>(
      scalars.begin() + i*num_points,
      scalars.begin() + (i+1)*num_points
    );

    output_local = sh::ProjectSparseSamples(
      K, //order
      samp, //directions
      sub_vec //scalars
    );

#pragma omp critical
{
    if(only_even){
      int coef_counter = 0;
      for(int k = 0; k < K+1; k += 2){
        for(int l = 0; l < 2*k + 1; ++l){
          output[ num_coef * i + coef_counter ] = (*output_local)[k*k + l];
          coef_counter += 1;
        }
      }
    } else {

      std::copy(output_local->begin(),output_local->end(),output.begin() + num_coef*i);
      //output.insert(output.end(),output_local->begin(),output_local->end());

    }
}

    //samp.clear();
    //sub_vec.clear();
  }


  return 0;
}

//
//
//
int SH_Wrapper::run_all_same_vecs(
  int K,
  std::vector< double > vecs,
  std::vector< double > scalars,
  std::vector< int > vecs_shape,
  std::vector< int > scalars_shape,
  std::vector< double >& output,
  bool only_even
  ){

  //TODO: generalize to "same vecs for all"
  if(vecs_shape.size() != 2 && scalars_shape.size() != 2){
    std::cout << "[shmpy/sh_wrap] Bad numpy dim pass. Aborting."
      << std::endl
      << "\tvecs dim is " << vecs_shape.size() << " should be 2"
      << std::endl
      << "\tscalars dim is " << scalars_shape.size() << " should be 2"
      << std::endl;
    return 1;
  }

  int num_vox = scalars_shape[0];
  int num_points = vecs_shape[0];
  int sphere_dim = vecs_shape[1];
  int mat_size = vecs_shape[0]*vecs_shape[1];

  if(sphere_dim != 3){
    std::cout << "[shmpy/sh_wrap] vecs_shape passed without 3rd dim = 3."
      " Aborting" << std::endl;
    return 1;
  }

  if(verbose > 0){
    std::cout << "Dimensions\t"
      << num_vox << "\t"
      << num_points << "\t"
      << sphere_dim << "\t"
      << std::endl;
  }

  //
  // set up outputs
  //

  output.clear();
  if(only_even){
    output.reserve( check_num_even_coef(K) * num_vox );
  } else {
    output.reserve( check_num_coef(K) * num_vox );
  }

  //
  // iterate through voxels and fit
  //
  std::vector<Eigen::Vector3d> samp;
  std::unique_ptr< std::vector<double> > output_local;

  for( int pt_idx = 0; pt_idx < num_points; ++pt_idx )
    samp.push_back(Eigen::Vector3d(vecs.data() + pt_idx*3));

  output_local = sh::ProjectSparseSamplesSeq(
    K, //order
    samp, //directions
    scalars //scalars
  );

  int num_coefs = check_num_coef(K);

  if(only_even){
    for(int i = 0; i < num_vox; ++i){
      for(int k = 0; k < K+1; k += 2){
        for(int l = 0; l < 2*k + 1; ++l)
          output.push_back( (*output_local)[i*num_coefs + k*k + l] );
      }
    }
  } else {
    output.insert(output.end(),output_local->begin(),output_local->end());
  }

  return 0;
}


//
//
//
int SH_Wrapper::run_even(
  int K,
  std::vector< double > vecs,
  std::vector< double > scalars,
  std::vector< int > vecs_shape,
  std::vector< int > scalars_shape,
  std::vector< double >& output,
  bool all_same_vecs
  ){

  if( all_same_vecs ){
    return run_all_same_vecs(
      K, vecs, scalars, vecs_shape, scalars_shape, output, true);
  } else {
    return run(K, vecs, scalars, vecs_shape, scalars_shape, output, true);
  }

}

//
//
//
int SH_Wrapper::project_from_SH(
  int K,
  std::vector<double> sh_coefs,
  std::vector<double> vecs,
  std::vector<int> vecs_shape,
  std::vector<double>& output
  ){

  if(vecs.size() % 3 != 0){
    std::cout << "[shmpy/sh_wrap:project_from_SH] bad vec size. Aborting."
      << std::endl;
    return 1;
  }

  //
  //Prep output
  output.clear();
  output.reserve(vecs.size() / 3);

  int s = vecs.size();
  std::vector<double> temp;
  temp.reserve(3);
  Eigen::Vector3d v;
  for(int i = 0; i < s; i += 3){
    temp.insert(temp.end(), vecs.begin()+i, vecs.begin() + i + 3);
    v = Eigen::Vector3d(
      temp.data()
    );
    output.push_back( sh::EvalSHSum(K, sh_coefs, v) );
    temp.clear();
  }

  return 0;
}

//
//
//
int SH_Wrapper::project_from_even_SH(
  int K,
  std::vector<double> sh_coefs,
  std::vector<double> vecs,
  std::vector<int> vecs_shape,
  std::vector<double>& output
  ){

  std::vector<double> full_sh_coefs;
  full_sh_coefs.reserve( K*K );

  for(int k = 0; k < K+1; ++k){
    if(k % 2 == 0){
      for(int l = 0; l < 2*k + 1; ++l){
        full_sh_coefs.push_back( sh_coefs[check_num_even_coef(k-1) + l] );
      }
    } else {
      for(int l = 0; l < 2*k + 1; ++l){
        full_sh_coefs.push_back(0);        
      }
    }

  }

  return project_from_SH(
      K,
      full_sh_coefs,
      vecs,
      vecs_shape,
      output
    );
}

//
//
//
int SH_Wrapper::projection_matrix(
  int K,
  std::vector<double> vecs,
  std::vector<int> vecs_shape,
  std::vector<double>& output,
  std::vector<int>& out_shape
  ){

  int num_vecs = vecs_shape[0];
  std::vector<double> current_v;
  current_v.reserve(3);

  for(int vec_idx = 0; vec_idx < num_vecs; ++vec_idx){

    current_v.insert(
      current_v.begin(),
      vecs.begin() + vec_idx*3, vecs.begin() + (vec_idx+1)*3
    );
    Eigen::Vector3d vec(current_v.data());

    for(int k = 0; k < K; ++k){
      for(int ell = -k; ell < k+1; ++ell){
        output.push_back( sh::EvalSH(k, ell, vec) );
      }
    }

    current_v.clear();
  }

  out_shape.clear();
  out_shape.push_back(num_vecs);
  out_shape.push_back(output.size() / num_vecs);

  return 0;
}

//
//
//
int SH_Wrapper::projection_matrix_even(
  int K,
  std::vector<double> vecs,
  std::vector<int> vecs_shape,
  std::vector<double>& output,
  std::vector<int>& out_shape
  ){

  int num_vecs = vecs_shape[0];
  std::vector<double> current_v;
  current_v.reserve(3);

  for(int vec_idx = 0; vec_idx < num_vecs; ++vec_idx){

    current_v.insert(
      current_v.begin(),
      vecs.begin() + vec_idx*3, vecs.begin() + (vec_idx+1)*3
    );
    Eigen::Vector3d vec(current_v.data());

    //This needs to run up through K, so that the "harmonic order"
    //returned matches K. (There is a 0th harmonic, so K should be max_harm)
    for(int k = 0; k <= K; ++k){
      if(k % 2 == 1)
        continue;

      for(int ell = -k; ell < k+1; ++ell){
        output.push_back( sh::EvalSH(k, ell, vec) );
      }
    }

    current_v.clear();
  }

  out_shape.clear();
  out_shape.push_back(num_vecs);
  out_shape.push_back(output.size() / num_vecs);

  return 0;
}

//
//  Exposes sh:: coef counter to python
//
int SH_Wrapper::check_num_coef(int K){
  if(K < 0)
    return -1;
  return sh::GetCoefficientCount(K);
}

int SH_Wrapper::check_num_even_coef(int K){
  int num_even_coef = 0;
  for(int k = 0; k < K+1; k += 2)
    num_even_coef += 2*k + 1;
  return num_even_coef;
}






