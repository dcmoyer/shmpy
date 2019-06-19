

from libcpp.vector cimport vector
from libcpp cimport bool
import numpy as np
cimport numpy as np

##
## import in the c++ func
##

cdef extern from "sh_wrap.hpp":
  cdef cppclass SH_Wrapper:
    SH_Wrapper(int)

    int run(
      int K,
      vector[double] vecs,
      vector[double] scalars,
      vector[int] vecs_shape,
      vector[int] scalars_shape,
      vector[double]& outputs
    )

    int run_even(
      int K,
      vector[double] vecs,
      vector[double] scalars,
      vector[int] vecs_shape,
      vector[int] scalars_shape,
      vector[double]& outputs,
      bool all_same
    )

    int check_num_coef(
      int K
    )

    int check_num_even_coef(
      int K
    )

    int project_from_SH(
      int K,
      vector[double] sh_coefs,
      vector[double] vecs,
      vector[int] vecs_shape,
      vector[double]& output
    )

    int project_from_even_SH(
      int K,
      vector[double] sh_coefs,
      vector[double] vecs,
      vector[int] vecs_shape,
      vector[double]& output
    )

    int projection_matrix(
      int K,
      vector[double] vecs,
      vector[int] vecs_shape,
      vector[double]& output,
      vector[int]& out_shape
    )

    int projection_matrix_even(
      int K,
      vector[double] vecs,
      vector[int] vecs_shape,
      vector[double]& output,
      vector[int]& out_shape
    )

##
##
##

##This, unforunately, will be where we handle memory alloc
##
cdef class Py_SH_Wrapper:

  #this is the c++ instance
  cdef SH_Wrapper* wrapper_ptr 

  def __cinit__(self,verbose=0):
    self.wrapper_ptr = new SH_Wrapper(<int>verbose)

  def __dealloc__(self):
    del self.wrapper_ptr

  #THIS IS WHERE THE REAL STUFF WILL GO DOWN
  def fit_sh(self, K, np.ndarray vecs, np.ndarray scalars,
    vector[int] vecs_shape, vector[int] scalars_shape):
    print("Hello from the Python Wrapper.")

    cdef vector[double] output

    self.wrapper_ptr.run(K, vecs, scalars, vecs_shape, scalars_shape, output)

    return output

  def check_num_coef(self, int K):
    n_coef = self.wrapper_ptr.check_num_coef(K)
    print("[shmpy/sh_wrap] K=%i n_coef=%i" % (K, n_coef))
    return n_coef

  def check_num_even_coef(self, int K):
    #n_coef = 0
    #for k in range(0,K+1,2):
    #  n_coef += 2*k + 1
    n_coef = self.wrapper_ptr.check_num_even_coef(K)
    print("[shmpy/sh_wrap] K=%i n_even_coef=%i" % (K, n_coef))
    return n_coef

  #secretly, this fits all coefs and then just outputs only evens
  def fit_even_sh(self, K, np.ndarray vecs, np.ndarray scalars,
    vector[int] vecs_shape, vector[int] scalars_shape, all_same=False):
    print("Hello from the Python Wrapper.")
    cdef vector[double] output
    self.wrapper_ptr.run_even(\
      K, vecs, scalars, vecs_shape, scalars_shape, output, all_same)
    return output

  def project_from_sh(self,
      int K,
      vector[double] sh_coefs,
      np.ndarray vecs,
      vector[int] vecs_shape
    ):
    cdef vector[double] output
    self.wrapper_ptr.project_from_SH(K, sh_coefs, vecs, vecs_shape, output)
    return output

  def project_from_even_sh(self,
      int K,
      vector[double] sh_coefs,
      np.ndarray vecs,
      vector[int] vecs_shape
    ):
    cdef vector[double] output
    self.wrapper_ptr.project_from_even_SH(K, sh_coefs, vecs, vecs_shape, output)
    return output

  def projection_matrix(self,
      int K,
      np.ndarray vecs,
      vector[int] vecs_shape
    ):
    cdef vector[double] output
    cdef vector[int] output_shape
    self.wrapper_ptr.projection_matrix(K, vecs, vecs_shape, output, output_shape)
    return (output, output_shape)

  def projection_matrix_even(self,
      int K,
      np.ndarray vecs,
      vector[int] vecs_shape
    ):
    cdef vector[double] output
    cdef vector[int] output_shape
    self.wrapper_ptr.projection_matrix_even(\
      K, vecs, vecs_shape, output, output_shape)
    return (output, output_shape)

  #this is the "printed" representation
  def __repr__(self):
    return("You printed a Py_SH_Wrapper! " +\
      "(idk why you'd do that but here we are.)")





