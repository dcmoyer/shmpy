
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np
 
setup(
  name = 'SH_Wrap',
  ext_modules=[ 
    Extension("shmpy", 
      # Note, you can link against a c++ lib instead of including the source
      sources=[\
        "shmpy/shmpy.pyx",\
        "shmpy/sh_wrap.cpp"
      ],
      include_dirs=[np.get_include(),"spherical_harmonics/"],
      library_dirs=["spherical_harmonics/bin/"],
      libraries=["spherical_harmonics","gomp"],
      language="c++",
      extra_compile_args=["-fopenmp",'-std=c++11']),
    ],
  cmdclass = {'build_ext': build_ext},
)








