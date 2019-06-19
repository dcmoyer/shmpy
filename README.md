# shmpy

This is a basic Cython wrapper of the (unofficial)
[Google Spherical Harmonics Library](https://github.com/google/spherical-harmonics).

It requires:
1. [OpenMP](https://gcc.gnu.org/onlinedocs/libgomp/), currently hardcoded for `gomp`.
2. [Eigen](http://eigen.tuxfamily.org)
3. A compiler that supports C++11 or newer, and Make.

It also requires the SH library, which I've helpfully included.

It's not really made for human consumption yet, and is unstable.

Use at your own risk.

## Compilation

Run the following in the root directory of the repo:

    make

The compiled library will be at `shmpy/shmpy.*.so` where the wildcard is the
compiling python version.

