#### Falicov-Kimball MC
##### Dependencies 
- c++11-compatible compiler
- triqs 1.0.0 release : https://github.com/TRIQS/triqs/tree/1.0.0
- all triqs deps : boost, cmake, hdf5 (with c++ bindings compiled in c11 mode), mpi
- eigen (with unsupported arpack module) : `hg clone https://bitbucket.org/eigen/eigen/`
- TCLAP command line parser : http://tclap.sourceforge.net/
- gtest : https://code.google.com/p/googletest/ (fetched automatically) 

##### Sample cmake compilation script : 
```
cmake \
-DCMAKE_INSTALL_PREFIX="$SOME_PATH" \
-DCMAKE_CXX_FLAGS="-ftemplate-depth=256" \
-DEIGEN3_INCLUDE_DIR=${HOME}/code/eigen3 \
-DTesting=ON \
-DBenchmark=OFF \
-DCMAKE_BUILD_TYPE="Release" \
-DLATTICES="triangular" \
..
```
