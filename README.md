## Falicov-Kimball MC
### Dependencies 
- c++11-compatible compiler
- ALPSCore (http://alpscore.org)
- ALPSCore requirements : cmake, hdf5, mpi.
- Boost: boost-mpi
- arpack (arpack-ng)
- eigen `http://eigen.tuxfamily.org/`
- gtest : https://github.com/google/googletest (fetched automatically) 

### Installation 
1. Get all dependencies except of triqs 
2. Download fk_mc: `git clone https://github.com/aeantipov/fk_mc.git`. The installation is similar to TRIQS. 
    1. Create a temporary build directory
    2. In it, run 
      ```
        cmake \
        -DCMAKE_INSTALL_PREFIX=FK_MC_INSTALL_DIR" \
        -DCMAKE_CXX_FLAGS="-ftemplate-depth=256" \
        -DTesting=ON \
        -DBenchmark=ON \
        -DCMAKE_BUILD_TYPE="Release" \
        -DALPSCore_DIR=path_to_ALPSCoreConfig.cmake_dir
        FK_MC_DOWNLOAD_DIR
      ```
      
      where FK_MC_DOWNLOAD_DIR is the location if the downloaded fk_mc source (with CMakeLists.txt), TRIQS_INSTALL_DIR is a location where you want to install fk_mc. 
    3. (optional for configuration) specify -DLATTICES="cubic2d;cubic3d;triangular" etc to compile the code for different lattice types and dimensions. By default all lattices are compiled.
    4. run `make`
    5. (optional) run `make test` and verify that the tests pass. This typically saves for most encountered errors.
    6. run `make install`

#### Running fk_mc
In FK_MC_INSTALL_DIR/bin you will find fk_mc_LATTICE executables, where LATTICE corresponds to the chosen lattice type such as cubic2d. 
The standard way to run is to do 
```
  mpirun fk_mc_cubic2d --beta 1.0 --U 1.0 --L 8
```

Here beta is the inverse temperature, U is the interaction strength, and L is the linear system size (volume = L^d). Run `fk_mc_cubic2d --help` to learn about all the parameters. 

The output is stored in the hdf5 archive (by default - output.h5). Specifying `--plaintext` option in the run dumps the value of observables to plaintext files. For example specific heat is obtained then in `cv_error.dat` and reads:
`5.000000e+01 -2.474323e-01 1.207943e-28 1.554312e-15`. First column is the number of samples used, second - the value of the observable, third is variance, and last is the standard error. The binning is done to avoid autocorrelations, so the number of samples is typically smaller than the number of measurements.

#### Authors & Contributors
- Andrey Antipov, *Andrey.E.Antipov[at]gmail.com*, 2013-now
- Andreas Herrmann, University of Fribourg

#### Using 
- please cite http://arxiv.org/abs/1605.01390 [Phys. Rev. Lett. 117, 146601 (2016)]
