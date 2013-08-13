#  Copyright Olivier Parcollet 2010, Andrey Antipov 2013.
#
# This module looks for tclap command option parser.
# It sets up : ARPACK_LIBRARIES
# 

find_package(PkgConfig)
pkg_check_modules(PC_ARPACK QUIET tclap)

SET(TRIAL_PATHS
 $ENV{ARPACK_ROOT}/lib
 ${ARPACK_ROOT}/lib
 /usr/lib
 /usr/local/lib
 /opt/local/lib
 /sw/lib
 )

find_library(ARPACK_LIBRARY NAMES arpack libarpack
             HINTS ${PC_ARPACK_LIBRARY_DIRS} ${TRIAL_PATHS})

set(ARPACK_LIBRARIES ${ARPACK_LIBRARY} )


find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARY)

mark_as_advanced(ARPACK_LIBRARIES)

