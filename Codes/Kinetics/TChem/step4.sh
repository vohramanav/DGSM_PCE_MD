./configure CC=/usr/bin/gcc                   \
                   CXX=/usr/bin/g++                  \
                   F77=/usr/bin/gfortran             \
                   --with-examples=yes                                      \
                   --with-sundials-dir=$PWD/lib/cvode-2.7.0-install         \
                   --with-dvode-dir=$PWD/lib                                \
                   --with-boost-dir=/usr/lib/x86_64-linux-gnu

