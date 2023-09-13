autoreconf -i
./configure --prefix=$PREFIX \
            --with-mpi=$BUILD_PREFIX \
            --with-hdf5=$BUILD_PREFIX \
            CFLAGS="-O2 -DNDEBUG" \
            CXXFLAGS="-O2 -DNDEBUG" \
            LIBS="-ldl -lz"
make install
