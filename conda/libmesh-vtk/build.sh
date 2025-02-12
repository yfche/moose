#!/bin/bash
set -eu
export PATH=/bin:$PATH

if [[ $mpi == "openmpi" ]]; then
  export OMPI_MCA_plm=isolated
  export OMPI_MCA_rmaps_base_oversubscribe=yes
  export OMPI_MCA_btl_vader_single_copy_mechanism=none
elif [[ $mpi == "moose-mpich" ]]; then
  export HYDRA_LAUNCHER=fork
fi
export CC=mpicc CXX=mpicxx
mkdir -p build
cd build
VTK_PREFIX=${PREFIX}/libmesh-vtk
cmake .. -G "Ninja" \
    -Wno-dev \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH:PATH=${VTK_PREFIX} \
    -DCMAKE_INSTALL_PREFIX:PATH=${VTK_PREFIX} \
    -DCMAKE_INSTALL_RPATH:PATH=${VTK_PREFIX}/lib \
    -DVTK_BUILD_DOCUMENTATION:BOOL=OFF \
    -DVTK_BUILD_TESTING:BOOL=OFF \
    -DVTK_BUILD_EXAMPLES:BOOL=OFF \
    -DBUILD_SHARED_LIBS:BOOL=ON \
    -DVTK_USE_MPI:BOOL=ON \
    -DVTK_GROUP_ENABLE_Rendering:STRING=DONT_WANT \
    -DVTK_GROUP_ENABLE_Qt::STRING=NO \
    -DVTK_GROUP_ENABLE_Views:STRING=NO \
    -DVTK_GROUP_ENABLE_Web:STRING=NO \
    -DCMAKE_OSX_SYSROOT=${CONDA_BUILD_SYSROOT}

CORES=$(echo "${CPU_COUNT:-2} / 2" | bc)
ninja install -v -j $CORES

# VTK 9.1 now places libs in lib64 when installed on linux, linking to the "expected" location of lib
if [[ $(uname) == Linux ]]; then
    ln -s ${VTK_PREFIX}/lib64 ${VTK_PREFIX}/lib
fi

# Set LIBMESH_DIR environment variable for those that need it
mkdir -p "${PREFIX}/etc/conda/activate.d" "${PREFIX}/etc/conda/deactivate.d"
cat <<EOF > "${PREFIX}/etc/conda/activate.d/activate_${PKG_NAME}.sh"
export VTKLIB_DIR=${VTK_PREFIX}/lib
export VTKINCLUDE_DIR=${VTK_PREFIX}/include/vtk-${SHORT_VTK_NAME}
EOF
cat <<EOF > "${PREFIX}/etc/conda/deactivate.d/deactivate_${PKG_NAME}.sh"
unset VTKLIB_DIR
unset VTKINCLUDE_DIR
EOF
