#!/bin/bash
# Install SuiteSparse
# Call with
# ./install.sh /path/to/install/dir

# Exit the script at the first failure
set -e

INSTALL_DIR="$1"
# Number of procs for building (default 4)
NPROCS=${NPROCS=4}
VERSION="5.4.0"
CHECKSUM="d9d62d539410d66550d0b795503a556830831f50087723cb191a030525eda770"

wget --no-verbose https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${VERSION}.tar.gz
# Verify checksum
if [ $CHECKSUM = `sha256sum v${VERSION}.tar.gz | awk '{print $1}'` ]
then
  echo "Checksum matches"
else
  echo "Checksum does not match"
  exit 1
fi

tar -xzf v${VERSION}.tar.gz
cd SuiteSparse-${VERSION}/
make -j${NPROCS} library BLAS=-lblas
make install INSTALL=${INSTALL_DIR} BLAS=-lblas
cd ../
rm -rf SuiteSparse*
# Need to specify metis location?
