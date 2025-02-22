#!/bin/bash

#
# This is a helper script for installing the AMD clang-based compiler
# needed to compile a GPU-enabled version of unifrac.
#
# Note: The script currently assumes Linux_x86_64 platform.
#

# Note:
# for the time being, it reuqires sudo access
# as we do not support static linking, and thus must have abspaths

#
# based on https://github.com/ROCm/aomp/blob/aomp-dev/docs/INSTALL.md
#

export AOMP_VERSION=20.0-1

export AOMP_BASEDIR=/usr/lib

if [ -f ${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/llvm/bin/amdclang++ ]; then
 echo "[INFO] Found existing ${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/llvm/bin/amdclang++"
else
 dist=`cat /etc/os-release |awk '/^NAME=/{split($1,a,"="); split(a[2],b,"\""); print b[2]}'`
 osver=`cat /etc/os-release |awk '/^VERSION_ID=/{split($1,a,"="); split(a[2],b,"\""); print b[2]}'`
 if [ "x$dist" = "xUbuntu" ]; then
  if [ "x$osver" == "x22.04" ]; then
    echo "[INFO] Fetching https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_${AOMP_VERSION}/aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb"
    curl -s -L https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_${AOMP_VERSION}/aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb -o aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb
    if [ $? -ne 0 ]; then
      echo "Failed to download aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb"
      exit 1
    fi
    echo "[INFO] Executing sudo dpkg -i aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb"
    echo "Please insert password, if you want to install aomp in /usr/lib"
    sudo apt install ./aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb
    if [ $? -ne 0 ]; then
      echo "[WARNING] System install failed, doing unprivileged install in $PWD/aomp_clang"
      mkdir $PWD/aomp_clang
      dpkg -x aomp_Ubuntu2204_${AOMP_VERSION}_amd64.deb $PWD/aomp_clang
      AOMP_BASEDIR=$PWD/aomp_clang/usr/lib
    fi
    rm -f aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb
  elif [ "x$osver" == "x24.04" ]; then
    echo "[INFO] Fetching https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_${AOMP_VERSION}/aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb"
    curl -s -L https://github.com/ROCm-Developer-Tools/aomp/releases/download/rel_${AOMP_VERSION}/aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb -o  aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb
    if [ $? -ne 0 ]; then
      echo "Failed to download aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb"
      exit 1
    fi
    echo "[INFO] Executing sudo dpkg -i aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb"
    echo "Please insert password, if you want to install aomp in /usr/lib"
    sudo apt install ./aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb
    if [ $? -ne 0 ]; then
      echo "[WARNING] System install failed, doing unprivileged install in $PWD/aomp_clang"
      mkdir $PWD/aomp_clang
      dpkg -x aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb $PWD/aomp_clang
      AOMP_BASEDIR=$PWD/aomp_clang/usr/lib
    fi
    rm -f aomp_Ubuntu2404_${AOMP_VERSION}_amd64.deb
  else
   echo "Unsupported $dist Linux version $osver"
   exit 1
  fi
 else
  echo "Unsupported Linux distribution $dist"
  exit 1
 fi
fi


if [ -f ${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/llvm/bin/amdclang++ ]; then

cat > setup_amd_compiler.sh  << EOF
export PATH=${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/:${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/bin/:${AOMP_BASEDIR}/aomp_${AOMP_VERSION}/llvm/bin:\$PATH

export AMD_CXX=amdclang++

# no special  flags needed in most cases
export AMD_CPPFLAGS=
export AMD_CXXFLAGS=
export AMD_CFLAGS=

export AMD_LDFLAGS=
EOF

echo "Setup script avaiabile in $PWD/setup_amd_compiler.sh"

else
  # something went wrong with the installation process above
  echo "Failed to install aomp_${AOMP_VERSION}"
  exit 1
fi

