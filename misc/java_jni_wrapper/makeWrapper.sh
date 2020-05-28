#!/bin/bash

#******************************************************************************
#* KaHIPWrapper.cpp
#*
#* Example wrapper for Java integration of KaHIP via JNI
#*
#*****************************************************************************/

echo "==================================================================" 
echo "run ./compile.sh in main project folder before running this script" 
echo "==================================================================" 

export JAVA_HOME=/usr/lib/jvm/default-java
export KAHIP_HOME=../../deploy/
export OPENMPI_HOME=/usr/lib64/openmpi
export INCLUDEDIRS="-I$JAVA_HOME/include -I$JAVA_HOME/include/linux -I$KAHIP_HOME -I$OPENMPI_HOME/include"
export LIBDIRS="-L$KAHIP_HOME -L$OPENMPI_HOME/lib "

set -x
javac KaHIPWrapper.java KaHIPWrapperResult.java
javah KaHIPWrapper
g++ $INCLUDEDIRS $LIBDIRS -fPIC -c KaHIPWrapper.cpp  
mpicxx $INCLUDEDIRS $LIBDIRS -shared -fPIC -o libwrapkahip.so -Wl,-soname,wrapkahip KaHIPWrapper.o -lkahip -lmpi -lmpi_cxx -lgomp
set +x
