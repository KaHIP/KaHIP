#!/bin/bash

#******************************************************************************
#* KaHIPWrapper.cpp
#*
#* Example wrapper for Java integration of KaHIP via JNI
#*
#******************************************************************************
#* Copyright (C) 2014 Uniserv GmbH
#*
#* This program is free software: you can redistribute it and/or modify it
#* under the terms of the GNU General Public License as published by the Free
#* Software Foundation, either version 2 of the License, or (at your option)
#* any later version.
#*
#* This program is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
#* more details.
#*
#* You should have received a copy of the GNU General Public License along with
#* this program.  If not, see <http://www.gnu.org/licenses/>.
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
