How to build PlotUtils using CMake:

Quick start in bash:
1. cvs co Ana/PlotUtils #You need a kerberos ticket if you're not on a GPVM
2. mkdir PlotUtilsInstallPrefix && cd PlotUtilsInstallPrefix #An out of source place for the binaries to go
3. mkdir debug && cd debug #We're going to do a build with debug symbols
4. mkdir build && cd build #The Makefiles and raw .o files will go here
5. cmake ../../../Ana/PlotUtils -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Debug
6. make install #Puts files in PlotUtilsInstallPrefix/debug/*
7. source /path/to/PlotUtilsInstallPrefix/debug/bin/setup_PlotUtils.sh
8. echo $PLOTUTILSROOT && echo $PATH && ls $PLOTUTILSROOT #Check that env is correct and libraries are there

Use with other CMake Packages:
- In your other package's CMakeLists.txt (usually the top level):
  find_package(PlotUtils REQUIRED)
- When "generating" Makefiles for your package for the first time (like 5. above), you might need to:
  export PlotUtils_DIR=/path/to/PlotUtilsInstallPrefix

Expert details:
- Nota Bene: Requires CMake 2.8.12.  This is the one you get from the MINERvA framework on the GPVMs.
- An optimized build without debug symbols.  Substitute -DCMAKE_BUILD_TYPE=Release in step 5.
- Portability for e.g. the grid:
  export MINERVA_PREFIX=$CONDOR_SOME_DIR/PlotUtilsInstallPrefix
  source $CONDOR_SOME_DIR/PlotUtilsInstallPrefix/bin/setup_PlotUtils.sh
- Works great with CMake's ExternalProject_Add().  My favorite so far looks like this:
  ExternalProject_Add(PlotUtils
                      CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                      CVS_MODULE Ana/PlotUtils
                      SOURCE_DIR ${CMAKE_SOURCE_DIR}/Ana/PlotUtils
                      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

- ExternalProject_Add() is a fickle beast though.  It has a bad reputation online because people usually
  put it directly in a project that uses e.g. PlotUtils.  Don't do that.  Instead, I create an "uber build"
  package that orchestrates my package and all of its dependencies.  This is what lcgcmake actually is!
  - You just need find_package() like above
  - Write another CMakeLists.txt for the uber-build that goes roughly like this:
    include(ExternalProject)
    ExternalProject_Add(PlotUtils
                    CVS_REPOSITORY minervacvs@cdcvs.fnal.gov:/cvs/mnvsoft
                    CVS_MODULE Ana/PlotUtils
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/Ana/PlotUtils
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})
    ExternalProject_Add(NucCCNeutrons
                    GIT_REPOSITORY https://github.com/MinervaExpt/NucCCNeutrons.git
                    GIT_TAG develop
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/NucCCNeutrons
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} -Dyaml-cpp_DIR:STRING=${CMAKE_INSTALL_PREFIX}/lib/yaml-cpp/cmake/yaml-cpp
                    UPDATE_COMMAND ""
                    PATCH_COMMAND ""
                    DEPENDS PlotUtils)

  Notice the DEPENDS on the last line.  It seems to set up PlotUtils_DIR for you!
