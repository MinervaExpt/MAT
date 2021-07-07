#An alternative build system for PlotUtils.
#Inspired by Heidi Schellman's CMakeLists.txt file in build.
#Debugged and maintained by Andrew Olivier aolivier@ur.rochester.edu
cmake_minimum_required(VERSION 2.8.12)
project(PlotUtils CXX C)

# Include module with function 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

#Set up different build modes with custom compiler flags.
#I could do this, but a c++17-dependent ROOT build overrides it anyway.
#set( CMAKE_CXX_STANDARD 11 )
#set( CMAKE_CXX_STANDARD_REQUIRED ON )

#I think Heidi disabled warnings about the form() function with -Wformat=0
set( GCC_Flags_For_CXX "-Wall -fPIC -pthread -Wformat=0 --std=c++0x" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_Flags_For_CXX}" )
set( CMAKE_CXX_FLAGS_DEBUG "-ggdb" )

#Required packages
find_package(ROOT REQUIRED COMPONENTS RIO Net OPTIONAL_COMPONENTS Reflex Cintex)

if(${ROOT_VERSION} VERSION_LESS 6 AND NOT ${ROOT_Reflex_FOUND})
  MESSAGE(FATAL_ERROR "Reflex is optional except when it's not.  ROOT 6 has Reflex "
                      "support built in, so it doesn't have a separate component "
                      "for Reflex.  Reflex was an experimental feature in ROOT 5, "
                      "so I have to require it as a component there.  You appear to "
                      "be using ${ROOT_VERSION}, so I can't let you get away with "
                      "skipping out on Reflex support!  I need Reflex to build the "
                      "ROOT dictionary for PlotUtils.")
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} $ENV{BOOSTDIR})
#Hack to dodge a bug I'm encountering with ROOT 5.34 on the MINERvA GPVMs but
#not develop of ROOT (6) from mid 2019.  ROOT is failing to set up a variable
#to find genreflex for some weird reason.  It's on PATH when I have the MINERvA
#framework setup, so I'll just go find it myself.
if(NOT ROOT_genreflex_cmd)
  MESSAGE("Hmm, ROOT failed to set up genreflex for some weird reason.  I'll "
          "go off and try to find it myself...")
  find_program(ROOT_genreflex_cmd genreflex)
  if(NOT ROOT_genreflex_cmd)
    MESSAGE(FATAL_ERROR "Failed to find genreflex on PATH!  ROOT is supposed to "
                        "set it up but doesn't seem to on the MINERvA GPVMs with "
                        "ROOT 6.  I need genreflex to build the ROOT dictionary "
                        "for PlotUtils.")
  endif()
endif()

#Ignore Cintex extensions when they've been superseded by ROOT 6
if(${ROOT_VERSION} VERSION_LESS 6)
  MESSAGE("I need Cintex on this platform, so making sure to compile against it...")
else()
  MESSAGE("I don't need Cintex on this platform, so skipping it...")
  add_definitions(-DNCINTEX)
endif()

#From Heidi: ---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#All components in this project can refer to other components by
#relative path from the root source code directory.
include_directories( ${CMAKE_SOURCE_DIR} )

#From Heidi: put some special flags here forgot what the DBUILD_SHARED_LIBS=OFF does..
add_definitions(-DFORM -DMNVROOT6 -DPLOTUTILS_STANDALONE "-DBUILD_SHARED_LIBS=OFF")

#From Heidi: do this or mac will make it dylib
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")


MESSAGE("boost:$ENV{BOOSTDIR}")
include_directories(PlotUtils $ENV{PLOTUTILSROOT} $ENV{BOOSTDIR}/include $ENV{BOOSTDIR})

add_subdirectory(PlotUtils)
add_subdirectory(data)
add_subdirectory(macros)

#Install setup script.  This adds the install location of PlotUtils to
#LD_LIBRARY_PATH and sets other useful environment variables.
configure_file(setup.sh.in setup_${PROJECT_NAME}.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setup_${PROJECT_NAME}.sh DESTINATION bin)

#Make the results of this build into a distributable package.  Designed to be distributed as a .tar.gz
#Learned to do this from http://agateau.com/2009/cmake-and-make-dist/
set( CPACK_PACKAGE_VERSION_MAJOR "0" )
set( CPACK_PACKAGE_VERSION_MINOR "0" )
set( CPACK_PACKAGE_VERSION_PATCH "0" )
set( CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
set( CPACK_SOURCE_GENERATOR "TGZ" )

include( CPack )
add_custom_target( dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source )

#Learned to set up CMake package from Clark McGrew's edep-sim package
# Write the 'PlotUtilsConfigVersion.cmake' file which can be used to
# check if a version meets the requested properties.
write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  COMPATIBILITY SameMajorVersion
  VERSION 1.0)

# Write the 'EDepSimConfig.cmake' file so that a user package can
# access this with find_package.
configure_package_config_file(
  PackageConfig.cmake.in
  ${PROJECT_NAME}Config.cmake
  PATH_VARS CMAKE_INSTALL_PREFIX
  INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})

# Install the config files.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})
