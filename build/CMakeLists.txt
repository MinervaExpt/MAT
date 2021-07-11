# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.

cmake_minimum_required(VERSION 3.1...3.14)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables

project(plotutils VERSION 1.0 LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} $ENV{BOOSTDIR})

set(CMAKE_CXX_FLAGS "-Wno-format -Wno-unused-value" )  # get rid of errors from Form and assert

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
link_directories(${ROOT_LIBRARY_DIR})

MESSAGE("ROOT_LIBRARY_DIR:${ROOT_LIBRARY_DIR}")

# put some special flags here forgot what the DBUILD_SHARED_LIBS=OFF does..

add_definitions( -DFORM -DMNVROOT6 -DPLOTUTILS_STANDALONE "-DBUILD_SHARED_LIBS=OFF")

# do this or mac will make it dylib
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")


# tell it about the source files
set(SRC_FILES
PlotUtils/AnaBinning.cxx
PlotUtils/ArachneUtils.cxx
PlotUtils/BaseUniverse.cxx
PlotUtils/ChainWrapper.cxx
PlotUtils/Cutter.cxx
PlotUtils/DefaultCVUniverse.cxx
PlotUtils/FluxReweighter.cxx
PlotUtils/FluxReweighterWithWiggleFit.cxx
PlotUtils/FluxSystematics.cxx
PlotUtils/GeantHadronSystematics.cxx
PlotUtils/GridCanvas.cxx
PlotUtils/genie_particle.cxx
PlotUtils/GenieSystematics.cxx
PlotUtils/AngleSystematics.cxx
PlotUtils/Hist2DWrapper.cxx
PlotUtils/HistWrapper.cxx
PlotUtils/HistogramUtils.cxx
PlotUtils/HyperDimLinearizer.cxx
PlotUtils/PhysicsVariables.cxx
PlotUtils/MacroUtil.cxx
PlotUtils/MinosEfficiencySystematics.cxx
PlotUtils/MinosMuonEfficiencyCorrection.cxx
#PlotUtils/MinosMuonPlusEfficiencyCorrection.cxx
PlotUtils/MnvApplication.cxx
PlotUtils/MnvColors.cxx
PlotUtils/MnvFluxConstraint.cxx
PlotUtils/MnvH1D.cxx
PlotUtils/MnvH2D.cxx
PlotUtils/MnvHadronReweight.cxx
PlotUtils/MnvLatErrorBand.cxx
PlotUtils/MnvLatErrorBand2D.cxx
PlotUtils/MnvNormalization.cxx
PlotUtils/MnvNuclearModelWeight.cxx
PlotUtils/MnvPlotter.cxx
PlotUtils/MnvRecoShifter.cxx
PlotUtils/MnvTuneSystematics.cxx
PlotUtils/MnvVertErrorBand.cxx
PlotUtils/MnvVertErrorBand2D.cxx
PlotUtils/MinervaUniverse.cxx
PlotUtils/MuonSystematics.cxx
PlotUtils/MuonResolutionSystematics.cxx
PlotUtils/NuclModUtils.cxx
PlotUtils/ParticleResponseDefaults.cxx
PlotUtils/POTCounter.cxx
PlotUtils/ResponseSystematics.cxx
PlotUtils/ROOTglob.cxx
# PlotUtils/Table.cxx  # this is included in the header
PlotUtils/TargetUtils.cxx
PlotUtils/TreeWrapper.cxx
PlotUtils/TrueProtonKECutSystematics.cxx
PlotUtils/RecoProtonKECutSystematics.cxx
PlotUtils/makeChainWrapper.cxx
PlotUtils/weightDIS.cxx
PlotUtils/weightCoherentPi.cxx
PlotUtils/weightGenieBodekRitchieClass.cxx
PlotUtils/weightLowQ2Pi.cxx
PlotUtils/weightMK.cxx
PlotUtils/weightRPA.cxx
PlotUtils/weightRemoveUnphysical2p2hExtendedEventsClass.cxx
PlotUtils/weightSusaGenieQEClass.cxx
PlotUtils/weightSusaValenciaClass.cxx
PlotUtils/weightZExp.cxx
PlotUtils/weight_2p2h.cxx
PlotUtils/weight_fsi.cxx
PlotUtils/weight_fsi_absorption.cxx
PlotUtils/weight_fsi_cai.cxx
#PlotUtils/MnvH1DToCSV.cxx
#PlotUtils/MnvH2DToCSV.cxx
)

# tell it about the header files

set(H_FILES
PlotUtils/AnaBinning.h
PlotUtils/AngleSystematics.h
PlotUtils/ArachneUtils.h
PlotUtils/BaseUniverse.h
PlotUtils/ChainWrapper.h
PlotUtils/Cut.h
PlotUtils/Cutter.h
PlotUtils/DefaultCVUniverse.h
PlotUtils/Exceptions.h
PlotUtils/FluxReweighter.h
PlotUtils/FluxReweighterWithWiggleFit.h
PlotUtils/FluxSystematics.h
PlotUtils/GeantHadronSystematics.h
PlotUtils/GridCanvas.h
PlotUtils/genie_particle.h
PlotUtils/GenieSystematics.h
PlotUtils/Hist2DWrapper.h
PlotUtils/HistWrapper.h
PlotUtils/HistogramUtils.h
PlotUtils/HyperDimLinearizer.h
PlotUtils/KillCintex.h
PlotUtils/PhysicsVariables.h
PlotUtils/MacroUtil.h
PlotUtils/MinosEfficiencySystematics.h
PlotUtils/MinosMuonEfficiencyCorrection.h
#PlotUtils/MinosMuonPlusEfficiencyCorrection.h
PlotUtils/MinervaUniverse.h
PlotUtils/MnvAnaTuple.h
PlotUtils/MnvApplication.h
PlotUtils/MnvColors.h
PlotUtils/MnvFluxConstraint.h
PlotUtils/MnvH1D.h
PlotUtils/MnvH2D.h
PlotUtils/MnvHadronReweight.h
PlotUtils/MnvLatErrorBand.h
PlotUtils/MnvLatErrorBand2D.h
PlotUtils/MnvNormalization.h
#PlotUtils/MnvNormalizerME.h
PlotUtils/MnvNuclearModelWeight.h
PlotUtils/MnvPlotter.h
PlotUtils/MnvRecoShifter.h
PlotUtils/MnvTuneSystematics.h
PlotUtils/MnvVertErrorBand.h
PlotUtils/MnvVertErrorBand2D.h
PlotUtils/MuonSystematics.h
PlotUtils/MuonResolutionSystematics.h
PlotUtils/NuclModUtils.h
PlotUtils/ParticleResponseDefaults.h
PlotUtils/POTCounter.h
PlotUtils/PlotUtilsPhysicalConstants.h
PlotUtils/ROOTglob.h
PlotUtils/ResponseSystematics.h
PlotUtils/Table.h
PlotUtils/TargetUtils.h
PlotUtils/TreeWrapper.h
PlotUtils/TrueProtonKECutSystematics.h
PlotUtils/RecoProtonKECutSystematics.h
PlotUtils/makeChainWrapper.h
PlotUtils/weightDIS.h
PlotUtils/weightGenieBodekRitchieClass.h
PlotUtils/weightLowQ2Pi.h
PlotUtils/weightMK.h
PlotUtils/weightRPA.h
PlotUtils/weightRemoveUnphysical2p2hExtendedEventsClass.h
PlotUtils/weightSusaGenieQEClass.h
PlotUtils/weightSusaValenciaClass.h
PlotUtils/weightZExp.h
PlotUtils/weight_2p2h.h
PlotUtils/weightCoherentPi.h
PlotUtils/weight_fsi.h
PlotUtils/weight_fsi_absorption.h
PlotUtils/weight_fsi_cai.h
PlotUtils/weight_minervaq2qe.h
)

# need to know about Boost on macs so check that you know where it is.

MESSAGE("boost:$ENV{BOOSTDIR}")
include_directories(PlotUtils $ENV{PLOTUTILSROOT} $ENV{BOOSTDIR}/include $ENV{BOOSTDIR})

# now generate the dictionary

ROOT_GENERATE_DICTIONARY(G__${PROJECT_NAME} ${H_FILES} LINKDEF ${CMAKE_SOURCE_DIR}/../PlotUtils/Linkdef.h)

# and then the shared library
# Output libname matches target name, with the usual extensions on your system

add_library(${PROJECT_NAME} SHARED ${SRC_FILES} G__${PROJECT_NAME}.cxx )
#add_library(${PROJECT_NAME} SHARED ${SRC_FILES} ${PROJECT_NAME}.cxx )

MESSAGE("ROOT_LIBRARIES:${ROOT_LIBRARIES}")

MESSAGE("PROJECT:${PROJECT_NAME}")

#target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} MathMore)
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES})

# Link each target with other targets or add options, etc.

# Adding something we can run - Output name matches target name
#add_executable(TransWarpExtraction macros/TransWarpExtraction.cxx)

# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
#target_link_libraries(TransWarpExtraction PRIVATE ${PROJECT_NAME} )
