add_executable(madd madd.cxx)
target_link_libraries(madd PlotUtils ${ROOT_LIBRARIES})
install(TARGETS madd DESTINATION bin)

#add_executable(TransWarpExtraction TransWarpExtraction.cxx)
#target_link_libraries(TransWarpExtraction ${ROOT_LIBRARIES} PlotUtils ${PROJECT_NAME})
#install(TARGETS TransWarpExtraction DESTINATION bin)

install(FILES ProcessMCSampleSizeScan.py PlotMCSampleSizeScan.py DESTINATION bin)
