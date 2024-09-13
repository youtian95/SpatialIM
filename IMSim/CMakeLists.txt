file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)

add_executable(IMSim ${srcs})
target_include_directories(IMSim PUBLIC include)

target_link_libraries(IMSim Eigen)

