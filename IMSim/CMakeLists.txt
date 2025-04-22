file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)

# 设置静态链接
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

add_executable(IMSim ${srcs})
target_include_directories(IMSim PUBLIC include)

target_link_libraries(IMSim Eigen)

