cmake_minimum_required(VERSION 3.5.0)
project(SpatialIM)

# 设置静态链接，确保生成独立可执行文件
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -s") # 减小可执行文件体积

add_subdirectory(ThirdParty)
add_subdirectory(IMSim)
