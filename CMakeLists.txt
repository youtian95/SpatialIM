cmake_minimum_required(VERSION 3.5.0)
project(SpatialIM)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加pybind11的支持
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.11.1  # 使用最新的稳定版本
)
FetchContent_MakeAvailable(pybind11)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(ThirdParty)
add_subdirectory(IMSim)

# 添加Python绑定模块
add_subdirectory(python_bindings)
