cmake_minimum_required(VERSION 3.5.0)

# 创建Python模块
pybind11_add_module(spatialim spatialim_bindings.cpp)

# 链接到IMSim库
target_include_directories(spatialim PRIVATE 
    ${CMAKE_SOURCE_DIR}/IMSim/include
    ${pybind11_INCLUDE_DIRS}  # 添加pybind11的头文件路径
)

# 链接IMSimLib和Eigen库
target_link_libraries(spatialim PRIVATE IMSimLib Eigen)

# 设置输出目录
set_target_properties(spatialim PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)

# 安装目标
install(TARGETS spatialim DESTINATION python)