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

# 添加自定义命令：构建完成后自动复制.pyd文件到python_package目录
add_custom_command(TARGET spatialim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:spatialim>
        ${CMAKE_SOURCE_DIR}/python_package/spatialim/
    COMMENT "Copying spatialim.pyd to python_package directory"
)

# 安装目标
install(TARGETS spatialim DESTINATION python)