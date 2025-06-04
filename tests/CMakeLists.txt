# 简单测试CMakeLists.txt
cmake_minimum_required(VERSION 3.5.0)

# 创建简单测试可执行文件
add_executable(
  simple_test
  simple_test.cpp
)

# 链接库
target_link_libraries(
  simple_test
  IMSimLib
  Eigen
)

# 包含头文件
target_include_directories(simple_test PRIVATE 
    ${CMAKE_SOURCE_DIR}/IMSim/include
)

# 添加自定义命令：构建完成后将tests目录下的txt文件复制到exe文件目录
add_custom_command(TARGET simple_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/EQSource.txt
        $<TARGET_FILE_DIR:simple_test>/EQSource.txt
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/SiteFile.txt
        $<TARGET_FILE_DIR:simple_test>/SiteFile.txt
    COMMENT "复制测试数据文件到exe目录"
)