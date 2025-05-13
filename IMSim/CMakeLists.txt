file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)

# 设置静态链接
if(MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MT")
else()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

# 创建一个库，供其他目标使用
add_library(IMSimLib STATIC ${srcs})
target_include_directories(IMSimLib PUBLIC include)
target_link_libraries(IMSimLib Eigen)

# 创建可执行文件
add_executable(IMSim ${srcs})
target_include_directories(IMSim PUBLIC include)
target_link_libraries(IMSim Eigen)

