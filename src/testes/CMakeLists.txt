cmake_minimum_required(VERSION 3.8)
project(testes)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
 add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(Eigen3 REQUIRED)
# Find OSQP and OsqpEigen packages
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)

include_directories(
 include
 ${EIGEN3_INCLUDE_DIRS}
 /usr/include/python3.12
)

add_library(trajectory 
 src/trajectory.cpp
)
add_library(pred_control 
 src/pred_control.cpp
)

add_executable(testosqp src/testosqp.cpp)

# Link against the libraries

target_link_libraries(pred_control
 OsqpEigen::OsqpEigen
 osqp::osqp
 ${EIGEN3_LIBRARIES}
)

target_link_libraries(testosqp
 pred_control
 trajectory
 ${Python3_LIBRARIES}
#  OsqpEigen::OsqpEigen
#  osqp::osqp
 ${EIGEN3_LIBRARIES}
)

ament_target_dependencies(testosqp
 rclcpp
#  osqp
#  OsqpEigen
 Eigen3)

target_compile_definitions(testosqp PRIVATE
 OSQP_EIGEN_OSQP_IS_V1
)

if(BUILD_TESTING)
 find_package(ament_lint_auto REQUIRED)
 set(ament_cmake_copyright_FOUND TRUE)
 set(ament_cmake_cpplint_FOUND TRUE)
 ament_lint_auto_find_test_dependencies()
endif()

install(
 TARGETS
 testosqp
 DESTINATION
 lib/${PROJECT_NAME}
)

ament_package()