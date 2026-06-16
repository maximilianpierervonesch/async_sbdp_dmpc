cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(distributed_mpc)

# Find LCM
find_package(PkgConfig REQUIRED)
pkg_check_modules(LCM REQUIRED lcm)


# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/messages)

# Source files
set(SOURCES
    src/Agent.cpp
)

# Create the agent executable
add_executable(agent src/Agent.cpp)
add_executable(simulator src/simulator.cpp)


# Link librariescd
target_link_libraries(agent PRIVATE lcm)
target_link_libraries(simulator PRIVATE lcm)
