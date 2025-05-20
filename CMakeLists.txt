cmake_minimum_required(VERSION 3.10)
project(stage_project)

set(CMAKE_CXX_STANDARD 17)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist Graf Gpad)

include_directories(${ROOT_INCLUDE_DIRS} inc)
link_directories(${ROOT_LIBRARY_DIR})

add_executable(MainApp
    src/main.cpp
)

target_link_libraries(MainApp ${ROOT_LIBRARIES})

