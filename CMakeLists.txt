cmake_minimum_required(VERSION 3.9)

project(VoidExample LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################
## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND VoidExample_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND VoidExample_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND VoidExample_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND VoidExample_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND VoidExample_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND VoidExample_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################

list(APPEND VoidExample_SOURCES ${VoidExample_sources})
list(APPEND VoidExample_HEADERS ${VoidExample_headers})
list(APPEND VoidExample_INCLUDE ${VoidExample_includes})

# Create executable
################################################################################
add_executable(${PROJECT_NAME}
    main.cpp
    ${VoidExample_SOURCES}
    ${VoidExample_HEADERS})
    

target_link_libraries(${PROJECT_NAME} ${VoidExample_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${VoidExample_INCLUDE})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)
