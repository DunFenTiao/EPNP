cmake_minimum_required( VERSION 2.8)

project(EPNP)
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )
set( OpenCV_DIR /Users/apple/Documents/opencv-3.4.3/release)
find_package( OpenCV REQUIRED) 


include_directories( 
    ${OpenCV_INCLUDE_DIRS} 
    "/usr/local/include/eigen3/"
)

add_executable(opencvtest test.cpp)
target_link_libraries(opencvtest
	${OpenCV_LIBS}
) 