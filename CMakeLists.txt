set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED true)

PROJECT(electronicmap)

set(CMAKE_BUILD_TYPE Debug)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

############# dependencies #################
# eigen
include_directories("/usr/include/eigen3")

AUX_SOURCE_DIRECTORY(src DIR_SRCS)
# add_executable(obsforecast ${DIR_SRCS})
# target_link_libraries(obsforecast pthread)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_subdirectory(src)
add_subdirectory(app)
