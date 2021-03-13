add_library(coors SHARED
        common.cpp
        coorconfig.cpp
        coordinates.cpp
        mapproj.cpp)

target_link_libraries(coors
        ${THIRD_PARTY_LIBS})