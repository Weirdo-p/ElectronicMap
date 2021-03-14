add_library(coors SHARED
        coorcommon.cpp
        coorconfig.cpp
        coordinates.cpp
        mapproj.cpp)

add_library(dijkstra SHARED
            dijksconfig.cpp
            dijkstra.cpp)

target_link_libraries(coors
        ${THIRD_PARTY_LIBS})

target_link_libraries(dijkstra
${THIRD_PARTY_LIBS})

