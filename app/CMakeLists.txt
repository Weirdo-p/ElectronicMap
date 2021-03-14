add_executable(coor coormain.cpp)
add_executable(dijks dijksmain.cpp)

target_link_libraries(coor coors ${THIRD_PARTY_LIBS})
target_link_libraries(dijks dijkstra ${THIRD_PARTY_LIBS})
