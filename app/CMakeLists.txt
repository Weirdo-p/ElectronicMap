add_executable(coortest coormain.cpp)
add_executable(dijkstest dijksmain.cpp)

target_link_libraries(coortest coors ${THIRD_PARTY_LIBS})
target_link_libraries(dijkstest ${THIRD_PARTY_LIBS})
