add_library(electronicmap SHARED
        common.cpp
        config.cpp
        coordinates.cpp
        mapproj.cpp)

target_link_libraries(electronicmap
        ${THIRD_PARTY_LIBS})