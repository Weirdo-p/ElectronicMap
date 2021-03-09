add_library(electronicmap SHARED
        common.cpp
        config.cpp
        coordinates.cpp)

target_link_libraries(electronicmap
        ${THIRD_PARTY_LIBS})