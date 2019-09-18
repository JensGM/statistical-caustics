cmake_minimum_required(VERSION 3.14)
project(ocean-bin LANGUAGES C CXX)

add_executable(echo echo.cpp)
target_link_libraries(echo ocean::ocean)

install(
    TARGETS echo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
