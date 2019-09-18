cmake_minimum_required(VERSION 3.14)
project(ocean-lib LANGUAGES C CXX)

add_library(ocean src/ocean.cpp)
add_library(ocean::ocean ALIAS ocean)
target_include_directories(ocean
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(ocean
    PUBLIC
        FFTW3::fftw3
        glm
)

# install and export
install(
    TARGETS ocean
    EXPORT ocean-export
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
    DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(
    EXPORT ocean-export
    NAMESPACE ocean::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/ocean/cmake
    FILE ocean-config.cmake
)
export(
    TARGETS ocean
    NAMESPACE ocean::
    FILE ocean-config.cmake
)
