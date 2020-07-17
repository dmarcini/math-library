cmake_minimum_required(VERSION 3.10)

project(math-library)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

set(SOURCES "src/null.cpp")

add_library(math STATIC ${SOURCES})

set(HEADERS "include/matrix.hpp;"
	    "include/math_operations.hpp")

set_target_properties(math PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_link_libraries(math utility)

install(TARGETS math
        LIBRARY DESTINATION /usr/local/lib
        PUBLIC_HEADER DESTINATION /usr/local/include)

