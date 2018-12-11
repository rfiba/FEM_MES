cmake_minimum_required(VERSION 3.10)
project(FEM)

set(CMAKE_CXX_STANDARD 17)

add_executable(FEM main.cpp Node.cpp Node.h Element.cpp Element.h Grid.cpp Grid.h functions.cpp functions.h)