cmake_minimum_required(VERSION 3.16)

project(linalg)


add_executable(

    example.out

    example.cpp

    source/cpp/linalg.cpp
    source/cpp/Vector.cpp
    source/cpp/Polynomial.cpp
    source/cpp/Matrix.cpp

)