cmake_minimum_required(VERSION 3.15)
project(image_compresse_c)

set(CMAKE_CXX_STANDARD 14)

add_executable(image_compresse_c main.cpp Block.cpp Block.h Image.cpp Image.h)