# 3D-Rendering-Library

This is 3D Rendering library built in C++ for the course CSCI 580.

The library can handle rasterization, transformations, lighting and shading, texture maps and anti-aliasing.

The project solution is for VS 2015.

The starter code was provided by Prof. Neumann for the course, but the entire library was implemented by students over a course of 2 months.

# Important Files

1. gz.h: Defines the various constants and datatypes. Like OpenGL has gl, this library starts everything with gz

2. rend.h: All the methods and different objects used in the library can be found here

3. rend.cpp: Contains the implementation of the main library. Comments provided as needed.

4. tex_fun.cpp: Provides functionality for image textures and procedural textures (Julia set).

5. Application5.cpp: The main application function that initialises and calls the renderer (with anti-aliasing). The types of camera (default/custom) and types of textures can be changed here.