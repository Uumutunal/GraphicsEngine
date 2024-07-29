# SDL Graphics Engine
This project is a simple 3D graphics engine using SDL2, designed to render and manipulate 3D objects like meshes and perform basic 3D transformations. This project is for learning purposes only, not intended as a fully functional engine.

## Prerequisites

- SDL2 library

- ## Files

- `main.cpp`: The main application file.
- `structs.h`: Header file containing struct definitions.
- `vectorCalc.h`: Header file for vector calculations.
- `mesh.h`: Header file for mesh-related functionalities.
- `importer.h`: Header file for importing 3D models.

- ## Features

- **Z-Buffering**: Implementation of Z-buffering to handle depth calculations and hidden surface removal.
- **Line Drawing**: Functionality to draw lines with depth checking.
- **Triangle Drawing and Filling**: Functions to draw and fill triangles.
- **Projection**: Function to project 3D vertices to 2D screen space.
- **Transformations**: Functions to translate and rotate meshes.
- **Normal Calculation**: Calculation of surface normals for lighting and rendering.

- ## Issues

- Black artifact when apllying shading.
- Need to implement camera clipping.
