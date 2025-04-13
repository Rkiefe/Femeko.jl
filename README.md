# FEMjl
A wrapper for Gmsh, that simplifies the code required to create a Finite Element environment. Create your own models with cuboids and spheres, or import CAD files directly. Extract all the required information into simple matrices such as element connectivity, node coordinates, surface triangles of the 3D tetrahedral mesh of both inside objects and the outer container and more.
It also allows for a simplified local refinement of the mesh. It automatically selects the volumes inside and creates a local mesh refinement set by the user.

To install, make sure you add Gmsh and include the gmsh_wrapper.jl. That's it.
![twoBalls](https://github.com/user-attachments/assets/3b9549ba-3968-40f1-94a4-5c21ce37ca9e)
