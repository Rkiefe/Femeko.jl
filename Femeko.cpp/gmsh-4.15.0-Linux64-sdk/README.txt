This is the binary Software Development Kit (SDK) for Gmsh 4.15.0:

  * Operating system: Linux64-sdk (Linux)
  * C++ compiler: /usr/local/bin/g++
  * C++ compiler ID: GNU
  * C++ compiler version: 6.5.0
  * C++ compiler flags:-fPIC -fopenmp -O2 -g -DNDEBUG
  * Build options: 64Bit ALGLIB[contrib] ANN[contrib] Bamg Blas[petsc] Blossom Cgns DIntegration Dlopen DomHex Eigen[contrib] Fltk Gmm[contrib] Hxt Jpeg Kbipack Lapack[petsc] LinuxJoystick MathEx[contrib] Med Mesh Metis[contrib] Mmg Mpeg Netgen Nii2mesh ONELAB ONELABMetamodel OpenCASCADE OpenCASCADE-CAF OpenGL OpenMP OptHom PETSc Parser Plugins Png Post QuadMeshingTools QuadTri Solver TetGen/BR TinyXML2[contrib] Untangle Voro++[contrib] WinslowUntangler Zlib tinyobjloader

Gmsh is distributed under the terms of the GNU General Public License: see
share/doc/gmsh/LICENSE.txt and share/doc/gmsh/CREDITS.txt. For additional Gmsh
resources, see https://gmsh.info.

SDK layout:

  * lib/*gmsh*.{so,dll,dylib}: shared Gmsh library
  * lib/gmsh.dll.lib: import library (Windows only)
  * lib/gmsh.py: Python module
  * lib/gmsh.jl: Julia module
  * include/gmsh.h: C++ API header
  * include/gmshc.h: C API header
  * include/gmsh.h_cwrap: C++ wrapper of the C API (see the `Notes' below)
  * include/gmsh.f90: Fortran module
  * bin/gmsh: gmsh executable (linked with the shared Gmsh library)
  * share/doc/gmsh/tutorials/c++: C++ API tutorials
  * share/doc/gmsh/tutorials/c: C API tutorials
  * share/doc/gmsh/tutorials/python: Python API tutorials
  * share/doc/gmsh/tutorials/julia: Julia API tutorials
  * share/doc/gmsh/tutorials/fortran: Fortran API tutorials
  * share/doc/gmsh/examples/api: Other API examples

Notes:

  * The C API should work with most compilers.

  * The C++ API will only work if your compiler has an Application Binary
    Interface (ABI) that is compatible with the ABI of the compiler used to
    build this SDK (see above for the compiler ID and version).

  * If your C++ compiler does not have a compatible ABI and if there are no
    compatibility flags available, you can rename `gmsh.h_cwrap' as `gmsh.h':
    this implementation redefines the C++ API in terms of the C API. Using this
    header will lead to (slightly) reduced performance compared to using the
    native Gmsh C++ API from the original `gmsh.h' header, as it entails
    additional data copies between this C++ wrapper, the C API and the native
    C++ code.

    For example, the Windows SDK is currently compiled using the GNU Compiler
    Collection (GCC). To compile a C++ example with Microsoft Visual Studio 2017
    in the Visual Studio shell and run it, you would do:

    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h gmsh.h_original
    C:\gmsh-git-Windows64-sdk> ren include\gmsh.h_cwrap gmsh.h
    C:\gmsh-git-Windows64-sdk> cl /Iinclude share\doc\gmsh\tutorials\c++\t1.cpp lib\gmsh.dll.lib
    C:\gmsh-git-Windows64-sdk> cd lib
    C:\gmsh-git-Windows64-sdk\lib> ..\t1.exe

  * To make it as portable and as easy-to-use as possible the shared Gmsh
    library embeds most Gmsh dependencies statically linked directly inside the
    share library: libgfortran, FreeType, OpenCASCADE, OpenBLAS, HDF5, MED,
    CGNS, FLTK, ... Linking your app with different versions of some of these
    dependencies can cause issues: in this case you should rebuild the Gmsh
    shared library from source and manage the dependencies consistently.

  * The shared Gmsh library also references shared system libraries. If some of
    these libraries are missing you will need to install them.
