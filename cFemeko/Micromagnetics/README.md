## C++ version of Femeko.jl Landau-Lifshitz micromagnetics solver

### How to use
- Edit the LL.h to set the solver parameters such as magnetic saturation and number of iteration steps, etc.
- Compile with `g++ -O3 -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp`
- run the `permalloy.jl` to solve for the input parameters in .h and the input magnetization in .jl