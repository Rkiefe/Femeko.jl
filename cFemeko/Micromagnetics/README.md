# C++ version of Femeko.jl Landau-Lifshitz micromagnetics solver
Expect about 2x the performance with C++.

## How to use
- Edit the LL.h to set the solver parameters such as magnetic saturation and number of iteration steps, etc.
- Compile with `g++ -O3 -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp`
- run the `permalloy.jl` to solve for the input parameters in .h and the input magnetization in .jl

### Demagnetizing field Julia vs C++
<img width="1920" height="1032" alt="demag_Julia_vs_C++" src="https://github.com/user-attachments/assets/9b62f606-3f05-4b40-8a09-0a160b687864" />

### Expected result with default parameters
<img width="602" height="482" alt="cFemeko_micromagnetics" src="https://github.com/user-attachments/assets/dd97c2d8-cbf8-48a6-8931-7ca95b922ca5" />
