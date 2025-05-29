## Mixed Finite-Elements
Open-boundary problems such as solving the magnetostatic field of a magnetic body require the open space around the object to be included in the simulation volume. Additionally, the bounding shell must be quite a few times larger than the magnetic body for accurate simulations.

### Motivation
To skip the details, the Boundary-Element Method coupled with the Finite-Element Method (mixed finite elements) removes the necessity for a outer volume: only the magnetic body is considered.

This section was pretty much only developed to improve the micromagnetic simulations, as FEM/BEM is shown in literature to highly improve the simulation accuracy (and it did!). If this was simply a scaling problem, I would stay with the simpler FEM and just increase the number of elements.

### Drawbacks
The FEM implementation of the magnetostatic field in FEMjl scales ~ O(n) (n - number of volume elements), the BEM matrices scale ~ O(m^2) (m - number of surface elements). However these matrices only need to be calculated once.

BEM is, in my opinion, quite a bit harder than FEM, both in concept and in implementation.

### Conclusion
Here is a implementation of the mixed FEM/BEM method, to calculate the magnetostatic field of a magnetized object, which is then used for micromagnetic simulations. The implementation of BEM follows this article: https://doi.org/10.1016/j.jmmm.2012.01.016, but addapted to this use case, as the article shows a material with some magnetic permeability, while my implementation calculates the magnetostatic field of a magnetized object.
