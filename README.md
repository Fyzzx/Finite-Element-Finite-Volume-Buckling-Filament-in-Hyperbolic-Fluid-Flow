This code creates an initially straight filament immersed in a fluid of specified viscosity and flow (inputs make it easy to create a hyperbolic shaped flow)
The filament then moves and buckles under the fluidic drag and thermodynamic forces. 

Code is written in MATLAB, but can easily be transferred to Python and can easily be extended to parallelized methods because this code runs quickly, therefore new runs can be performed on alternate cores/threads. 

This code is from our paper:
P. M. Ryan and C. W. Wolgemuth. A finite volume algorithm for the dynamics of filaments,
rods, and beams. Journal of Computational Physics, 466:111375, 2022a. ISSN 0021-9991.
doi: https://doi.org/10.1016/j.jcp.2022.111375.


And it matches data from the following papers:
Y.-N. Young, M.J. Shelley, Stretch-coil transition and transport of fibers in cellular flows, Phys. Rev. Lett. 99 (2007) 058303, https://doi.org/10.1103/PhysRevLett.99.058303.
V. Kantsler, R.E. Goldstein, Fluctuations, dynamics, and the stretch-coil transition of single actin filaments in extensional flows, Phys. Rev. Lett. 108 (2012) 038103, https://doi.org/10.1103/PhysRevLett.108.038103.
R. Chelakkot, R.G. Winkler, G. Gompper, Flow-induced helical coiling of semiflexible polymers in structured microchannels, Phys. Rev. Lett. 109 (2012) 178101, https://doi.org/10.1103/PhysRevLett.109.178101.
B. Chakrabarti, Y. Liu, J. LaGrone, R. Cortez, L. Fauci, O. du Roure, D. Saintillan, A. Lindner, Flexible filaments buckle into helicoidal shapes in strong compressional flows, Nat. Phys. 16 (2020) 689, https://doi.org/10.1038/s41567-020-0843-7.
H. Manikantan, D. Saintillan, Buckling transition of a semiflexible filament in extensional flow, Phys. Rev. E 92 (2015) 041002, https://doi.org/10.1103/PhysRevE.92.041002.


![image](https://github.com/Fyzzx/Finite-Element-Finite-Volume-Buckling-Filament-in-Hyperbolic-Fluid-Flow/assets/103218124/e9297c18-0629-485f-a81c-a85a0156d192)
