## FDTD simulation of the Maxwell-Schrodinger system in MATLAB

If this MATLAB script is useful to you, please consider citing the following
reference in your work:

C. J. Ryu, A. Y. Liu, W. E. I. Sha, and W. C. Chew,
["Finite-difference time-domain simulation of the Maxwell-Schrodinger system
,"](https://ieeexplore.ieee.org/abstract/document/7558203) *IEEE J. Multiscale
Multiphys. Computat. Techn.*, vol. 1, pp. 40–47, 2016.

You can also refer to my [MS thesis](https://www.ideals.illinois.edu/items/89312).

This MATLAB script was written in 2014 by C. J. Ryu. It has the following novel
FDTD simulation features implemented based on the vector-potential formulation
of Maxwell's equations.

- Perfectly matched layers
- Plane-wave excitation
- Coupled Maxwell-Schrodinger system

The main script to run is [MAX_SCHROv4_3.m](MAX_SCHROv4_3.m).

## Compatibility issues

This script may not run on your verison of MATLAB since MathWorks has updated
MATLAB over the years. The following are fixes that you can apply.

- [insertHOGM_MS.m](insertHOGM_MS.m) uses `hermite` which is now called `hermiteH`
which requires the Symbolic Math Toolbox. A great alternative is [this MathWorks
File Exchange submission](https://www.mathworks.com/matlabcentral/fileexchange/27746-hermite-polynomials).
- [updateMaxv2.m](updateMaxv2.m) uses `heaviside` which also requires the Symbolic
Math Toolbox. You can simply implement this function using an if statement.

        function out = heaviside(n)
        
        if n < 0
            out = 0;
        elseif n == 0
            out = 0.5;
        else
            out = 1;
        end
        
        end

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).