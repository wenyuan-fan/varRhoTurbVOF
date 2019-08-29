# varRhoTurbVOF


Under the VOF framework, the flow of the isothermal mixture belongs to the variable-density incompressible flow category. For such flows, VOF-based solvers of OpenFOAM fail to construct the correct governing equations for turbulence modeling. varRhoTurbVOF contains a set of newly designed VOF-based solvers which could use the desired governing equations for turbulence quantities.

Reference:

Wenyuan Fan and Henryk Anglart. "varRhoTurbVOF: A new set of volume of fluid solvers for turbulent isothermal multiphase flows in OpenFOAM." Computer Physics Communications (2019), doi: 10.1016/j.cpc.2019.106876.

## Installation:

Implementations for multiple recent OpenFOAM versions are provided. In order to install the new solvers for a given version of OpenFOAM, one needs to enter the corresponding folder, load the environment variable, and execute:

./Allwmake

