# varRhoTurbVOF


Under the VOF framework, the flow of the isothermal mixture belongs to the variable-density incompressible flow category. For such flows, VOF-based solvers of OpenFOAM fail to construct the correct governing equations for turbulence modeling. varRhoTurbVOF contains a set of newly designed VOF-based solvers which could use the desired governing equations for turbulence quantities.

Reference:

Wenyuan Fan and Henryk Anglart (2020). "varRhoTurbVOF: A new set of volume of fluid solvers for turbulent isothermal multiphase flows in OpenFOAM." Computer Physics Communications, 247, 106876.

## Installation:

Implementations for multiple versions of OpenFOAM are provided. In order to install the new solvers for a given version of OpenFOAM, one needs to enter the corresponding folder, load the environment variable, and execute:

./Allwmake

# Note

No overall update will be provided to this repository in the future since these versions are quite old. Implementations for new versions (starting from OpenFOAM-7 and OpenFOAM-1912) are provided in [varRhoTurbVOF_2](https://github.com/wenyuan-fan/varRhoTurbVOF_2), where the variable-density effect is selectable. In addition, [varRhoTurbVOF_2](https://github.com/wenyuan-fan/varRhoTurbVOF_2) provides the implementation for an advanced interfacial turbulence damping model. Please feel free to contact me should you have questions regarding these two repositories.
