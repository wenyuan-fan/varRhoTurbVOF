/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "VarRhoIncompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"


#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeTurbulenceModelTypes
(
    geometricOneField,
    volScalarField,
    varRhoIncompressibleTurbulenceModel,
    VarRhoIncompressibleTurbulenceModel,
    transportModel
);

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    varRhoIncompressibleTurbulenceModel,
    VarRhoIncompressibleTurbulenceModel,
    transportModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedTurbulenceModel                                               \
    (transportModelVarRhoIncompressibleTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelVarRhoIncompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelVarRhoIncompressibleTurbulenceModel, LES, Type)


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "generalizedNewtonian.H"
makeLaminarModel(generalizedNewtonian);

#include "Maxwell.H"
makeLaminarModel(Maxwell);

#include "Giesekus.H"
makeLaminarModel(Giesekus);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "SpalartAllmaras.H"
makeRASModel(SpalartAllmaras);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "RNGkEpsilon.H"
makeRASModel(RNGkEpsilon);

#include "realizableKE.H"
makeRASModel(realizableKE);

#include "LaunderSharmaKE.H"
makeRASModel(LaunderSharmaKE);

#include "kOmega.H"
makeRASModel(kOmega);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSAS.H"
makeRASModel(kOmegaSSTSAS);

#include "kOmegaSSTLM.H"
makeRASModel(kOmegaSSTLM);

#include "v2f.H"
makeRASModel(v2f);

#include "LRR.H"
makeRASModel(LRR);

#include "SSG.H"
makeRASModel(SSG);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "WALE.H"
makeLESModel(WALE);

#include "kEqn.H"
makeLESModel(kEqn);

#include "dynamicKEqn.H"
makeLESModel(dynamicKEqn);

#include "dynamicLagrangian.H"
makeLESModel(dynamicLagrangian);

#include "SpalartAllmarasDES.H"
makeLESModel(SpalartAllmarasDES);

#include "SpalartAllmarasDDES.H"
makeLESModel(SpalartAllmarasDDES);

#include "SpalartAllmarasIDDES.H"
makeLESModel(SpalartAllmarasIDDES);

#include "DeardorffDiffStress.H"
makeLESModel(DeardorffDiffStress);

#include "kOmegaSSTDES.H"
makeLESModel(kOmegaSSTDES);


// ************************************************************************* //
