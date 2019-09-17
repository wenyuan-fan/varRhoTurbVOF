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

Class
    Foam::selectableIncompressibleTransportModel

Description
    Transport model selection class for the compressibleInterFoam family of
    solvers.

    By default the standard mixture transport modelling approach is used in
    which a single momentum stress model (laminar, non-Newtonian, LES or RAS) is
    constructed for the mixture.  However if the \c simulationType in
    constant/turbulenceProperties is set to \c twoPhaseTransport the alternative
    Euler-Euler two-phase transport modelling approach is used in which separate
    stress models (laminar, non-Newtonian, LES or RAS) are instantiated for each
    of the two phases allowing for different modeling for the phases.

SourceFiles
    selectableIncompressibleTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef selectableIncompressibleTransportModel_H
#define selectableIncompressibleTransportModel_H

#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "varRhoTurbulentTransportModel.H"
#include "selectableIncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class selectableIncompressibleTransportModel Declaration
\*---------------------------------------------------------------------------*/
class selectableIncompressibleTransportModel
{
    // Private Data

        //- Switch to select two-phase or mixture transport modelling
        Switch varRhoIncompressible_;

        //- Two-phase mixture
        const immiscibleIncompressibleTwoPhaseMixture& mixture_;


        //- Mixture transport model (constructed for mixture transport)
        autoPtr<varRhoIncompressible::turbulenceModel> turbulenceVarRho_;

        autoPtr<incompressible::turbulenceModel> turbulenceStrict_;


public:

    // Constructors

        //- Construct from components
        selectableIncompressibleTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const surfaceScalarField& rhoPhi,
            const immiscibleIncompressibleTwoPhaseMixture& mixture
        );

        //- Disallow default bitwise copy construction
        selectableIncompressibleTransportModel
        (
            const selectableIncompressibleTransportModel&
        );


    // Member Functions


        //- Return the effective momentum stress divergence
        tmp<fvVectorMatrix> divDevRhoReff(volScalarField& rho, volVectorField& U) const;

        //- Correct the phase or mixture transport models
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const selectableIncompressibleTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
