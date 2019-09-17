/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "selectableIncompressibleTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::selectableIncompressibleTransportModel::selectableIncompressibleTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    varRhoIncompressible_(true),
    mixture_(mixture)
{
    {
        IOdictionary turbulenceProperties
        (
            IOobject
            (
                turbulenceModel::propertiesName,
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word incompressibleType
        (
            turbulenceProperties.lookup("incompressibleType")
        );

        if (incompressibleType == "strictIncompressible")
        {
            varRhoIncompressible_ = false;
            Info << "Strict incompressible turbulence modeling"
                 << endl;
        }

        else if (incompressibleType == "varRhoIncompressible")
        {
            Info << "Variable-density incompressible turbulence modeling"
                 << endl;
        }

        else
        {
                FatalError
                    << "Valid types for incompressibleType are" << nl
                    << "strictIncompressible and varRhoIncompressible" << nl
                    << exit(FatalError);
        }
    }

    if (varRhoIncompressible_)
    {
        turbulenceVarRho_ = varRhoIncompressible::turbulenceModel::New
        (
            rho,
            U,
            rhoPhi,
            phi,
            mixture
        );

        turbulenceVarRho_->validate();
    }
    else
    {
        turbulenceStrict_ = incompressible::turbulenceModel::New
        (
            U,
            phi,
            mixture
        );

        turbulenceStrict_->validate();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::fvVectorMatrix>
Foam::selectableIncompressibleTransportModel::divDevRhoReff
(
    volScalarField& rho,
    volVectorField& U
) const
{
    if (varRhoIncompressible_)
    {
        return turbulenceVarRho_->divDevRhoReff(U);
    }
    else
    {
        return turbulenceStrict_->divDevRhoReff(rho, U);
    }
}


void Foam::selectableIncompressibleTransportModel::correct()
{
    if (varRhoIncompressible_)
    {
        turbulenceVarRho_->correct();
    }
    else
    {
        turbulenceStrict_->correct();
    }
}


// ************************************************************************* //
