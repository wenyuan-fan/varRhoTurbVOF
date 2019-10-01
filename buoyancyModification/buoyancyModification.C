/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "buoyancyModification.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(buoyancyModification, 0);

    addToRunTimeSelectionTable
    (
        option,
        buoyancyModification,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::buoyancyModification::buoyancyModification
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    sigmaB_(coeffs_.lookupOrDefault<scalar>("sigmaB", 0.85)),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    coeffs_.lookup("fields") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::buoyancyModification::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& nut_ =
        mesh().lookupObject<volScalarField>("nut");

    const volScalarField& k_ =
        mesh().lookupObject<volScalarField>(fieldNames_[fieldi]);

    volScalarField::Internal Gb = -nut_/sigmaB_*(fvc::grad(rho)&g_);

    eqn += fvm::Sp(Gb/k_, k_);
}

bool Foam::fv::buoyancyModification::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}

// ************************************************************************* //
