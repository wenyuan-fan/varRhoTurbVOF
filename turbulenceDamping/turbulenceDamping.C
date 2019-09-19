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

#include "turbulenceDamping.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulenceDamping, 0);

    addToRunTimeSelectionTable
    (
        option,
        turbulenceDamping,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

volScalarField::Internal Foam::fv::turbulenceDamping::calculateSource
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
            // read properties for each phase
            const dictionary& transportProperties = mesh().lookupObject<IOdictionary>
            (
             "transportProperties"
            );

            const word phase1Name_(wordList(transportProperties.lookup("phases"))[0]);
            const word phase2Name_(wordList(transportProperties.lookup("phases"))[1]);

            const dictionary& phase1_ = transportProperties.subDict(phase1Name_);
            const dictionary& phase2_ = transportProperties.subDict(phase2Name_);

            const dimensionedScalar& rho1_ = phase1_.lookup("rho");
            const dimensionedScalar& rho2_ = phase2_.lookup("rho");

            const dimensionedScalar& nu1_ = phase1_.lookup("nu");
            const dimensionedScalar& nu2_ = phase2_.lookup("nu");


            const volScalarField& Alpha = mesh().lookupObject<volScalarField>("alpha." + phase1Name_);
            const volVectorField grad_Alpha = fvc::grad(Alpha);
            const volScalarField grad_Alpha_mag = mag(grad_Alpha);

            scalar beta = 0.075;


            // calculate interfacial area density
            volScalarField::Internal A1 = 2.0*Alpha*grad_Alpha_mag;
            volScalarField::Internal A2 = 2.0*(1.0-Alpha)*grad_Alpha_mag;


            // calculate the inverse of the length scale
            const volScalarField::Internal& V = mesh_.V();

            volScalarField oneByDn
            (
            IOobject
            (
            "oneByDn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("oneByDn",dimensionSet(0,-1,0,0,0,0,0),0.0)
            );

            if (lengthScale_ == "FA")
            {
                forAll(oneByDn, celli)
                {
                    if (grad_Alpha_mag[celli] > SMALL)
                    {
                        scalar projectedArea = 0.0;
                        const cell& c = mesh_.cells()[celli];
                        forAll(c, facej)
                        {
                            projectedArea += mag(grad_Alpha[celli] & mesh_.Sf()[c[facej]]);
                        }

                        oneByDn[celli] = 0.5*projectedArea/V[celli]/grad_Alpha_mag[celli];
                    }
                }
            }

            else if (lengthScale_ == "cubicRoot")
            {
                oneByDn.ref() = pow(V,-1.0/3.0);
            }


            // calculate separate damping terms
            volScalarField::Internal coeffs = 36.0*sqr(B_)/beta*pow(oneByDn, 3.0);
            volScalarField::Internal source1 = coeffs*A1*rho1_*sqr(nu1_);
            volScalarField::Internal source2 = coeffs*A2*rho2_*sqr(nu2_);

            // calculate the total damping term
            dimensionedScalar heavy("heavy", dimless, 0.0);

            volScalarField::Internal source = 0.0 * source1;

            if (dampingTreatment_ == "heavyNegative")
            {
                if (rho1_ > rho2_)
                {
                    heavy = - rho2_/rho1_*sqr(nu2_)/sqr(nu1_);

                    source = source1*heavy + source2;
                }

                else
                {
                    heavy = - rho1_/rho2_*sqr(nu1_)/sqr(nu2_);

                    source = source1 + source2*heavy;
                }   
            }

            else if (dampingTreatment_ == "heavyZero")
            {
                if (rho1_ > rho2_)
                {
                    source = source2;
                }

                else
                {
                    source = source1;
                }   
            }

            else if (dampingTreatment_ == "symmetric")
            {
                source = source1 + source2;   
            }

            return source;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulenceDamping::turbulenceDamping
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    B_("B", dimless, coeffs_.lookup("B")),
    lengthScale_(coeffs_.lookupOrDefault<word>("lengthScale", "FA")),
    dampingTreatment_(coeffs_.lookupOrDefault<word>("dampingTreatment", "heavyNagative"))
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

void Foam::fv::turbulenceDamping::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& Rho = mesh().lookupObject<volScalarField>("rho");
    const volScalarField& omega = mesh().lookupObject<volScalarField>(fieldNames_[fieldi]);
    volScalarField::Internal source = calculateSource(eqn, fieldi)/omega/Rho;
    eqn += fvm::Sp(source, omega);
}


void Foam::fv::turbulenceDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const volScalarField& omega = mesh().lookupObject<volScalarField>(fieldNames_[fieldi]);
    volScalarField::Internal source = calculateSource(eqn, fieldi)/omega;
    eqn += fvm::Sp(source, omega);
}

bool Foam::fv::turbulenceDamping::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}

// ************************************************************************* //
