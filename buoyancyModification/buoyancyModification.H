/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::fv::buoyancyModification

Description
    Introduces a buoyancy correction term to the k equation.

    Reference:
    \verbatim
        Devolder, B., Rauwoens, P., & Troch, P. (2017).
        Application of a buoyancy-modified k-ω SST turbulence model to
        simulate wave run-up around a monopile subjected to regular waves
        using OpenFOAM.
        Coastal Engineering, 125, 81–94.
    \endverbatim


Usage
    Example usage:
    \verbatim
    fields             (k);   // Name of k field

    sigmaB             0.85;  // Optional coefficient, set to 0.85 by default
    \endverbatim

Author
    Wenyuan Fan

SourceFiles
    buoyancyModification.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyModification_H
#define buoyancyModification_H
#include "fvCFD.H"
#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class buoyancyModification Declaration
\*---------------------------------------------------------------------------*/

class buoyancyModification
:
    public option
{
    // Private Data

        scalar sigmaB_;

        uniformDimensionedVectorField g_;

public:

    //- Runtime type information
    TypeName("buoyancyModification");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyModification
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        buoyancyModification(const buoyancyModification&) = delete;


    // Member Functions

        // Evaluate

            //- Add turbulence damping to full-form k equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyancyModification&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
