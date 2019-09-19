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
    Foam::fv::turbulenceDamping

Description
    Calculates and applies the turbulence damping term to the omega equation.
    The current implementation is based on the following work:
    \verbatim
        Fan, W. & Anglart, H. (2019).
        Progress in Phenomenological Modeling of Turbulence Damping around a
        Two-Phase Interface.
        Fluids, 4(3), 136.
    \endverbatim
    where modifications are made to the model proposed by
    \verbatim
        Egorov, Y. (2004).
        Validation of CFD codes with PTS-relevant test cases.
        5th Euratom Framework Programme ECORA project, pp. 91–116.
    \endverbatim


Usage
    Example usage:
    \verbatim
    fields             (omega);   // Name of omega field

    B                  0;         // Damping coefficient, 0 means no damping

    lengthScale        FA;        // Optional parameter to specify the method
                                  // to calculate the length scale for damping.
                                  // The alternative is "cubicRoot".

    dampingTreatment   heavyNegative; // Optional parameter to specify the 
                                      // treatment for the heavier phase.
                                      // Alternatives are "heavyZero" and
                                      // "symmetric".
    \endverbatim

Author
    Wenyuan Fan

SourceFiles
    turbulenceDamping.C

\*---------------------------------------------------------------------------*/

#ifndef turbulenceDamping_H
#define turbulenceDamping_H
#include "fvCFD.H"
#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class turbulenceDamping Declaration
\*---------------------------------------------------------------------------*/

class turbulenceDamping
:
    public option
{
    // Private Data

        dimensionedScalar B_;

        word lengthScale_;

        word dampingTreatment_;

        //- Source term to omega equation
        volScalarField::Internal calculateSource
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

public:

    //- Runtime type information
    TypeName("turbulenceDamping");


    // Constructors

        //- Construct from explicit source name and mesh
        turbulenceDamping
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        turbulenceDamping(const turbulenceDamping&) = delete;


    // Member Functions

        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
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
        void operator=(const turbulenceDamping&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
