/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2026
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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
    Foam::massTransferModels::SchnerrSauer

Description
    Schnerr-Sauer cavitation mass transfer model, adapted from the
    incompressible VoF formulation used by interPhaseChangeFoam to the
    compressible, pressure-non-equilibrium two-fluid Euler-Euler framework
    used by blastEulerFoam.

    The model must be specified for an ordered phase pair, e.g.

    \verbatim
        massTransfer
        (
            (vapor in liquid)
            {
                type    SchnerrSauer;

                //- Saturation vapour pressure
                pSat    2300;

                //- Bubble number density [1/m3]
                n       1.6e13;

                //- Nucleation site diameter [m]
                dNuc    2e-6;

                //- Condensation rate coefficient
                Cc      1;

                //- Vaporisation rate coefficient
                Cv      1;
            }
        );
    \endverbatim

    where "vapor" is the dispersed phase and "liquid" is the continuous
    phase of the pair. The model returns the net mass transfer rate into
    the vapour (dispersed) phase, K() = mDotVaporisation - mDotCondensation,
    with dimensions of density/time; this is consumed directly by
    phaseSystem::mDot/mDotU/mDotE/mDotByRho, which already add the
    corresponding sources to the continuity, momentum and energy equations
    of both phases (see phaseModel::solveAlphaRho and phaseModel::solve).

    Reference:
    \verbatim
        Schnerr, G. H., And Sauer, J.,
        "Physical and Numerical Modeling of Unsteady Cavitation Dynamics",
        Proc. 4th International Conference on Multiphase Flow,
        New Orleans, U.S.A., 2001.
    \endverbatim

SourceFiles
    SchnerrSauer.C

\*---------------------------------------------------------------------------*/

#ifndef SchnerrSauer_H
#define SchnerrSauer_H

#include "massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace massTransferModels
{

/*---------------------------------------------------------------------------*\
                        Class SchnerrSauer Declaration
\*---------------------------------------------------------------------------*/

class SchnerrSauer
:
    public massTransferModel
{
    // Private Data

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        //- Bubble number density
        dimensionedScalar n_;

        //- Nucleation site diameter
        dimensionedScalar dNuc_;

        //- Condensation rate coefficient
        dimensionedScalar Cc_;

        //- Vaporisation rate coefficient
        dimensionedScalar Cv_;


    // Private Member Functions

        //- Nucleation site volume-fraction
        dimensionedScalar alphaNuc() const;

        //- Reciprocal bubble radius as a function of the limited liquid
        //  (continuous phase) volume fraction
        tmp<volScalarField> rRb(const volScalarField& limitedAlphaL) const;


public:

    //- Runtime type information
    TypeName("SchnerrSauer");


    // Constructors

        //- Construct from a dictionary and a phase pair
        SchnerrSauer
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~SchnerrSauer();


    // Member Functions

        //- Return the net mass transfer rate into the dispersed
        //  (vapour) phase
        virtual tmp<volScalarField> K() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
