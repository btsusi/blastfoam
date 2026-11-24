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

\*---------------------------------------------------------------------------*/

#include "SchnerrSauer.H"
#include "phasePair.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(SchnerrSauer, 0);
    addToRunTimeSelectionTable(massTransferModel, SchnerrSauer, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::SchnerrSauer::SchnerrSauer
(
    const dictionary& dict,
    const phasePair& pair
)
:
    massTransferModel(dict, pair),
    pSat_("pSat", dimPressure, dict),
    n_("n", dimless/dimVolume, dict),
    dNuc_("dNuc", dimLength, dict),
    Cc_("Cc", dimless, dict),
    Cv_("Cv", dimless, dict)
{
    if (!pair.ordered())
    {
        FatalErrorInFunction
            << "The SchnerrSauer mass transfer model must be specified for "
            << "an ordered phase pair, e.g. (vapor in liquid), where the "
            << "first phase is the vapour (dispersed) phase and the second "
            << "is the liquid (continuous) phase."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::SchnerrSauer::~SchnerrSauer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::massTransferModels::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc(n_*constant::mathematical::pi*pow3(dNuc_)/6);
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::SchnerrSauer::rRb
(
    const volScalarField& limitedAlphaL
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlphaL/(1.0 + alphaNuc() - limitedAlphaL),
        1.0/3.0
    );
}


Foam::tmp<Foam::volScalarField>
Foam::massTransferModels::SchnerrSauer::K() const
{
    const phaseModel& liquid = pair_.continuous();
    const phaseModel& vapor = pair_.dispersed();

    // Liquid volume fraction, limited to [0, 1]
    const volScalarField limitedAlphaL
    (
        "SchnerrSauer:limitedAlphaL",
        min(max(liquid, scalar(0)), scalar(1))
    );

    const volScalarField& rhoL = liquid.rho();
    const volScalarField& rhoV = vapor.rho();

    // Local mixture density of the pair
    const volScalarField rho
    (
        "SchnerrSauer:rho",
        limitedAlphaL*rhoL + (scalar(1) - limitedAlphaL)*rhoV
    );

    // Use the liquid (continuous) phase pressure to determine the
    // degree of cavitation
    const volScalarField& p = liquid.p();

    const dimensionedScalar p0(pSat_.dimensions(), 0);

    // Coefficient common to both the condensation and vaporisation rates
    const volScalarField pCoeff
    (
        "SchnerrSauer:pCoeff",
        (3*rhoL*rhoV)*sqrt(2.0/(3*rhoL))
       *rRb(limitedAlphaL)
       /(rho*sqrt(mag(p - pSat_) + 0.01*pSat_))
    );

    // Condensation: vapour -> liquid, active when p > pSat.
    // Proportional to the liquid already present.
    tmp<volScalarField> mDotCondensation
    (
        Cc_*limitedAlphaL*pCoeff*max(p - pSat_, p0)
    );

    // Vaporisation: liquid -> vapour, active when p < pSat.
    // Proportional to the liquid available to evaporate (regularised by
    // the nucleation-site volume fraction so that it remains bounded as
    // the liquid fraction tends to zero).
    tmp<volScalarField> mDotVaporisation
    (
        Cv_*(1.0 + alphaNuc() - limitedAlphaL)*pCoeff*max(pSat_ - p, p0)
    );

    // Net mass transfer rate into the vapour (dispersed) phase
    return volScalarField::New
    (
        "SchnerrSauer:mDot",
        mDotVaporisation - mDotCondensation
    );
}


// ************************************************************************* //
