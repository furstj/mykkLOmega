/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "mykkLOmega.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mykkLOmega, 0);
addToRunTimeSelectionTable(RASModel, mykkLOmega, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> mykkLOmega::fv(const volScalarField& Ret) const
{
    return(1.0 - exp(-sqrt(Ret)/Av_));
}


tmp<volScalarField> mykkLOmega::fINT() const
{
    return
    (
        min
        (
	 kt_/(Cint_*max(kl_ + kt_, kMin_)),
            dimensionedScalar("1.0", dimless, 1.0)
        )
    );
}


tmp<volScalarField> mykkLOmega::fSS(const volScalarField& Omega) const
{
  return(exp(-sqr(Css_*nu()*Omega/max(kt_,kMin_))));
}


tmp<volScalarField> mykkLOmega::Cmu(const volScalarField& S) const
{
    return(1.0/(A0_ + As_*(S/omega_)));
}


  tmp<volScalarField> mykkLOmega::BetaTS(const volScalarField& Rew, 
					 const volScalarField& L) const
{
  return(scalar(1) - exp(-sqr(max(Rew - CTSCrit(L), scalar(0)))/Ats_));
}


tmp<volScalarField> mykkLOmega::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -CtauL_*ktL
          /
            (
                sqr
                (
		 max(lambdaEff*Omega,
		     dimensionedScalar
		     (
		      "ROOTVSMALL",
		      dimLength*inv(dimTime),
		      ROOTVSMALL
		      )
		     )
		 )
	     )
	 )
     );
}
  

tmp<volScalarField> mykkLOmega::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}


tmp<volScalarField> mykkLOmega::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
            -0.41
           * pow4
            (
                lambdaEff
		/ max(lambdaT, 
		      dimensionedScalar
		      (
		       "ROTVSMALL",
		       lambdaT.dimensions(),
		       ROOTVSMALL
		       )
		      )
	     )
	 )
     );
}
  

tmp<volScalarField> mykkLOmega::gammaBP(const volScalarField& Omega) const
{
    if (medinaPhiBP_)
        return
            (
                max
                (
                    sqrt(kt_)*y_/nu()
                    - CbpCrit_,
                    scalar(0)
                )
            );

    else
        return
            (
                max
                (
                    kt_/nu()
                    / max(
                        Omega,
                        dimensionedScalar("ROOTVSMALL", Omega.dimensions(), ROOTVSMALL)
                    )
                    - CbpCrit_,
                    scalar(0)
                )
            );
}


tmp<volScalarField> mykkLOmega::gammaNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit,
    const volScalarField& L
) const
{
    return
        (
            max
            (
                ReOmega
                - CnatCrit(L)
                / max(
                    fNatCrit, dimensionedScalar("ROTVSMALL", dimless, ROOTVSMALL)
                ),
                scalar(0)
            )
        );
}

tmp<volScalarField> mykkLOmega::L(const volScalarField& Rew) const
{
  const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
  
  dimensionedScalar pTot("pTot", p.dimensions(),
			 gMax( volScalarField( p + 0.5*magSqr(U_) ) ) );

  dimensionedScalar uMin("uMin", dimVelocity, VSMALL);

  volScalarField Ue = sqrt( 2.0 * (pTot - p) ); 

  volScalarField dpdx = (fvc::grad(p) & U_) / max( mag(U_), uMin); 
  
  volScalarField K = - nu() / pow3(max(Ue,uMin)) * dpdx;

  return sqr(Rew) * K;
}

  
tmp<volScalarField> mykkLOmega::CTSCrit(const volScalarField& L) const
{

  if (furstPG_) {
    return 536.40 / (1.0 - CtsApg_ * min(L, 0.0));
  } 
  else 
    return CtsCrit_ * min( max(L, 1.0), 1.0); 
}

tmp<volScalarField> mykkLOmega::CnatCrit(const volScalarField& L) const
{

  if (furstPG_) {
    return CnatCrit_ / (1.0 - CnatApg_ * min(L, 0.0));
  } 
  else 
    return CnatCrit_ * min( max(L, 1.0), 1.0); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mykkLOmega::mykkLOmega
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
    :
    RASModel(modelName, U, phi, transport, turbulenceModelName),
    
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.04
        )
    ),
    As_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "As",
            coeffDict_,
            2.12
        )
    ),
    Av_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Av",
            coeffDict_,
            6.75
        )
    ),
    Abp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Abp",
            coeffDict_,
            0.6
        )
    ),
    Anat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anat",
            coeffDict_,
            200
        )
    ),
    Ats_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ats",
            coeffDict_,
            200
        )
    ),
    CbpCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CbpCrit",
            coeffDict_,
            1.2
        )
    ),
    Cnc_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cnc",
            coeffDict_,
            0.1
        )
    ),
    CnatCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatCrit",
            coeffDict_,
            1250
        )
    ),
    CnatApg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CnatApg",
            coeffDict_,
            8.963
        )
    ),
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            0.75
        )
    ),
    CtsCrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsCrit",
            coeffDict_,
            1000
        )
    ),
    CtsApg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtsApg",
            coeffDict_,
            8.963
        )
    ),
    CrNat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CrNat",
            coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            coeffDict_,
            0.12
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            coeffDict_,
            0.035
        )
    ),
    Css_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Css",
            coeffDict_,
            1.5
        )
    ),
    CtauL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CtauL",
            coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.92
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            0.3
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            coeffDict_,
            1.5
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            coeffDict_,
            2.495
        )
    ),
    CmuStd_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuStd",
            coeffDict_,
            0.09
        )
    ),
    Prtheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prtheta",
            coeffDict_,
            0.85
        )
    ),
    Sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmak",
            coeffDict_,
            1
        )
    ),
    Sigmaw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Sigmaw",
            coeffDict_,
            1.17
        )
    ),

    kMax_("kMax", sqr(dimVelocity), HUGE),

    omegaMax_("omegaMax", inv(dimTime), HUGE),

    kt_
    (
        IOobject
        (
            "kt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("kt", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    kl_
    (
        IOobject
        (
            "kl",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("kl", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),

    y_(mesh_),

    medinaPhiBP_
    (
        Switch::lookupOrAddToDict
        (
            "medinaPhiBP",
            coeffDict_,
            false
        )
     ),

    furstPG_
    (
        Switch::lookupOrAddToDict
        (
            "furstPG",
            coeffDict_,
            false
        )
     ),

    CDT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDT",
            coeffDict_,
            1.0
        )
     ),

    viscosityRatioLimit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "viscosityRatioLimit",
            coeffDict_,
            1.0e5
        )
     ),

    fv_
    (
        IOobject
        (
            "fv",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
     ),

    fINT_
    (
        IOobject
        (
            "fINT",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
     ),

    fSS_
    (
        IOobject
        (
            "fSS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    Cmu_
    (
        IOobject
        (
            "Cmu",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    BetaTS_
    (
        IOobject
        (
            "BetaTS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    BetaNAT_
    (
        IOobject
        (
            "BetaNAT",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    BetaBP_
    (
        IOobject
        (
            "BetaBP",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    fTaul_
    (
        IOobject
        (
            "fTaul",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    fOmega_
    (
        IOobject
        (
            "fOmega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    gammaBP_
    (
        IOobject
        (
            "gammaBP",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    gammaNAT_
    (
        IOobject
        (
            "gammaNAT",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    ReOmega_
    (
        IOobject
        (
            "ReOmega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0
    ),

    RBP_
    (
        IOobject
        (
            "RBP",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("zero", kt_.dimensions()/dimTime, 0)
    ),

    RNAT_
    (
        IOobject
        (
            "RNAT",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("zero", kt_.dimensions()/dimTime, 0)
     ),

    Pkt_
    (
        IOobject
        (
            "Pkt",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("zero", kt_.dimensions()/dimTime, 0)
     ),

    Pkl_
    (
        IOobject
        (
            "Pkl",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("zero", kt_.dimensions()/dimTime, 0)
   ),

    L_
    (
        IOobject
        (
            "L",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("zero", dimless, 0)
     )

{
    kMin_ = dimensionedScalar("kMin", sqr(dimVelocity), ROOTVSMALL);
    kMax_.readIfPresent(*this);
    omegaMin_ = dimensionedScalar("omegaMin", inv(dimTime), ROOTVSMALL);
    omegaMax_.readIfPresent(*this);

    boundMinMax(kt_, kMin_, kMax_);
    boundMinMax(kl_, kMin_, kMax_);
    boundMinMax(omega_, omegaMin_, omegaMax_);

    nut_ = kt_/(omega_ + omegaMin_);
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> mykkLOmega::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*(kt_) - nut_*twoSymm(fvc::grad(U_)),
            kt_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> mykkLOmega::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> mykkLOmega::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> mykkLOmega::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool mykkLOmega::read()
{
    if (RASModel::read())
    {
        A0_.readIfPresent(coeffDict());
        As_.readIfPresent(coeffDict());
        Av_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Anat_.readIfPresent(coeffDict());
        Abp_.readIfPresent(coeffDict());
        Ats_.readIfPresent(coeffDict());
        CbpCrit_.readIfPresent(coeffDict());
        Cnc_.readIfPresent(coeffDict());
        CnatCrit_.readIfPresent(coeffDict());
        CnatApg_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CtsCrit_.readIfPresent(coeffDict());
        CtsApg_.readIfPresent(coeffDict());
        CrNat_.readIfPresent(coeffDict());
        C11_.readIfPresent(coeffDict());
        C12_.readIfPresent(coeffDict());
        CR_.readIfPresent(coeffDict());
        CalphaTheta_.readIfPresent(coeffDict());
        Css_.readIfPresent(coeffDict());
        CtauL_.readIfPresent(coeffDict());
        Cw1_.readIfPresent(coeffDict());
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        CwR_.readIfPresent(coeffDict());
        Clambda_.readIfPresent(coeffDict());
        CmuStd_.readIfPresent(coeffDict());
        Prtheta_.readIfPresent(coeffDict());
        Sigmak_.readIfPresent(coeffDict());
        Sigmaw_.readIfPresent(coeffDict());

        medinaPhiBP_.readIfPresent("medinaPhiBP", coeffDict());
        furstPG_.readIfPresent("furstPG", coeffDict());
        CDT_.readIfPresent(coeffDict());
        viscosityRatioLimit_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

  void mykkLOmega::boundMinMax    // taken from foam-3.0
(
    volScalarField& vsf,
    const dimensionedScalar& vsf0,
    const dimensionedScalar& vsf1
) const
{
    scalar minVsf = min(vsf).value();
    scalar maxVsf = max(vsf).value();

    if (minVsf < vsf0.value() || maxVsf > vsf1.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << maxVsf
            << " average: " << gAverage(vsf.internalField())
            << endl;
    }

    if (minVsf < vsf0.value())
    {
        vsf.internalField() = max
        (
            max
            (
                vsf.internalField(),
                fvc::average(max(vsf, vsf0))().internalField()
                *pos(vsf0.value() - vsf.internalField())
            ),
            vsf0.value()
        );
        Info<< "new min: " << gMin(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        vsf.boundaryField() = max(vsf.boundaryField(), vsf0.value());
    }
    
    if (maxVsf > vsf1.value())
    {
        vsf.internalField() = min
        (
            min
            (
                vsf.internalField(),
                fvc::average(min(vsf, vsf1))().internalField()
                *neg(vsf1.value() - vsf.internalField())
                // This is needed when all values are above max
                // HJ, 18/Apr/2009
              + pos(vsf1.value() - vsf.internalField())*vsf1.value()
            ),
            vsf1.value()
        );
        Info<< "new max: " << gMax(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        vsf.boundaryField() = min(vsf.boundaryField(), vsf1.value());
    }
}


void mykkLOmega::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        y_.boundaryField() = max(y_.boundaryField(), VSMALL);
    }

    const volScalarField kT(kt_ + kl_);

    const volScalarField lambdaT( sqrt(kt_)/max(omega_,omegaMin_) );

    const volScalarField lambdaEff(min(Clambda_*y_, lambdaT));

    const volScalarField fw
    (
        pow(
            lambdaEff/max(lambdaT,dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
            2./3.)
    );

    const volTensorField gradU(fvc::grad(U_));

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2.0*magSqr(symm(gradU)));

    ReOmega_ = sqr(y_)*Omega/nu();
    L_ = L(ReOmega_);

    fSS_ = fSS(Omega);
    fv_ = fv(sqr(fw)*kt_/nu()/omega_);
    fINT_ = fINT();
    Cmu_ = Cmu(sqrt(S2));

    const volScalarField ktS(fSS_*fw*kt_);

    const volScalarField nuts
    (
        fv_
        *fINT_
        *Cmu_*sqrt(ktS)*lambdaEff
    );
    Pkt_ = nuts*S2;

    const volScalarField ktL(kt_ - ktS);

    BetaTS_ = BetaTS(ReOmega_, L_);
    fTaul_  = fTaul(lambdaEff, ktL, Omega);

    volScalarField nutl
    (
	 C11_*fTaul_*Omega*sqr(lambdaEff)
          * sqrt(ktL)*lambdaEff/nu()
	 + C12_*BetaTS_*ReOmega_*sqr(y_)*Omega
    );

    if (!furstPG_)
      nutl = min(nutl, 0.5*(kl_ + ktL)/sqrt(S2));

    Pkl_ = nutl*S2;

    if (furstPG_)
      nutl = min(nutl, 0.5*(kl_ + ktL)/sqrt(S2));
    

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv_, ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("SMALL", dimless, ROOTVSMALL);
    gammaBP_ = gammaBP(Omega);
    BetaBP_ = 1.0 - exp(-gammaBP_/Abp_); 
    const volScalarField Rbp
    (
        CR_* BetaBP_ *omega_
	/ max(fw,fwMin)
    );
    RBP_ = Rbp * kl_;

    const volScalarField fNatCrit(1.0 - exp(-Cnc_*sqrt(kl_)*y_/nu()));
    // Natural source term divided by kl_
    gammaNAT_ = gammaNAT(ReOmega_, fNatCrit, L_);
    BetaNAT_ = 1.0 - exp(-gammaNAT_/Anat_);
    const volScalarField Rnat
    (
        CrNat_ * BetaNAT_ * Omega
    );
    RNAT_ = Rnat * kl_;

    // Turbulence specific dissipation rate equation
    fOmega_ = fOmega(lambdaEff, lambdaT);
    omega_.boundaryField().updateCoeffs();
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian
        (
            DomegaEff(alphaTEff),
            omega_,
            "laplacian(alphaTEff,omega)"
        )
     ==
        Cw1_*Pkt_*omega_/(kt_ + kMin_)
      - fvm::SuSp
        (
            (1.0 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)/(kt_ + kMin_)
          , omega_
        )
        - fvm::Sp(Cw2_*omega_*sqr(fw), omega_)
      + Cw3_*fOmega_*alphaTEff*sqr(fw)*sqrt(kt_)/pow3(y_)
    );


    omegaEqn().relax();
    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    boundMinMax(omega_, omegaMin_, omegaMax_);


    // Laminar kinetic energy equation
    const volScalarField Dl(CDT_ * nu() * magSqr(fvc::grad(sqrt(kl_))));
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(kl_)
      + fvm::div(phi_, kl_)
      - fvm::Sp(fvc::div(phi_), kl_)
      - fvm::laplacian(nu(), kl_, "laplacian(nu,kl)")
     ==
        Pkl_
        - fvm::Sp(Rbp, kl_)
        - fvm::Sp(Rnat, kl_)
        - fvm::Sp(Dl/max(kl_,kMin_), kl_)
    );

    klEqn().relax();
    klEqn().boundaryManipulate(kl_.boundaryField());

    solve(klEqn);
    boundMinMax(kl_, kMin_, kMax_);

    // Turbulent kinetic energy equation
    const volScalarField Dt(CDT_ * nu() * magSqr(fvc::grad(sqrt(kt_))));
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(kt_)
      + fvm::div(phi_, kt_)
      - fvm::Sp(fvc::div(phi_), kt_)
      - fvm::laplacian(DkEff(alphaTEff), kt_, "laplacian(alphaTEff,kt)")
     ==
        Pkt_
        + (Rbp + Rnat)*kl_
        - fvm::Sp(Dt/max(kt_,kMin_), kt_)
        - fvm::Sp(omega_, kt_)
    );

    ktEqn().relax();
    ktEqn().boundaryManipulate(kt_.boundaryField());

    solve(ktEqn);
    boundMinMax(kt_, kMin_, kMax_);



    // Re-calculate viscosity
    nut_ = nuts + nutl;
    nut_.correctBoundaryConditions();
    dimensionedScalar nuRef("nuRef", dimensionSet(0,2,-1,0,0,0,0), nu()()[0]);
    boundMinMax(nut_, 0.0*nuRef, viscosityRatioLimit_*nuRef);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
