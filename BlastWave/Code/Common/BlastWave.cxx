#include "BlastWave.h"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

BlastWave::BlastWave(double aT,
		     double aRho0,
		     double aR, double aAs,
		     double aTau, double aDeltat):mPar(0),mStat(0){
  setParameters(aT,aRho0,aR,aAs,aTau,aDeltat);
}

BlastWave::BlastWave(double aT,
		     double aRho0, double aRho2, double aRho4, 
		     double aRx, double aRy, double aAs,
		     double aTau, double aDeltat):mPar(0),mStat(0){
  setParameters(aT,aRho0,aRho2,aRho4,aRx,aRy,aAs,aTau,aDeltat);
}

void BlastWave::setParameters(double aT,
			      double aRho0,
			      double aR, double aAs,
			      double aTau, double aDeltat){
  if(!mPar) mPar = new BWParameters;
  mPar->T = aT;
  mPar->Rho0 = aRho0;
  mPar->Rho2 = 0.;
  mPar->Rho4 = 0.;
  mPar->Rx = aR;
  mPar->Ry = aR;
  mPar->As = aAs;
  mPar->Tau = aTau;
  mPar->Deltat = aDeltat;
  mNewPar=1;
  mNewHisto=1;
  mNewAvSepHisto=1;
  mNewAvSepVsPhiPHisto=1;  
  mNewPhiPDepHisto=1;
  mUseR=0; // cannot be set to 1 here. Need to fix minuit Ry as well.
  mSpace =0;
}

void BlastWave::setParameters(double aT,
			      double aRho0, double aRho2, double aRho4,
			      double aRx, double aRy, double aAs,
			      double aTau, double aDeltat){
  if(!mPar) mPar = new BWParameters;
  mPar->T = aT;
  mPar->Rho0 = aRho0;
  mPar->Rho2 = aRho2;
  mPar->Rho4 = aRho4;
  mPar->Rx = aRx;
  mPar->Ry = aRy;
  mPar->As = aAs;
  mPar->Tau = aTau;
  mPar->Deltat = aDeltat;
  mNewPar = 1;
  mNewHisto=1;
  mNewAvSepHisto=1;
  mNewAvSepVsPhiPHisto=1;  
  mNewPhiPDepHisto=1;
  mUseR=0;
}
	
void BlastWave::calc(double aMass, double aPt){
  if(mPar->Ry!=mPar->Rx || mPar->Rho2 !=0 || mPar->Rho4!=0){
    int tNBinPhi = 24;
    double tPhiP;
    double tdNdptpt = 0.;
    double tV2 = 0.;
    double tV4 = 0.;
    double tV6 = 0.;
    double tV8 = 0.;
    double tROutSquare = 0.;
    double tRSideSquare = 0.;
    double tROutSideSquare = 0.;
    double tRLongSquare = 0.;
    double tROutSquareOsc = 0.;
    double tRSideSquareOsc = 0.;
    double tROutSideSquareOsc = 0.;
    double tRLongSquareOsc = 0.;
    double tAvROutStar = 0.;
    double tAvRSide = 0.;
    double tAvRLongStar = 0.;
    double tCos2PhiP = 0.;
    for(int ti=0; ti<tNBinPhi; ti++){
      tPhiP = ti*1./tNBinPhi*TMath::Pi()*2.;
      calc(aMass,aPt,tPhiP);
      tCos2PhiP = cos(2*tPhiP);
      tV2 += ( mResult.dNdptpt  * tCos2PhiP);
      tV4 += ( mResult.dNdptpt  * cos(4.*tPhiP));
      tV6 += ( mResult.dNdptpt  * cos(6.*tPhiP));
      tV8 += ( mResult.dNdptpt  * cos(8.*tPhiP));
      tdNdptpt += mResult.dNdptpt;
      tROutSquare += mResult.ROutSquare;
      tRSideSquare += mResult.RSideSquare;
      tROutSideSquare += mResult.ROutSideSquare;
      tRLongSquare += mResult.RLongSquare;
      tROutSquareOsc += (mResult.ROutSquare * tCos2PhiP);
      tRSideSquareOsc += (mResult.RSideSquare* tCos2PhiP);
      tROutSideSquareOsc += (mResult.ROutSideSquare* sin(2*tPhiP));
      tRLongSquareOsc += (mResult.RLongSquare* tCos2PhiP);
      tAvROutStar += mResult.AvROutStar;
      tAvRSide += mResult.AvRSide;
      tAvRLongStar += mResult.AvRLongStar;
      //      cout << tPhiP << " " << mResult.dNdptpt << endl;
    }
    mResult.dNdptpt = tdNdptpt/tNBinPhi;
    mResult.V2 = tV2/tdNdptpt;
    mResult.V4 = tV4/tdNdptpt;
    mResult.V6 = tV6/tdNdptpt;
    mResult.V8 = tV8/tdNdptpt;
    mResult.ROutSquare = tROutSquare/tNBinPhi;
    mResult.RSideSquare = tRSideSquare/tNBinPhi;
    mResult.ROutSideSquare = tROutSideSquare/tNBinPhi;
    mResult.RLongSquare = tRLongSquare/tNBinPhi;

    mResult.ROutSquareOsc = tROutSquareOsc/tNBinPhi;

    mResult.RSideSquareOsc = tRSideSquareOsc/tNBinPhi;
    mResult.ROutSideSquareOsc = tROutSideSquareOsc/tNBinPhi;
    mResult.RLongSquareOsc = tRLongSquareOsc/tNBinPhi;
    mResult.AvROutStar = tAvROutStar/tNBinPhi;
    mResult.AvRSide = tAvRSide/tNBinPhi;
    mResult.AvRLongStar = tAvRLongStar/tNBinPhi;
  }
  else{
    calc(aMass,aPt,0.);
    mResult.V2=0.;
    mResult.V4=0.;
    mResult.V6=0.;
    mResult.V8=0.;
    mResult.ROutSquareOsc = 0.;
    mResult.RSideSquareOsc = 0.;
    mResult.ROutSideSquareOsc = 0.;
    mResult.RLongSquareOsc = 0.;
  }

}

void BlastWave::calcR(double aMass, double aPt, double aPhiP){
  int tNBinPhiS = 48;
  double tPhiSStep = TMath::Pi()*2./tNBinPhiS;
  double tPhiS;
  double tRTilda;
  double tNBinR = 50;
  double tRMax = mPar->Rx>mPar->Ry?  
    mPar->Rx*(1.+10.*mPar->As) : 
    mPar->Ry*(1.+10.*mPar->As);
  double tRStep = tRMax/1./tNBinR;
  double tRMin = tRStep/2.;

  double tTemp = fabs(mPar->T);

   double tMt = sqrt(aPt*aPt+aMass*aMass);
   double tSum = 0.;
   double tX = 0.;
   double tX2 = 0.;
   double tY = 0.;
   double tY2 = 0.;
   double tXY = 0.;
   double tZ2 = 0.;
   double tT = 0.;
   double tT2 = 0.;
   double tXT = 0.;
   double tYT = 0.;
   double tR;
   double tPhiB;
   double tRho; 
   double tBeta; 
   double tK0Beta; 
   double tK1Beta; 
   double tK2Beta; 
   double tBetaWeight; 
   double tAlphaWeight;
   double tH1OverH0 = (mPar->Deltat*mPar->Deltat+
		       mPar->Tau*mPar->Tau)/mPar->Tau;
   double tH2OverH0 = (3*mPar->Deltat*mPar->Deltat+mPar->Tau*mPar->Tau);
   double tSpaceWeight;

   double tRXSquare = mPar->Rx*mPar->Rx;
   double tRYSquare = mPar->Ry*mPar->Ry;
   double tRxOverRySquare = tRXSquare/tRYSquare;
   double tCosPhiS;
   double tSinPhiS;
   double tEdgeSize;
   double tLX;
   double tLY;

   int tNStat = TMath::Abs(mStat)*3+1; // 4 bins
   double tStatSign = 1.;
   double tStatSignCoef = mStat<0.? -1. : 1.;

   //   if(fabs(aPt-0.990678)<0.001 && aPhiP==0.) cout << " >>> " << aPhiP << " " << aPt << " ";
   for(int tStat=1; tStat<=tNStat; tStat++){ 
     tStatSign *= tStatSignCoef;
     for(double tR=tRMin; tR<=tRMax; tR+=tRStep){
       for(int tiPhiS=0; tiPhiS<tNBinPhiS; tiPhiS++){
	 tPhiS = tiPhiS*tPhiSStep;
	 tCosPhiS = TMath::Cos(tPhiS);
	 tSinPhiS = TMath::Sin(tPhiS);
	 tLX = tCosPhiS*tR;
	 tLY = tSinPhiS*tR;
	 tRTilda = sqrt(tLX*tLX/tRXSquare+tLY*tLY/tRYSquare);
	 tSpaceWeight = mPar->As==0.? tRTilda<1. : 
	   1./(1.+exp((tRTilda-1.)/mPar->As));
	 tPhiB = TMath::ATan2(tSinPhiS*tRxOverRySquare,tCosPhiS);
	 tPhiB = tPhiB<0.? 2.*TMath::Pi()+tPhiB : tPhiB;
	 tRho = tRTilda*(mPar->Rho0+
			 mPar->Rho2*cos(2.*tPhiB)+
			 mPar->Rho4*cos(4.*tPhiB));	
	 tBeta = tMt/tTemp*TMath::CosH(tRho)*tStat;
	 tK1Beta = TMath::BesselK1(tBeta);
	 tBetaWeight = 2.*tK1Beta; 
	 tAlphaWeight = exp(tStat*aPt/tTemp*TMath::SinH(tRho)*
			    cos(aPhiP-tPhiB))*tSpaceWeight*tR;
	 tSum += (tAlphaWeight*tBetaWeight);
	 if(mSpace){
	   tK0Beta = TMath::BesselK0(tBeta);
	   tK2Beta = 2./tBeta*tK1Beta+tK0Beta;
	   tX   += (tAlphaWeight*tBetaWeight * tLX);
	   tX2  += (tAlphaWeight*tBetaWeight * tLX * tLX);
	   tY   += (tAlphaWeight*tBetaWeight * tLY);
	   tY2  += (tAlphaWeight*tBetaWeight * tLY * tLY);
	   tXY  += (tAlphaWeight*tBetaWeight * tLX * tLY);
	   tBetaWeight = 2.*(tK1Beta/tBeta+tK0Beta);
	   tT   += (tAlphaWeight*tBetaWeight * tH1OverH0);
	   tXT  += (tAlphaWeight*tBetaWeight * tLX * tH1OverH0);
	   tYT  += (tAlphaWeight*tBetaWeight * tLY * tH1OverH0);
	   tBetaWeight = 2./tBeta*tK2Beta;
	   tZ2  += (tAlphaWeight*tBetaWeight * tH2OverH0); 
	   tBetaWeight += (2*tK1Beta);
	   tT2  += (tAlphaWeight*tBetaWeight * tH2OverH0);
	 }
	 //	 if(fabs(aPt-0.990678)<0.1 && aPhiP==0.) cout << tPhiS << " " << tR << " " << tRTilda << endl;
       }
     }
   }
   mResult.dNdptpt = tSum*tMt*tStatSign;   
   if(mSpace){
     tX  /= tSum;
     tX2 /= tSum; tX2 -= (tX*tX);
     tY  /= tSum;
     tY2 /= tSum; tY2 -= (tY*tY);
     tXY /= tSum; tXY -= (tX*tY);
     tZ2 /= tSum;
     tT  /= tSum;
     tT2 /= tSum; tT2 -= (tT*tT);
     tXT /= tSum; tXT -= (tX*tT);
     tYT /= tSum; tYT -= (tY*tT);
    
     double tBetaT = aPt/tMt;
     double tGammaT = tMt/aMass;
     double tCosPhiP = cos(aPhiP);
     double tSinPhiP = sin(aPhiP);
     
     mResult .AvROutStarSpace = tGammaT * (tX * tCosPhiP + tY * tSinPhiP);
     mResult .AvROutStarTime = tT * -1. * tGammaT * tBetaT ;
     mResult .AvROutStar =  mResult .AvROutStarSpace + mResult .AvROutStarTime;
     mResult.ROutSquare = 0.5*(tX2+tY2)+0.5*(tX2-tY2)*cos(2*aPhiP)
       + tXY*sin(2*aPhiP)
       - 2*tBetaT*(tXT*tCosPhiP+tYT*tSinPhiP)+tBetaT*tBetaT*tT2;
     mResult.AvRSide = tY * tCosPhiP - tX * tSinPhiP;
     mResult.RSideSquare = 0.5*(tX2+tY2)-0.5*(tX2-tY2)*cos(2*aPhiP)
       -tXY*sin(2*aPhiP);
     mResult.ROutSideSquare = tXY*cos(2*aPhiP)-0.5*(tX2-tY2)*sin(2*aPhiP)
       + tBetaT*(tXT*tSinPhiP-tYT*tCosPhiP);
     mResult.AvRLongStar = 0; // Not implemented yet
     mResult.RLongSquare = tZ2;   
   }
}

// _______________________
void BlastWave::calcRTilda(double aMass, double aPt, double aPhiP){
  int tNBinPhiS = 24;
  double tPhiSStep = TMath::Pi()*2./tNBinPhiS;
  double tPhiS;
  int tNBinRTilda = 20;
  double tRTildaMax =  (1.+10.*mPar->As);
  double tRTildaStep = 1./tNBinRTilda;
  double tRTildaMin = tRTildaStep/2.;
  double tRTilda;
  double tTemp = fabs(mPar->T);
   double tMt = sqrt(aPt*aPt+aMass*aMass);
   double tSum = 0.;
   double tX = 0.;
   double tX2 = 0.;
   double tY = 0.;
   double tY2 = 0.;
   double tXY = 0.;
   double tZ2 = 0.;
   double tT = 0.;
   double tT2 = 0.;
   double tXT = 0.;
   double tYT = 0.;
   double tR;
   double tPhiB;
   double tRho; 
   double tBeta; 
   double tK0Beta; 
   double tK1Beta; 
   double tK2Beta; 
   double tBetaWeight; 
   double tAlphaWeight;
   double tWeight;
   double tH1OverH0 = (mPar->Deltat*mPar->Deltat+
		       mPar->Tau*mPar->Tau)/mPar->Tau;
   double tH2OverH0 = (3*mPar->Deltat*mPar->Deltat+mPar->Tau*mPar->Tau);
   double tSpaceWeight;

   double tRXSquare = mPar->Rx*mPar->Rx;
   double tRYSquare = mPar->Ry*mPar->Ry;
   double tRxOverRySquare = tRXSquare/tRYSquare;
   double tCosPhiS;
   double tSinPhiS;
   double tEdgeSize;
   double tLX;
   double tLY;

   int tNStat = TMath::Abs(mStat)*3+1; // 4 bins
   double tStatSign = 1.;
   double tStatSignCoef = mStat<0.? -1. : 1.;
   for(int tStat=1; tStat<=tNStat; tStat++){ 
     tStatSign *= tStatSignCoef;
     for(double tRTilda=tRTildaMin;
	 tRTilda<=tRTildaMax;
	 tRTilda+=tRTildaStep){
       tSpaceWeight = (mPar->As==0.? 
	 ((tRTilda<1.) ? 1. : 0.) : 
	 1./(1.+exp((tRTilda-1.)/mPar->As)))*tRTildaStep;
       if(tSpaceWeight>0.){
	 tSpaceWeight *= tStatSign;
	 for(int tiPhiS=0; tiPhiS<tNBinPhiS; tiPhiS++){
	   tPhiS = tiPhiS*tPhiSStep;
	   tCosPhiS = TMath::Cos(tPhiS);
	   tSinPhiS = TMath::Sin(tPhiS);
	   tSpaceWeight = mPar->As==0.? tRTilda<1. : 
	     1./(1.+exp((tRTilda-1.)/mPar->As));
	   tEdgeSize = mPar->Rx*mPar->Ry/sqrt(tCosPhiS*tCosPhiS*tRYSquare+
					      tSinPhiS*tSinPhiS*tRXSquare);
	   tR = tRTilda*tEdgeSize;
	   tPhiB = TMath::ATan2(tSinPhiS*tRxOverRySquare,tCosPhiS);
	   tPhiB = tPhiB<0.? 2.*TMath::Pi()+tPhiB : tPhiB;
	   tRho = tRTilda*(mPar->Rho0+
			   mPar->Rho2*cos(2.*tPhiB)+
			   mPar->Rho4*cos(4.*tPhiB));	
	   tBeta = tMt/tTemp*TMath::CosH(tRho)*tStat;
	   tK1Beta = TMath::BesselK1(tBeta);
	   tBetaWeight = 2.*tK1Beta; 
	   tAlphaWeight = exp(tStat*aPt/tTemp*TMath::SinH(tRho)*	
			      cos(aPhiP-tPhiB))*tSpaceWeight*tR*tEdgeSize;
	   tWeight = tAlphaWeight*tBetaWeight*tPhiSStep;
	   tSum += tWeight;
	   if(mSpace){
	     tK0Beta = TMath::BesselK0(tBeta);
	     tK2Beta = 2./tBeta*tK1Beta+tK0Beta;
	     tLX = tCosPhiS*tR;
	     tLY = tSinPhiS*tR;
	     tX   += (tWeight * tLX);
	     tX2  += (tWeight * tLX * tLX);
	     tY   += (tWeight * tLY);
	     tY2  += (tWeight * tLY * tLY);
	     tXY  += (tWeight * tLX * tLY);
	     tBetaWeight = 2.*(tK1Beta/tBeta+tK0Beta);
	     tWeight = tAlphaWeight*tBetaWeight;//(1./tBeta+tK0Beta/tK1Beta);
	     tT   += (tWeight * tH1OverH0);
	     tXT  += (tWeight * tLX * tH1OverH0);
	     tYT  += (tWeight * tLY * tH1OverH0);
	     tBetaWeight = 2./tBeta*tK2Beta;
	     tZ2  += (tAlphaWeight*tBetaWeight * tH2OverH0); 
	     tBetaWeight += (2*tK1Beta);
	     tT2  += (tAlphaWeight*tBetaWeight * tH2OverH0);
	   }	
	 }
       }
     }
   }
   mResult.dNdptpt = tSum*tMt;
   
   if(mSpace){
     tX  /= tSum;
     tX2 /= tSum; tX2 -= (tX*tX);
     tY  /= tSum;
     tY2 /= tSum; tY2 -= (tY*tY);
     tXY /= tSum; tXY -= (tX*tY);
     tZ2 /= tSum;
     tT  /= tSum;
     tT2 /= tSum; tT2 -= (tT*tT);
     tXT /= tSum; tXT -= (tX*tT);
     tYT /= tSum; tYT -= (tY*tT);
    
     double tBetaT = aPt/tMt;
     double tGammaT = tMt/aMass;
     double tCosPhiP = cos(aPhiP);
     double tSinPhiP = sin(aPhiP);
     
     mResult .AvROutStarSpace = tGammaT * (tX * tCosPhiP + tY * tSinPhiP);
     mResult .AvROutStarTime = tT * -1. * tGammaT * tBetaT ;
     mResult .AvROutStar =  mResult .AvROutStarSpace + mResult .AvROutStarTime;
     mResult.ROutSquare = 0.5*(tX2+tY2)+0.5*(tX2-tY2)*cos(2*aPhiP)
       + tXY*sin(2*aPhiP)
       - 2*tBetaT*(tXT*tCosPhiP+tYT*tSinPhiP)+tBetaT*tBetaT*tT2;
     mResult.AvRSide = tY * tCosPhiP - tX * tSinPhiP;
     mResult.RSideSquare = 0.5*(tX2+tY2)-0.5*(tX2-tY2)*cos(2*aPhiP)
       -tXY*sin(2*aPhiP);
     mResult.ROutSideSquare = tXY*cos(2*aPhiP)-0.5*(tX2-tY2)*sin(2*aPhiP)
       + tBetaT*(tXT*tSinPhiP-tYT*tCosPhiP);
     mResult.AvRLongStar = 0; // Not implemented yet
     mResult.RLongSquare = tZ2;   
   }
}
// _______________________
void BlastWave::calcXY(double aMass, double aPt, double aPhiP){
  int tNBinX = 100;
  double tXMax = mPar->Rx*(1.+10.*mPar->As);
  double tXStep = tXMax/1./tNBinX;
  double tXMin = -1.*tXMax;
  int tNBinY = 100;
  double tYMax = mPar->Ry*(1.+10.*mPar->As);
  double tYStep = tYMax/1./tNBinY;
  double tYMin = -1.*tYMax;
  double tR;
  double tRTilda;
  double tTemp = fabs(mPar->T);

   double tMt = sqrt(aPt*aPt+aMass*aMass);
   double tSum = 0.;
   double tX = 0.;
   double tX2 = 0.;
   double tY = 0.;
   double tY2 = 0.;
   double tXY = 0.;
   double tZ2 = 0.;
   double tT = 0.;
   double tT2 = 0.;
   double tXT = 0.;
   double tYT = 0.;

   double tPhiB;
   double tRho; 
   double tBeta; 
   double tK0Beta; 
   double tK1Beta; 
   double tK2Beta; 
   double tBetaWeight; 
   double tAlphaWeight;
   double tH1OverH0 = (mPar->Deltat*mPar->Deltat+
		       mPar->Tau*mPar->Tau)/mPar->Tau;
   double tH2OverH0 = (3*mPar->Deltat*mPar->Deltat+mPar->Tau*mPar->Tau);
   double tSpaceWeight;

   double tRXSquare = mPar->Rx*mPar->Rx;
   double tRYSquare = mPar->Ry*mPar->Ry;
   double tRxOverRySquare = tRXSquare/tRYSquare;
   double tPhiS;
   double tCosPhiS;
   double tSinPhiS;
   double tLX;
   double tLY;

   int tNStat = TMath::Abs(mStat)*3+1; // 4 bins
   double tStatSign = 1.;
   double tStatSignCoef = mStat<0.? -1. : 1.;
   for(int tStat=1; tStat<=tNStat; tStat++){ 
     tStatSign *= tStatSignCoef;
     for(tLX = tXMin; tLX<=tXMax; tLX+=tXStep){
       for(tLY = tYMin; tLY<=tYMax; tLY+=tYStep){
	 tR = sqrt(tLX*tLX+tLY*tLY);
	 tRTilda = sqrt(tLX*tLX/tRXSquare+tLY*tLY/tRYSquare);
	 tSpaceWeight = mPar->As==0.? tRTilda<1. : 
	   1./(1.+exp((tRTilda-1.)/mPar->As));	 
	 tPhiS = TMath::ATan2(tLY,tLX);
	 tCosPhiS = TMath::Cos(tPhiS);
	 tSinPhiS = TMath::Sin(tPhiS);
	 if(tSpaceWeight>0.){
	   tPhiB = TMath::ATan2(tSinPhiS*tRxOverRySquare,tCosPhiS);
	   tPhiB = tPhiB<0.? 2.*TMath::Pi()+tPhiB : tPhiB;
	   tRho = tRTilda*(mPar->Rho0+
			   mPar->Rho2*cos(2.*tPhiB)+
			   mPar->Rho4*cos(4.*tPhiB));	
	   tBeta = tMt/tTemp*TMath::CosH(tRho)*tStat;
	   tK1Beta = TMath::BesselK1(tBeta);
	   tBetaWeight = 2.*tK1Beta; 
	   tAlphaWeight = exp(tStat*aPt/tTemp*TMath::SinH(tRho)*
			      cos(aPhiP-tPhiB))*tSpaceWeight;
	   tSum += (tAlphaWeight*tBetaWeight);
	   if(mSpace){
	     tK0Beta = TMath::BesselK0(tBeta);
	     tK2Beta = 2./tBeta*tK1Beta+tK0Beta;
	     tX   += (tAlphaWeight*tBetaWeight * tLX);
	     tX2  += (tAlphaWeight*tBetaWeight * tLX * tLX);
	     tY   += (tAlphaWeight*tBetaWeight * tLY);
	     tY2  += (tAlphaWeight*tBetaWeight * tLY * tLY);
	     tXY  += (tAlphaWeight*tBetaWeight * tLX * tLY);
	     tBetaWeight = 2.*(tK1Beta/tBeta+tK0Beta);
	     tT   += (tAlphaWeight*tBetaWeight * tH1OverH0);
	     tXT  += (tAlphaWeight*tBetaWeight * tLX * tH1OverH0);
	     tYT  += (tAlphaWeight*tBetaWeight * tLY * tH1OverH0);
	     tBetaWeight = 2./tBeta*tK2Beta;
	     tZ2  += (tAlphaWeight*tBetaWeight * tH2OverH0); 
	     tBetaWeight += (2*tK1Beta);
	     tT2  += (tAlphaWeight*tBetaWeight * tH2OverH0);
	   }
	 }
       }
     }
   }
   mResult.dNdptpt = tSum*tMt*tStatSign;
   
   if(mSpace){
     tX  /= tSum;
     tX2 /= tSum; tX2 -= (tX*tX);
     tY  /= tSum;
     tY2 /= tSum; tY2 -= (tY*tY);
     tXY /= tSum; tXY -= (tX*tY);
     tZ2 /= tSum;
     tT  /= tSum;
     tT2 /= tSum; tT2 -= (tT*tT);
     tXT /= tSum; tXT -= (tX*tT);
     tYT /= tSum; tYT -= (tY*tT);
    
     double tBetaT = aPt/tMt;
     double tGammaT = tMt/aMass;
     double tCosPhiP = cos(aPhiP);
     double tSinPhiP = sin(aPhiP);
     
     mResult .AvROutStarSpace = tGammaT * (tX * tCosPhiP + tY * tSinPhiP);
     mResult .AvROutStarTime = tT * -1. * tGammaT * tBetaT ;
     mResult .AvROutStar =  mResult .AvROutStarSpace + mResult .AvROutStarTime;
     mResult.ROutSquare = 0.5*(tX2+tY2)+0.5*(tX2-tY2)*cos(2*aPhiP)
       + tXY*sin(2*aPhiP)
       - 2*tBetaT*(tXT*tCosPhiP+tYT*tSinPhiP)+tBetaT*tBetaT*tT2;
     mResult.AvRSide = tY * tCosPhiP - tX * tSinPhiP;
     mResult.RSideSquare = 0.5*(tX2+tY2)-0.5*(tX2-tY2)*cos(2*aPhiP)
       -tXY*sin(2*aPhiP);
     mResult.ROutSideSquare = tXY*cos(2*aPhiP)-0.5*(tX2-tY2)*sin(2*aPhiP)
       + tBetaT*(tXT*tSinPhiP-tYT*tCosPhiP);
     mResult.AvRLongStar = 0; // Not implemented yet
     mResult.RLongSquare = tZ2;   
   }
}
//_____

void BlastWave::calcHisto(const char* aNamePref, double aMass){
  if(mNewHisto || mPrevMass!=aMass){
    mNewHisto=0;
    mPrevMass = aMass;
    double tNorm=0.;
    int tNBin =59;
    double tPtMin =0.05;
    double tPtMax =3.05;
    char tHSpectraName[50];
    sprintf(tHSpectraName,"%sSpectra",aNamePref);
    mHSpectra = (TH1D*) gROOT->FindObject(tHSpectraName);
    if(mHSpectra){
      mHSpectra->Reset("ICE");
    }
    else{
      mHSpectra = new TH1D(tHSpectraName,tHSpectraName,
			   tNBin,tPtMin,tPtMax);
    }

    char tHV2Name[50];
    sprintf(tHV2Name,"%sV2",aNamePref);
    mHV2 = (TH1D*) gROOT->FindObject(tHV2Name);
    if(mHV2){
      mHV2->Reset("ICE");
    }
    else{
      mHV2 = new TH1D(tHV2Name,tHV2Name,
		      tNBin,tPtMin,tPtMax);
    }

    char tHV4Name[50];
    sprintf(tHV4Name,"%sV4",aNamePref);
    mHV4 = (TH1D*) gROOT->FindObject(tHV4Name);
    if(mHV4){
      mHV4->Reset("ICE");
    }
    else{
      mHV4 = new TH1D(tHV4Name,tHV4Name,
		      tNBin,tPtMin,tPtMax);
    }

    char tHV6Name[50];
    sprintf(tHV6Name,"%sV6",aNamePref);
    mHV6 = (TH1D*) gROOT->FindObject(tHV6Name);
    if(mHV6){
      mHV6->Reset("ICE");
    }
    else{
      mHV6 = new TH1D(tHV6Name,tHV6Name,
		      tNBin,tPtMin,tPtMax);
    }

    char tHV8Name[50];
    sprintf(tHV8Name,"%sV8",aNamePref);
    mHV8 = (TH1D*) gROOT->FindObject(tHV8Name);
    if(mHV8){
      mHV8->Reset("ICE");
    }
    else{
      mHV8 = new TH1D(tHV8Name,tHV8Name,
		      tNBin,tPtMin,tPtMax);
    }

    char tHROutName[50];
    sprintf(tHROutName,"%sROut",aNamePref);
    mHROut = (TH1D*) gROOT->FindObject(tHROutName);
    if(mHROut){
      mHROut->Reset("ICE");
    }
    else{
      mHROut = new TH1D(tHROutName,tHROutName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRSideName[50];
    sprintf(tHRSideName,"%sRSide",aNamePref);
    mHRSide = (TH1D*) gROOT->FindObject(tHRSideName);
    if(mHRSide){
      mHRSide->Reset("ICE");
    }
    else{
      mHRSide = new TH1D(tHRSideName,tHRSideName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRLongName[50];
    sprintf(tHRLongName,"%sRLong",aNamePref);
    mHRLong = (TH1D*) gROOT->FindObject(tHRLongName);
    if(mHRLong){
      mHRLong->Reset("ICE");
    }
    else{
      mHRLong = new TH1D(tHRLongName,tHRLongName,
			tNBin,tPtMin,tPtMax);
    }


    double tMtMin = sqrt(aMass*aMass+mHSpectra->GetBinLowEdge(1)*
			 mHSpectra->GetBinLowEdge(1));
    double tMtMax =  sqrt(aMass*aMass+
			  mHSpectra->GetBinLowEdge(mHSpectra->GetNbinsX()+1)*
			  mHSpectra->GetBinLowEdge(mHSpectra->GetNbinsX()+1));

    char tHROutVsMtName[50];
    sprintf(tHROutVsMtName,"%sROutVsMt",aNamePref);
    mHROutVsMt = (TH1D*) gROOT->FindObject(tHROutVsMtName);
    if(mHROutVsMt){
      mHROutVsMt->Reset("ICE");
    }
    else{
      mHROutVsMt = new TH1D(tHROutVsMtName,tHROutVsMtName,
			tNBin,tMtMin,tMtMax);
    }

    char tHRSideVsMtName[50];
    sprintf(tHRSideVsMtName,"%sRSideVsMt",aNamePref);
    mHRSideVsMt = (TH1D*) gROOT->FindObject(tHRSideVsMtName);
    if(mHRSideVsMt){
      mHRSideVsMt->Reset("ICE");
    }
    else{
      mHRSideVsMt = new TH1D(tHRSideVsMtName,tHRSideVsMtName,
			tNBin,tMtMin,tMtMax);
    }

    char tHRLongVsMtName[50];
    sprintf(tHRLongVsMtName,"%sRLongVsMt",aNamePref);
    mHRLongVsMt = (TH1D*) gROOT->FindObject(tHRLongVsMtName);
    if(mHRLongVsMt){
      mHRLongVsMt->Reset("ICE");
    }
    else{
      mHRLongVsMt = new TH1D(tHRLongVsMtName,tHRLongVsMtName,
			tNBin,tMtMin,tMtMax);
    }

    char tHROutSquareName[50];
    sprintf(tHROutSquareName,"%sROut2",aNamePref);
    mHROutSquare = (TH1D*) gROOT->FindObject(tHROutSquareName);
    if(mHROutSquare){
      mHROutSquare->Reset("ICE");
    }
    else{
      mHROutSquare = new TH1D(tHROutSquareName,tHROutSquareName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRSideSquareName[50];
    sprintf(tHRSideSquareName,"%sRSide2",aNamePref);
    mHRSideSquare = (TH1D*) gROOT->FindObject(tHRSideSquareName);
    if(mHRSideSquare){
      mHRSideSquare->Reset("ICE");
    }
    else{
      mHRSideSquare = new TH1D(tHRSideSquareName,tHRSideSquareName,
			tNBin,tPtMin,tPtMax);
    }

    char tHROutSideSquareName[50];
    sprintf(tHROutSideSquareName,"%sROutSide2",aNamePref);
    mHROutSideSquare = (TH1D*) gROOT->FindObject(tHROutSideSquareName);
    if(mHROutSideSquare){
      mHROutSideSquare->Reset("ICE");
    }
    else{
      mHROutSideSquare = new TH1D(tHROutSideSquareName,tHROutSideSquareName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRLongSquareName[50];
    sprintf(tHRLongSquareName,"%sRLong2",aNamePref);
    mHRLongSquare = (TH1D*) gROOT->FindObject(tHRLongSquareName);
    if(mHRLongSquare){
      mHRLongSquare->Reset("ICE");
    }
    else{
      mHRLongSquare = new TH1D(tHRLongSquareName,tHRLongSquareName,
			tNBin,tPtMin,tPtMax);
    }

    char tHROutSquareOscName[50];
    sprintf(tHROutSquareOscName,"%sROut2Osc",aNamePref);
    mHROutSquareOsc = (TH1D*) gROOT->FindObject(tHROutSquareOscName);
    if(mHROutSquareOsc){
      mHROutSquareOsc->Reset("ICE");
    }
    else{
      mHROutSquareOsc = new TH1D(tHROutSquareOscName,tHROutSquareOscName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRSideSquareOscName[50];
    sprintf(tHRSideSquareOscName,"%sRSide2Osc",aNamePref);
    mHRSideSquareOsc = (TH1D*) gROOT->FindObject(tHRSideSquareOscName);
    if(mHRSideSquareOsc){
      mHRSideSquareOsc->Reset("ICE");
    }
    else{
      mHRSideSquareOsc = new TH1D(tHRSideSquareOscName,tHRSideSquareOscName,
			tNBin,tPtMin,tPtMax);
    }

    char tHROutSideSquareOscName[50];
    sprintf(tHROutSideSquareOscName,"%sROutSide2Osc",aNamePref);
    mHROutSideSquareOsc = (TH1D*) gROOT->FindObject(tHROutSideSquareOscName);
    if(mHROutSideSquareOsc){
      mHROutSideSquareOsc->Reset("ICE");
    }
    else{
      mHROutSideSquareOsc = new TH1D(tHROutSideSquareOscName,tHROutSideSquareOscName,
			tNBin,tPtMin,tPtMax);
    }

    char tHRLongSquareOscName[50];
    sprintf(tHRLongSquareOscName,"%sRLong2Osc",aNamePref);
    mHRLongSquareOsc = (TH1D*) gROOT->FindObject(tHRLongSquareOscName);
    if(mHRLongSquareOsc){
      mHRLongSquareOsc->Reset("ICE");
    }
    else{
      mHRLongSquareOsc = new TH1D(tHRLongSquareOscName,tHRLongSquareOscName,
			tNBin,tPtMin,tPtMax);
    }

    char tHAvOutName[50];
    sprintf(tHAvOutName,"%sHAvOut",aNamePref);
    mHAvOut = (TH1D*) gROOT->FindObject(tHAvOutName);
    if(mHAvOut){
      mHAvOut->Reset("ICE");
    }
    else{
      mHAvOut = new TH1D(tHAvOutName,tHAvOutName,
			tNBin,tPtMin,tPtMax);
    }

    char tHAvTimeName[50];
    sprintf(tHAvTimeName,"%sHAvTime",aNamePref);
    mHAvTime = (TH1D*) gROOT->FindObject(tHAvTimeName);
    if(mHAvTime){
      mHAvTime->Reset("ICE");
    }
    else{
      mHAvTime = new TH1D(tHAvTimeName,tHAvTimeName,
			tNBin,tPtMin,tPtMax);
    }

    char tHAvSideName[50];
    sprintf(tHAvSideName,"%sHAvSide",aNamePref);
    mHAvSide = (TH1D*) gROOT->FindObject(tHAvSideName);
    if(mHAvSide){
      mHAvSide->Reset("ICE");
    }
    else{
      mHAvSide = new TH1D(tHAvSideName,tHAvSideName,
			tNBin,tPtMin,tPtMax);
    }

    mSpace = 1;
    for(int ti=1; ti<=mHSpectra->GetNbinsX(); ti++){
      calc(aMass, mHSpectra->GetBinCenter(ti));    
      mHSpectra->SetBinContent(ti,mResult.dNdptpt);
      tNorm += mResult.dNdptpt;    
      mHV2->SetBinContent(ti, mResult.V2);
      mHV4->SetBinContent(ti, mResult.V4);
      mHV6->SetBinContent(ti, mResult.V6);
      mHV8->SetBinContent(ti, mResult.V8);
      mHROut->SetBinContent(ti, sqrt(mResult.ROutSquare));
      mHRSide->SetBinContent(ti, sqrt(mResult.RSideSquare));
      mHRLong->SetBinContent(ti, sqrt(mResult.RLongSquare));
      mHROutVsMt->SetBinContent(ti, sqrt(mResult.ROutSquare));
      mHRSideVsMt->SetBinContent(ti, sqrt(mResult.RSideSquare));
      mHRLongVsMt->SetBinContent(ti, sqrt(mResult.RLongSquare));
      mHROutSquare->SetBinContent(ti, mResult.ROutSquare);
      mHRSideSquare->SetBinContent(ti, mResult.RSideSquare);
      mHRLongSquare->SetBinContent(ti, mResult.RLongSquare);
      mHROutSideSquare->SetBinContent(ti, mResult.ROutSideSquare);
      mHROutSquareOsc->SetBinContent(ti, mResult.ROutSquareOsc);
      mHRSideSquareOsc->SetBinContent(ti, mResult.RSideSquareOsc);
      mHRLongSquareOsc->SetBinContent(ti, mResult.RLongSquareOsc);
      mHROutSideSquareOsc->SetBinContent(ti, mResult.ROutSideSquareOsc);

      double tPt = mHSpectra->GetBinCenter(ti);
      double tMt =sqrt(aMass*aMass+tPt*tPt);
      double tBetaT = tPt/tMt;
      double tGammaT = tMt/aMass;      
      mHAvOut->SetBinContent(ti,mResult .AvROutStarSpace/tGammaT);
      mHAvTime->SetBinContent(ti,-1.*mResult .AvROutStarTime/tGammaT/tBetaT);
      mHAvSide->SetBinContent(ti,mResult .AvRSide);								      
    }
    mHSpectra->Scale(1./tNorm);
  }
}

TH1D* BlastWave::spectra(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHSpectra; 
}
TH1D* BlastWave::rOut(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROut; 
}
TH1D* BlastWave::rSide(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRSide; 
}
TH1D* BlastWave::rLong(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRLong; 
}
TH1D* BlastWave::rOutVsMt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROutVsMt; 
}
TH1D* BlastWave::rSideVsMt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRSideVsMt; 
}
TH1D* BlastWave::rLongVsMt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRLongVsMt; 
}
TH1D* BlastWave::rOutSquare(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROutSquare; 
}
TH1D* BlastWave::rSideSquare(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRSideSquare; 
}
TH1D* BlastWave::rLongSquare(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRLongSquare; 
}
TH1D* BlastWave::rOutSideSquare(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROutSideSquare; 
}
TH1D* BlastWave::rOutSquareOsc(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROutSquareOsc; 
}
TH1D* BlastWave::rSideSquareOsc(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRSideSquareOsc; 
}
TH1D* BlastWave::rLongSquareOsc(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHRLongSquareOsc; 
}
TH1D* BlastWave::rOutSideSquareOsc(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHROutSideSquareOsc; 
}

TH1D* BlastWave::v2VsPt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHV2; 
}
TH1D* BlastWave::v4VsPt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHV4; 
}
TH1D* BlastWave::v6VsPt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHV6; 
}
TH1D* BlastWave::v8VsPt(const char* aNamePref, double aMass){
  calcHisto(aNamePref, aMass);
  return mHV8; 
}
TH1D* BlastWave::avOut(const char* aNamePref, double aMass){
	calcHisto(aNamePref, aMass);		
	return mHAvOut;
}
TH1D* BlastWave::avTime(const char* aNamePref, double aMass){
	calcHisto(aNamePref, aMass);		
	return mHAvTime;	
}
TH1D* BlastWave::avSide(const char* aNamePref, double aMass){
	calcHisto(aNamePref, aMass);		
	return mHAvSide;	
}


void BlastWave::calcAvSepHisto(const char* aNamePref, double aMassLight, 
			       double aMassHeavy){
  if(aMassLight>aMassHeavy){
    double tMass = aMassLight;
    aMassLight = aMassHeavy;
    aMassHeavy = tMass;
  }
  static double tPrevMassLight;
  static double tPrevMassHeavy;
  if(mNewAvSepHisto || aMassLight!=tPrevMassLight || 
     aMassHeavy!=tPrevMassHeavy){
    mNewAvSepHisto=0;
    tPrevMassLight = aMassLight;
    tPrevMassHeavy = aMassHeavy;
    char tHOutName[50];
    sprintf(tHOutName,"%sOut",aNamePref);
    mHAvSepOut = new TH1D(tHOutName,tHOutName,59,0.,1.5);
    char tHOutSpaceName[50];
    sprintf(tHOutSpaceName,"%sOutSpace",aNamePref);
    mHAvSepOutSpace = new TH1D(tHOutSpaceName,tHOutSpaceName,59,0.,1.5);
    char tHOutTimeName[50];
    sprintf(tHOutTimeName,"%sOutTime",aNamePref);
    mHAvSepOutTime = new TH1D(tHOutTimeName,tHOutTimeName,59,0.,1.5);
    char tHSideName[50];
    sprintf(tHSideName,"%sSide",aNamePref);
    mHAvSepSide = new TH1D(tHSideName,tHSideName,59,0.,1.5);
    char tHLongName[50];
    sprintf(tHLongName,"%sLong",aNamePref);
    mHAvSepLong = new TH1D(tHLongName,tHLongName,59,0.,1.5);
    mSpace=1;
    for(int ti=1; ti<=mHAvSepOut->GetNbinsX(); ti++){
      double tPtHeavy = mHAvSepOut->GetBinCenter(ti);      
      calc(aMassLight,aMassLight/aMassHeavy*tPtHeavy);
      double tOut = mResult.AvROutStar;
      double tOutSpace = mResult.AvROutStarSpace;
      double tOutTime = mResult.AvROutStarTime;
      double tSide = mResult.AvRSide;
      double tLong = mResult.AvRLongStar;
      calc(aMassHeavy,tPtHeavy);
      mHAvSepOut->SetBinContent(ti,tOut-mResult.AvROutStar);
      mHAvSepOutSpace->SetBinContent(ti,tOutSpace-mResult.AvROutStarSpace);
      mHAvSepOutTime->SetBinContent(ti,tOutTime-mResult.AvROutStarTime);
      mHAvSepSide->SetBinContent(ti,tSide-mResult.AvRSide);
      mHAvSepLong->SetBinContent(ti,tLong-mResult.AvRLongStar);
    }
  }
}

TH1D* BlastWave::avSepOut(const char* aNamePref, double aMassLight, double aMassHeavy){
  calcAvSepHisto(aNamePref, aMassLight, aMassHeavy);
  return mHAvSepOut;
}
TH1D* BlastWave::avSepOutSpace(const char* aNamePref, double aMassLight, double aMassHeavy){
  calcAvSepHisto(aNamePref, aMassLight, aMassHeavy);
  return mHAvSepOutSpace;
}
TH1D* BlastWave::avSepOutTime(const char* aNamePref, double aMassLight, double aMassHeavy){
  calcAvSepHisto(aNamePref, aMassLight, aMassHeavy);
  return mHAvSepOutTime;
}
TH1D* BlastWave::avSepSide(const char* aNamePref, 
			   double aMassLight, double aMassHeavy){
  calcAvSepHisto(aNamePref, aMassLight, aMassHeavy);
  return mHAvSepSide;
};
TH1D* BlastWave::avSepLong(const char* aNamePref, 
			   double aMassLight, double aMassHeavy){
  calcAvSepHisto(aNamePref, aMassLight, aMassHeavy);
  return mHAvSepLong;
};

void BlastWave::calcAvSepVsPhiPHisto(const char* aNamePref, double aMassLight, 
			       double aMassHeavy, double aPtHeavy){
  if(aMassLight>aMassHeavy){
    double tMass = aMassLight;
    aMassLight = aMassHeavy;
    aMassHeavy = tMass;
  }
  static double tPrevMassLight;
  static double tPrevMassHeavy;
    static double tPrevPtHeavy;
  if(mNewAvSepVsPhiPHisto || aMassLight!=tPrevMassLight || 
     aMassHeavy!=tPrevMassHeavy || aPtHeavy!=tPrevPtHeavy){
    mNewAvSepVsPhiPHisto=0;
    tPrevMassLight = aMassLight;
    tPrevMassHeavy = aMassHeavy;
    tPrevPtHeavy = aPtHeavy;
    char tHOutName[50];
    sprintf(tHOutName,"%sOut",aNamePref);
    mHAvSepOutVsPhiP = new TH1D(tHOutName,tHOutName,29,0.,2*TMath::Pi());
    char tHOutSpaceName[50];
    sprintf(tHOutSpaceName,"%sOutSpace",aNamePref);
    mHAvSepOutSpaceVsPhiP = new TH1D(tHOutSpaceName,tHOutSpaceName,29,0.,
				     2*TMath::Pi());
    char tHOutTimeName[50];
    sprintf(tHOutTimeName,"%sOutTime",aNamePref);
    mHAvSepOutTimeVsPhiP = new TH1D(tHOutTimeName,tHOutTimeName,
				    29,0.,2*TMath::Pi());
    char tHSideName[50];
    sprintf(tHSideName,"%sSide",aNamePref);
    mHAvSepSideVsPhiP = new TH1D(tHSideName,tHSideName,29,0.,2*TMath::Pi());
    char tHLongName[50];
    sprintf(tHLongName,"%sLong",aNamePref);
    mHAvSepLongVsPhiP = new TH1D(tHLongName,tHLongName,29,0.,2*TMath::Pi());
    mSpace=1;
    for(int ti=1; ti<=mHAvSepOutVsPhiP->GetNbinsX(); ti++){
      double tPhiP = mHAvSepOutVsPhiP->GetBinCenter(ti);      
      calc(aMassLight,aMassLight/aMassHeavy*aPtHeavy, tPhiP);
      double tOut = mResult.AvROutStar;
      double tOutSpace = mResult.AvROutStarSpace;
      double tOutTime = mResult.AvROutStarTime;
      double tSide = mResult.AvRSide;
      double tLong = mResult.AvRLongStar;
      calc(aMassHeavy,aPtHeavy, tPhiP);
      mHAvSepOutVsPhiP->SetBinContent(ti,tOut-mResult.AvROutStar);
      mHAvSepOutSpaceVsPhiP->SetBinContent(ti,tOutSpace-mResult.AvROutStarSpace);
      mHAvSepOutTimeVsPhiP->SetBinContent(ti,tOutTime-mResult.AvROutStarTime);
      mHAvSepSideVsPhiP->SetBinContent(ti,tSide-mResult.AvRSide);
      mHAvSepLongVsPhiP->SetBinContent(ti,tLong-mResult.AvRLongStar);
    }
  }
}

TH1D* BlastWave::avSepOutVsPhiP(const char* aNamePref, double aMassLight, double aMassHeavy, double aPtHeavy){
  calcAvSepVsPhiPHisto(aNamePref, aMassLight, aMassHeavy, aPtHeavy);
  return mHAvSepOutVsPhiP;
}
TH1D* BlastWave::avSepOutSpaceVsPhiP(const char* aNamePref, double aMassLight, double aMassHeavy, double aPtHeavy){
  calcAvSepVsPhiPHisto(aNamePref, aMassLight, aMassHeavy, aPtHeavy);
  return mHAvSepOutSpaceVsPhiP;
}
TH1D* BlastWave::avSepOutTimeVsPhiP(const char* aNamePref, double aMassLight, double aMassHeavy, double aPtHeavy){
  calcAvSepVsPhiPHisto(aNamePref, aMassLight, aMassHeavy, aPtHeavy);
  return mHAvSepOutTimeVsPhiP;
}
TH1D* BlastWave::avSepSideVsPhiP(const char* aNamePref, 
			   double aMassLight, double aMassHeavy, double aPtHeavy){
  calcAvSepVsPhiPHisto(aNamePref, aMassLight, aMassHeavy, aPtHeavy);
  return mHAvSepSideVsPhiP;
};
TH1D* BlastWave::avSepLongVsPhiP(const char* aNamePref, 
			   double aMassLight, double aMassHeavy, double aPtHeavy){
  calcAvSepVsPhiPHisto(aNamePref, aMassLight, aMassHeavy, aPtHeavy);
  return mHAvSepLongVsPhiP;
}

void BlastWave::calcPhiPDependentHisto(const char* aNamePref, 
				       double aMass, double aPt){
  static double tPrevPt;
  if(mNewPhiPDepHisto || tPrevPt!=aPt){
    mNewPhiPDepHisto=0;
    tPrevPt = aPt;

    char tHROutName[50];
    sprintf(tHROutName,"%sROut2VsPhiP",aNamePref);
    mHROutSquareVsPhiP = new TH1D(tHROutName,tHROutName,
			    20,0.,180.);
    char tHRSideName[50];
    sprintf(tHRSideName,"%sRSide2VsPhi",aNamePref);
    mHRSideSquareVsPhiP = new TH1D(tHRSideName,tHRSideName,
			     20,0.,180.);
    char tHRLongName[50];
    sprintf(tHRLongName,"%sRLong2VsPhiP",aNamePref);
    mHRLongSquareVsPhiP = new TH1D(tHRLongName,tHRLongName,
			     20,0.,180.);
    char tHROutSideName[50];
    sprintf(tHROutSideName,"%sROutSide2VsPhiP",aNamePref);
    mHROutSideSquareVsPhiP = new TH1D(tHROutSideName,
				      tHROutSideName,
				      20,0.,180.);
    mSpace = 1;
    for(int ti=1; ti<=mHROutSquareVsPhiP->GetNbinsX(); ti++){
      calc(aMass, aPt, mHROutSquareVsPhiP->GetBinCenter(ti)/180.*TMath::Pi());
      mHROutSquareVsPhiP->SetBinContent(ti, mResult.ROutSquare);
      mHRSideSquareVsPhiP->SetBinContent(ti, mResult.RSideSquare);
      mHRLongSquareVsPhiP->SetBinContent(ti, mResult.RLongSquare);
      mHROutSideSquareVsPhiP->SetBinContent(ti, mResult.ROutSideSquare);
    }
  }
}
TH1D* BlastWave::rOutSquareVsPhiP(const char* aNamePref, double aMass, double aPt){
  calcPhiPDependentHisto(aNamePref, aMass, aPt);
  return mHROutSquareVsPhiP;
}
TH1D* BlastWave::rSideSquareVsPhiP(const char* aNamePref, double aMass, double aPt){
  calcPhiPDependentHisto(aNamePref, aMass, aPt);
  return mHRSideSquareVsPhiP;
};
TH1D* BlastWave::rLongSquareVsPhiP(const char* aNamePref, double aMass, double aPt){
  calcPhiPDependentHisto(aNamePref, aMass, aPt);
  return mHRLongSquareVsPhiP;
};
TH1D* BlastWave::rOutSideSquareVsPhiP(const char* aNamePref, double aMass, 
				double aPt){
  calcPhiPDependentHisto(aNamePref, aMass, aPt);
  return mHROutSideSquareVsPhiP;
}



TH2D* BlastWave::emPosDist(double aMass){
  double tPtMax = 2.05;
  double tPtMin = -2.05;
  double tPtStep = 0.1;
  //double tPtStep = 0.1;
  //int tNPtBin = 20;
  //int tNPhiBin = 20;
  //double tPhiStep = 2*TMath::Pi()/tNPhiBin;
  double tStep = 0.1; // fm
  double tMax = mPar->Rx>mPar->Ry? mPar->Rx : mPar->Ry;
  tMax+=mPar->As*50.;
  int tNBinX = (int) floor(2*tMax/tStep) +1;
  int tNBinY = (int) floor(2*tMax/tStep) +1;
  TH2D* tH = new TH2D("emPosDist","emPosDist",
		     tNBinX, -tMax, tMax,
		     tNBinY, -tMax, tMax);
  int tii, tjj;
  double tAmp;
  double tAmpSum=0.;
  double tPx, tPy;
  for(int ti=1; ti <= (tNBinX/2+1); ti++){
    double tX = tH->GetXaxis()->GetBinCenter(ti);
    tii = tNBinX-ti+1;
    for(int tj=1; tj <= (tNBinY/2+1); tj++){
      double tY = tH->GetYaxis()->GetBinCenter(tj);
      double tR = sqrt(tX*tX+tY*tY);
      tjj = tNBinY-tj+1;
      tAmp =0.;
      for(tPx = tPtMin; tPx<tPtMax; tPx+=tPtStep){
	for(tPy = tPtMin; tPy<tPtMax; tPy+=tPtStep){	  
	  tAmp += getEmissionProba(aMass, tX, tY, 
				   sqrt(tPx*tPx+tPy*tPy),
				   TMath::ATan2(tPy,tPx));
	}
      }
      tAmpSum += tAmp;
      tH->SetCellContent(ti,tj,tAmp);
      tH->SetCellContent(ti,tjj,tAmp);
      tH->SetCellContent(tii,tj,tAmp);
      tH->SetCellContent(tii,tjj,tAmp);
    }
  }
  tH->Scale(1./tAmpSum);
  return tH;
}
TH2D* BlastWave::emMomDist(double aMass){
  double tStep = 0.1;
  TH2D* tH = new TH2D("emMomDist","emMomDist",41,-1.025,1.025,41,-1.025,1.025);
  double tAmp;
  double tAmpSum=0.;
  double tNBin =20;
  double tStepX = 2.*mPar->Rx/tNBin;
  double tStepY = 2.*mPar->Ry/tNBin;
  double tMinX = -mPar->Rx+tStepX/2.;
  double tMinY = -mPar->Ry+tStepX/2.;
  for(int ti=1; ti<= (tH->GetNbinsX()/2+1); ti++){
    double tPx = tH->GetXaxis()->GetBinCenter(ti);
    int tii = tH->GetNbinsX()-ti+1;
    for(int tj=1; tj<= (tH->GetNbinsY()/2+1); tj++){
      double tPy = tH->GetYaxis()->GetBinCenter(tj);
      double tPt = sqrt(tPx*tPx+tPy*tPy);
      double tPhi = TMath::ATan2(tPy,tPx);
      int tjj = tH->GetNbinsY()-tj+1;
      tAmp = 0.;
      for(double tX = tMinX; tX<mPar->Rx; tX+= tStepX){
	for(double tY = tMinY; tY<mPar->Ry; tY+= tStepY){
	  tAmp += getEmissionProba(aMass, tX, tY, tPt, tPhi);
	}
      }
      tAmpSum += tAmp;
      tH->SetCellContent(ti,tj,tAmp);
      tH->SetCellContent(ti,tjj,tAmp);
      tH->SetCellContent(tii,tj,tAmp);
      tH->SetCellContent(tii,tjj,tAmp);
    }
  }
  tH->Scale(1./tAmpSum);
  return tH;
}
TH2D* BlastWave::emMomDist(double aMass, double aX, double aY){
  return 0;
}
TH2D* BlastWave::emPosDist(double aMass, double aPt, double aPhiP){
	static int tiHist;
  double tStep = 0.1; // fm
  double tMax = mPar->Rx>mPar->Ry? mPar->Rx : mPar->Ry;
  tMax+=mPar->As*50.+tStep;
  int tNBinX = (int) floor(2*tMax/tStep) +1;
  int tNBinY = (int) floor(2*tMax/tStep) +1;
  char tHName[50];
  sprintf(tHName,"emPosDist%i",tiHist);
  tiHist++;
  TH2D* tH = new TH2D(tHName,tHName,
		     tNBinX, -tMax-tStep/2., tMax+tStep/2.,
		     tNBinY, -tMax-tStep/2., tMax+tStep/2.);
  int tii, tjj;
  double tAmp;
  double tAmpSum=0.;
  double tPx, tPy;
  for(int ti=1; ti <= tNBinX; ti++){
    double tX = tH->GetXaxis()->GetBinCenter(ti);
    tii = tNBinX-ti+1;
    for(int tj=1; tj <= tNBinY; tj++){
      double tY = tH->GetYaxis()->GetBinCenter(tj);
      //double tR = sqrt(tX*tX+tY*tY);
      //tjj = tNBinY-tj+1;
      tAmp =getEmissionProba(aMass, tX, tY, aPt, aPhiP);
      tAmpSum += tAmp;
      tH->SetCellContent(ti,tj,tAmp);
    }
  }
  tH->Scale(1./tAmpSum);
  return tH;
}

//TH2D* TimeVsPtDist(double aMass){
//  
//}
TH1D* BlastWave::ZDist(double aMass, double aPt) {
  static int tiHist;
  char tHName[50];
  sprintf(tHName,"ZDist%i",tiHist);
  tiHist++;
  int tNBinZ=100;
  TH1D* tH = new TH1D(tHName,tHName,tNBinZ,-50.,50.);

  double tStep = 0.5; // fm
  double tMax = mPar->Rx>mPar->Ry? mPar->Rx : mPar->Ry;
  tMax+=mPar->As*50.;
  tMax = floor(tMax/tStep)*tStep;

  double tAmpSum=0.;
  for(int tiZ=1; tiZ<=tNBinZ; tiZ++){
    double tEta = TMath::ASinH(tH->GetBinCenter(tiZ)/mPar->Tau);
    double tAmp = 0.;
    for(double tX=0.; tX <= tMax; tX+=tStep){
      for(double tY=0.; tY <= tMax; tY+=tStep){
	tAmp += getEmissionProba(aMass, tX, tY, aPt, 0, tEta);
      }
    }
    tH->SetBinContent(tiZ,tAmp);
    tAmpSum += tAmp;
  }
  tH->Scale(1./tAmpSum);
  return tH;
}

void BlastWave::coutParameters(){
  cout << "Blast wave paremeters" << endl
       << "T (MeV) = " << mPar->T*1000. << endl
       << "rho0 = " << mPar->Rho0 << endl
       << "rhoa = " << mPar->Rho2 << endl
       << "Rx (fm) = " << mPar->Rx << endl
       << "Ry (fm) = " << mPar->Ry << endl
       << "Tau (fm/c) = " << mPar->Tau << endl
       << "Deltat (fm/c) = " << fabs(mPar->Deltat) << endl;
}

void BlastWave::testForMike(const char* aFileName){
  ifstream tFMike(aFileName);
  char tDum[50];
  for(int ti=0; ti<57; ti++) tFMike >> tDum;
  cout << tDum << endl;
  double aPt;
  double aPhiP;
  double aN;
  double aX2;
  double aY2;
  double aXY;
  double aXT;
  double aYT;
  double aT2;
  double aZ2;
  double aX;
  double aY;
  double aT;
  double aRO2;
  double aRS2;
  double aROS2;
  double aRL2;	
  double aMass = 0.139;

  double tT = fabs(mPar->T);
		
  int tNBinX = 40;
  double tMinX = -1.*mPar->Rx;
  double tMaxX = mPar->Rx;
  double tStepX = (tMaxX-tMinX)/(tNBinX-1.);
  tMinX += tStepX/2.; 
  int tNBinY = 40;
  double tMinY = -1.*mPar->Ry;
  double tMaxY = mPar->Ry;
  double tStepY = (tMaxY-tMinY)/(tNBinY-1.);
  tMinY += tStepY/2.; 
	
  tFMike >> aPt;
  double tPrevPt = -1.;
  while(!tFMike.eof()){
    tFMike >> tDum >> aN >> aX2 >> aY2 >> aXY >> aXT >> aYT >> aT2 
	   >> aZ2 >> aX >> aY >> aT >> aRO2 >> aRS2 >> aROS2 >> aRL2;
    
    if(tPrevPt==aPt){
      aPhiP+= TMath::Pi()/8.;		
    }
    else{
      aPhiP = 0.;
    }        
    tPrevPt = aPt;
    
    double tMt = sqrt(aPt*aPt+aMass*aMass);
    double tSum = 0.;
    double tX = 0.;
    double tX2 = 0.;
    double tY = 0.;
    double tY2 = 0.;
    double tXY = 0.;
    double tZ2 = 0.;
    double tT = 0.;
    double tT2 = 0.;
    double tXT = 0.;
    double tYT = 0.;
    double tR;
    double tPhiB;
    double tRho; 
    double tBeta; 
    double tK0Beta; 
    double tK1Beta; 
    double tK2Beta; 
    double tBetaWeight; 
    double tAlphaWeight;
    double tH1OverH0 = (mPar->Deltat*mPar->Deltat+mPar->Tau*mPar->Tau)/mPar->Tau;
    double tH2OverH0 = (3*mPar->Deltat*mPar->Deltat+mPar->Tau*mPar->Tau);
    
    for(double tLX =tMinX; tLX <= tMaxX; tLX += tStepX){
      for(double tLY =tMinY; tLY <= tMaxY; tLY += tStepY){
	tR = sqrt(tLX*tLX/mPar->Rx/mPar->Rx+tLY*tLY/mPar->Ry/mPar->Ry);
	if(tR<1.){
	  tPhiB = TMath::ATan2(tLY*mPar->Rx*mPar->Rx/mPar->Ry/mPar->Ry,tLX);
	  tRho = tR*(mPar->Rho0+
		     mPar->Rho2*cos(2*tPhiB)+
		     mPar->Rho4*cos(4*tPhiB));	
	  tBeta = tMt/tT*TMath::CosH(tRho);
	  tK0Beta = TMath::BesselK0(tBeta);
	  tK1Beta = TMath::BesselK1(tBeta);
	  tK2Beta = 2./tBeta*tK1Beta+tK0Beta;
	  tBetaWeight = tK1Beta;
	  tAlphaWeight = exp(aPt/tT*TMath::SinH(tRho)*cos(aPhiP-tPhiB));	  
	  tSum += (tAlphaWeight*tBetaWeight);
	  tX   += (tAlphaWeight*tBetaWeight * tLX);
	  tX2  += (tAlphaWeight*tBetaWeight * tLX * tLX);
	  tY   += (tAlphaWeight*tBetaWeight * tLY);
	  tY2  += (tAlphaWeight*tBetaWeight * tLY * tLY);
	  tXY  += (tAlphaWeight*tBetaWeight * tLX * tLY);
	  tBetaWeight = tK1Beta/tBeta+tK0Beta;
	  tT   += (tAlphaWeight*tBetaWeight * tH1OverH0);
	  tXT  += (tAlphaWeight*tBetaWeight * tLX * tH1OverH0);
	  tYT  += (tAlphaWeight*tBetaWeight * tLY * tH1OverH0);
	  tBetaWeight = 1./tBeta*tK2Beta;
	  tZ2  += (tAlphaWeight*tBetaWeight * tH2OverH0); 
	  tBetaWeight += tK1Beta;
	  tT2  += (tAlphaWeight*tBetaWeight * tH2OverH0);
	  
	}
      }
    }
    mResult.dNdptpt = tSum*tMt;

    tX  /= tSum;
    tX2 /= tSum; tX2 -= (tX*tX);
    tY  /= tSum;
    tY2 /= tSum; tY2 -= (tY*tY);
    tXY /= tSum; tXY -= (tX*tY);
    tZ2 /= tSum;
    tT  /= tSum;
    tT2 /= tSum; tT2 -= (tT*tT);
    tXT /= tSum; tXT -= (tX*tT);
    tYT /= tSum; tYT -= (tY*tT);
    
    double tBetaT = aPt/tMt;
    double tGammaT = tMt/aMass;
    double tCosPhiP = cos(aPhiP);
    double tSinPhiP = sin(aPhiP);

    mResult .AvROutStarSpace = tGammaT * 
      (tX * cos(aPhiP) + tY * sin(aPhiP));
    mResult .AvROutStarTime = - tGammaT * tBetaT * tT;
    mResult .AvROutStar =  mResult .AvROutStarSpace + mResult .AvROutStarTime;
    mResult.ROutSquare = 0.5*(tX2+tY2)+0.5*(tX2-tY2)*cos(2*aPhiP)
      + tXY*sin(2*aPhiP)
      - 2*tBetaT*(tXT*cos(aPhiP)+tYT*sin(aPhiP))+tBetaT*tBetaT*tT2;
    mResult.AvRSide = tY * cos(aPhiP) - tX * sin(aPhiP);
    mResult.RSideSquare = 0.5*(tX2+tY2)-0.5*(tX2-tY2)*cos(2*aPhiP)
      -tXY*sin(2*aPhiP);
    mResult.ROutSideSquare = tXY*cos(2*aPhiP)-0.5*(tX2-tY2)*sin(2*aPhiP)
      + tBetaT*(tXT*sin(aPhiP)-tYT*cos(aPhiP));
    mResult.AvRLongStar = 0; // Not implemented yet
    mResult.RLongSquare = tZ2;
    
    cout << aPt << " " << aPhiP << " "
    << aY << " " << tY << endl;
           /*<<  aN << " " << tSum << " " //(aN-tSum)/(aN+tSum)*2. << " "
           << (aX2-tX2)/(aX2+tX2)*2. << " " 
           << (aY2-tY2)/(aY2+tY2)*2. << " " 
           << ((aXY<0.001 && tXY<0.001) ? 0 : (aXY-tXY)/(aXY+tXY)*2.)  << " "
           <<  ((aXT<0.001 && tXT<0.001) ? 0 : (aXT-tXT)/(aXT+tXT)*2.) << " " 
           << ((aYT<0.001 && tYT<0.001) ? 0 : (aYT-tYT)/(aYT+tYT)*2.) << " " 
           << (aT2-tT2)/(aT2+tT2)*2. <<  " "
           << (aZ2-tZ2)/(aZ2+tZ2)*2. <<  " " 
           << ((aX<0.001 && tX<0.001) ? 0 : (aX-tX)/(aX+tX)*2.) <<  " "
           << ((aY<0.001 && tY<0.001) ? 0 : (aY-tY)/(aY+tY)*2.) <<  " "
           << (aT-tT)/(aT+tT)*2. <<  " "
           << (aRO2-mResult.ROutSquare)/(aRO2+mResult.ROutSquare)*2. <<  " "
           << (aRS2-mResult.RSideSquare)/(aRS2+mResult.RSideSquare)*2. <<  " "
           << ((aROS2<0.001 && mResult.ROutSideSquare<0.001) ? 0 : (aROS2-mResult.ROutSideSquare)/(aROS2+mResult.ROutSideSquare)*2.) <<  " "
           << (aRL2-mResult.RLongSquare)/(aRL2+mResult.RLongSquare)*2. << " "
  			<< endl;*/
           	tFMike >> aPt;
  }  
}

TH1D* BlastWave::meanPtVsPhi(const char* aNamePref, double aMass){
  calcEStructHisto(aNamePref,aMass);
  return mHMeanPtVsPhi;
}
TH1D* BlastWave::autoCorrVsPhi(const char* aNamePref, double aMass){
  calcEStructHisto(aNamePref,aMass);
  return mHAutoCorrVsPhi;
}

double BlastWave::meanPt(double aPhip, double aMass,
			 double aMinPt, double aMaxPt){
  int tNBin=20;
  double tStep=(aMaxPt-aMinPt)/tNBin;
  double tPtSum =0.;
  double tSum=0.;
  for(double tPt=aMinPt+tStep/2.; tPt<aMaxPt; tPt+=tStep){
    calc(aMass,tPt,aPhip);
    tPtSum+=mResult.dNdptpt*tPt*tPt;
    tSum+=mResult.dNdptpt*tPt;
  }
  return tPtSum/tSum;
}

//6046652627 daryl

void BlastWave::calcEStructHisto(const char* aNamePref, double aMass){
  static double tMass;
  if(aMass!=tMass){
    mSpace=0;
    tMass=aMass;
    int tNPhiBin=24;
    char tHName[50];
    sprintf(tHName,"%sMeanPtVsPhi", aNamePref);
    mHMeanPtVsPhi = new TH1D(tHName,"",tNPhiBin,
			     -1./tNPhiBin*TMath::Pi(),
			     (1.-1./tNPhiBin/2.)*2*TMath::Pi());
    char tHName1[50];
    sprintf(tHName1,"%sAutoCorrVsPhi", aNamePref);
    mHAutoCorrVsPhi = new TH1D(tHName1,"",
			       tNPhiBin,
			       -0.5/(tNPhiBin-1.)*TMath::Pi(),
			       (1+0.5/(tNPhiBin-1.))*TMath::Pi());
    double tMeanPt;
    double tMeanPtSum=0.;
    double tPhiP;
    double tMinPt = 0.15;
    double tMaxPt = 2.;
    int tNPtBin=20;
    double tPtStep=(tMaxPt-tMinPt)/(tNPtBin-1.);
    double tPt;
    double tPtP2;
    double tdNdptpt1;
    double tPtWeight2;
    double tPtWeight1;
    double tPtWeight12;
    double tPtWeightRef;
    double tYield1;
    double tYield2;
    double tYield12;
    double tYieldRef;
    double tPtSum=0.;
    double tPt2Sum=0.;
    double tSum=0.;
    double tStep = (mHMeanPtVsPhi->GetBinWidth(1)*tPtStep);
    tStep *= tStep;
    for(int ti=1; ti<=tNPhiBin; ti++){
      tPhiP = mHMeanPtVsPhi->GetBinCenter(ti);
      tMeanPt = meanPt(tPhiP,aMass,0.15,2.);
      mHMeanPtVsPhi->SetBinContent(ti,tMeanPt);
      for(tPt=tMinPt; tPt<tMaxPt; tPt+=tPtStep){
	calc(aMass,tPt,tPhiP);
	tPt2Sum += mResult.dNdptpt*tPt*tPt*tPt;
	tPtSum += mResult.dNdptpt*tPt*tPt;
	tSum += mResult.dNdptpt*tPt;
      }
    }
    tPtSum/=tNPhiBin;
    tSum/=tNPhiBin;
    cout << tPtSum << " " << tSum << endl;
    //tPt2Sum/=tSum;
    double tPhiDelta;
    double tMeanPtProduct;
    double tAv=0.;
    for(int ti=1; ti<=tNPhiBin; ti++){
      tPhiDelta = mHAutoCorrVsPhi->GetBinCenter(ti);
      tMeanPtProduct=0.;

      tYield12=0.;
      tPtWeight12=0.;
      tYieldRef =0.;
      tPtWeightRef =0.;
      for(int tj=1; tj<=tNPhiBin; tj++){
	tPhiP = mHMeanPtVsPhi->GetBinCenter(tj);
	tPtWeight1 =0.;
	tPtWeight2 =0.;
	tYield1=0.;
	tYield2=0.;
	for(tPt=tMinPt; tPt<tMaxPt; tPt+=tPtStep){
	  tPtP2 = tPt*tPt;
	  calc(aMass,tPt,tPhiP);
	  tYield1 += (mResult.dNdptpt*tPt);
	  tPtWeight1 += (mResult.dNdptpt*tPtP2);
	  calc(aMass,tPt,tPhiP+tPhiDelta);
	  tYield2 += (mResult.dNdptpt*tPt);
	  tPtWeight2 += (mResult.dNdptpt*tPtP2);
	}
	tPtWeight12 += tPtWeight1*tPtWeight2;
	tYield12 += tYield1*tYield2;
	tPtWeightRef += tPtWeight1*tPtSum;
	tYieldRef+= tYield1*tSum;
      }
      tPtWeightRef *= tStep;
      tYieldRef *= tStep;
      tPtWeight12 *= tStep;
      cout << tPhiDelta << " " 
	   << tPtWeight12 << " " 
	   << tPtWeightRef << " " 
	   << tYieldRef << " "
	   << tPtWeight12-tPtWeightRef << " "
	   << (tPtWeight12-tPtWeightRef)/sqrt(tYieldRef) << " "
	   << endl;
	//<< tYieldRef << " " << tPt2Sum << " " 
	//<< tPt2Sum-tPtSum*tPtSum << " " 
	//   << tPtSum*tPtSum << endl;
      //tAv += tPtWeight12;
      mHAutoCorrVsPhi->SetBinContent(ti,(tPtWeight12-tPtWeightRef)/sqrt(tYieldRef));
    }
    /*    tAv/=tNPhiBin;
    for(int ti=1; ti<=tNPhiBin; ti++){
      mHAutoCorrVsPhi->SetBinContent(ti,(mHAutoCorrVsPhi->GetBinContent(ti)-tAv)/sqrt(tAv));
    }
    */
  }
}
