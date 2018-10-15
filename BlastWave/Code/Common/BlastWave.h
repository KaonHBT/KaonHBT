#ifndef BlastWave_h
#define BlastWave_h
#include "TMath.h"
class TH1D;
class TH2D;

#include <iomanip>
#include <cmath>
using namespace std;


struct BlastWaveResult{
  double ROutSquare;
  double RSideSquare;
  double ROutSideSquare;
  double RLongSquare;
  double ROutSquareOsc;
  double RSideSquareOsc;
  double ROutSideSquareOsc;
  double RLongSquareOsc;
  double AvROutStarSpace;
  double AvROutStarTime;
  double AvROutStar;
  double AvRSide;
  double AvRLongStar;
  double AvROutStarOsc;
  double AvRSideOsc;
  double AvRLongStarOsc;
  double V2;
  double V4;
  double V6;
  double V8;
  double dNdptpt;
};

struct BWParameters{
  double T;
  double Rho0;
  double Rho2;
  double Rho4;
  double Rx;
  double Ry;
  double As;
  double Tau; 
  double Deltat;
};

class BlastWave{
 public:
  BlastWave(double aT,
	    double aRho0,
	    double aR, double aAs, 
	    double aTau, double aDeltat);
  BlastWave(double aT,
	    double aRho0, double aRho2, double aRho4,
	    double aRx, double aRy, double aAs,

	    double aTau, double aDeltat);
  //~BlastWave(){}
  void setParameters(double aT,
		     double aRho0,
		     double aR, double aAs,
		     double aTau, double aDeltat);
  void setParameters(double aT,
		     double aRho0, double aRho2, double aRho4,
		     double aRx, double aRy, double aAs,
		     double aTau, double aDeltat);
  //void set(TMinuit* aMinuit, double aFixParameters);
  void setBoseEinsteinStatistics();
  void setFermiDiracStatistics();
  void setBoltzmanStatistics();
  void setStat(int aStat);

  // ______________________________________________________________
  // --- Internal 
  void setParameters(double* aPar); // to be used only by BlastWaveFitter
  void useR();
  void doNotUseR();
  void setSpace(int aSpace){mSpace=aSpace;}
  const BlastWaveResult& getBlastWaveResult(double aMass,double aPt,
					    double aPhiP);
  const BlastWaveResult& getBlastWaveResult(double aMass,double aPt);
  const BlastWaveResult& getBlastWaveNoSpaceResult(double aMass,double aPt);

  double getEmissionProba(double aMass, double aX, double aY,
			  double aPt, double aPhiP, double aEta);
  double getEmissionProba(double aMass, double aX, double aY,
			  double aPt, double aPhiP);

  // _______________________________________________________________
  // --- Output
  TH2D* emMomDist(double aMass, double aX, double aY);
  TH2D* emPosDist(double aMass, double aPt, double aPhiP);
  TH2D* emMomDist(double aMass);
  TH2D* emPosDist(double aMass);
//  TH2D* TimeVsPtDist(double aMass);  
  TH1D* ZDist(double aMass, double aPt);  

  TH1D* spectra(const char* aNamePref, double aMass);
  TH1D* v2VsPt(const char* aNamePref, double aMass);
  TH1D* v4VsPt(const char* aNamePref, double aMass);
  TH1D* v6VsPt(const char* aNamePref, double aMass);
  TH1D* v8VsPt(const char* aNamePref, double aMass);

  TH1D* rOut(const char* aNamePref, double aMass);
  TH1D* rSide(const char* aNamePref, double aMass);
  TH1D* rLong(const char* aNamePref, double aMass);
  TH1D* rOutVsMt(const char* aNamePref, double aMass);
  TH1D* rSideVsMt(const char* aNamePref, double aMass);
  TH1D* rLongVsMt(const char* aNamePref, double aMass);

  TH1D* rOutSquare(const char* aNamePref, double aMass);
  TH1D* rSideSquare(const char* aNamePref, double aMass);
  TH1D* rLongSquare(const char* aNamePref, double aMass);
  TH1D* rOutSideSquare(const char* aNamePref, double aMass);

  TH1D* rOutSquareOsc(const char* aNamePref, double aMass);
  TH1D* rSideSquareOsc(const char* aNamePref, double aMass);
  TH1D* rLongSquareOsc(const char* aNamePref, double aMass);
  TH1D* rOutSideSquareOsc(const char* aNamePref, double aMass);

  TH1D* avOut(const char* aNamePref, double aMass);  
  TH1D* avTime(const char* aNamePref, double aMass);  
  TH1D* avSide(const char* aNamePref, double aMass);      
  
  TH1D* avSepOut(const char* aNamePref, double aMass1, double aMass2);
  TH1D* avSepOutSpace(const char* aNamePref, double aMass1, double aMass2);
  TH1D* avSepOutTime(const char* aNamePref, double aMass1, double aMass2);
  TH1D* avSepSide(const char* aNamePref, double aMass1, double aMass2);
  TH1D* avSepLong(const char* aNamePref, double aMass1, double aMass2);

  TH1D* avSepOutVsPhiP(const char* aNamePref, double aMass1, double aMass2, double aPt2);
  TH1D* avSepOutSpaceVsPhiP(const char* aNamePref, double aMass1, double aMass2, double aPt2);
  TH1D* avSepOutTimeVsPhiP(const char* aNamePref, double aMass1, double aMass2, double aPt2);
  TH1D* avSepSideVsPhiP(const char* aNamePref, double aMass1, double aMass2, double aPt2);
  TH1D* avSepLongVsPhiP(const char* aNamePref, double aMass1, double aMass2, double aPt2);  
  
  
  TH1D* rOutSquareVsPhiP(const char* aNamePref, double aMass, double aPt);
  TH1D* rSideSquareVsPhiP(const char* aNamePref, double aMass, double aPt);
  TH1D* rLongSquareVsPhiP(const char* aNamePref, double aMass, double aPt);
  TH1D* rOutSideSquareVsPhiP(const char* aNamePref, double aMass, double aPt);

  TH1D* meanPtVsPhi(const char* aNamePref, double aMass);
  TH1D* autoCorrVsPhi(const char* aNamePref, double aMass);
  

  BWParameters* par();

  void coutParameters();

  void testForMike(const char* aFileName)  ;
  
 private:
  BWParameters* mPar;

  int mSpace; // turn on spatial calculations
  int mNewPar;
  int mStat;

  int mUseR;

  double mPrevMass;
  double mPrevPt;
  double mPrevPhiP;
  BlastWaveResult mResult;
  void calc(double aMass,double aPt,double aPhiP);
  void calc(double aMass,double aPt);
  // different way of doing the integral
  void calcR(double aMass,double aPt,double aPhiP); 
  void calcRTilda(double aMass,double aPt,double aPhiP);
  void calcXY(double aMass,double aPt,double aPhiP); 

  void calcHisto(const char* aNamePref, double aMass);
  int mNewHisto;
  TH1D* mHSpectra;
  TH1D* mHV2;
  TH1D* mHV4;
  TH1D* mHV6;
  TH1D* mHV8;

  TH1D* mHROut;
  TH1D* mHRSide;
  TH1D* mHRLong;
  TH1D* mHROutVsMt;
  TH1D* mHRSideVsMt;
  TH1D* mHRLongVsMt;
  TH1D* mHROutSquare;
  TH1D* mHRSideSquare;
  TH1D* mHRLongSquare;
  TH1D* mHROutSideSquare;
  TH1D* mHROutSquareOsc;
  TH1D* mHRSideSquareOsc;
  TH1D* mHRLongSquareOsc;
  TH1D* mHROutSideSquareOsc;

  TH1D* mHAvOut;
  TH1D* mHAvTime;
  TH1D* mHAvSide;    
  
  void calcAvSepHisto(const char* aNamePref, double aMass1, double aMass2);
  int mNewAvSepHisto; 
  TH1D* mHAvSepOut;
  TH1D* mHAvSepOutTime;
  TH1D* mHAvSepOutSpace;
  TH1D* mHAvSepSide;
  TH1D* mHAvSepLong;

  void calcAvSepVsPhiPHisto(const char* aNamePref, double aMass1, double aMass2, double aPt2);
  int mNewAvSepVsPhiPHisto; 
  TH1D* mHAvSepOutVsPhiP;
  TH1D* mHAvSepOutTimeVsPhiP;
  TH1D* mHAvSepOutSpaceVsPhiP;
  TH1D* mHAvSepSideVsPhiP;
  TH1D* mHAvSepLongVsPhiP;  
  
  
  void calcPhiPDependentHisto(const char* aNamePref, double aMass, double aPt);
  int mNewPhiPDepHisto;
  TH1D* mHROutSquareVsPhiP;
  TH1D* mHRSideSquareVsPhiP;
  TH1D* mHRLongSquareVsPhiP;
  TH1D* mHROutSideSquareVsPhiP;

  double meanPt(double aPhip, double aMass,
		double aMinPt, double aMaxPt);
  void calcEStructHisto(const char* aNamePref, double aMass);
  TH1D* mHMeanPtVsPhi;
  TH1D* mHAutoCorrVsPhi;
};

inline void BlastWave::calc(double aMass, double aPt, double aPhiP){
  //calcRTilda(aMass,aPt,aPhiP);
  calcR(aMass,aPt,aPhiP);
}

inline void BlastWave::setParameters(double* aPar){
  mNewHisto=1;
  mPar = (BWParameters*) aPar;
  mPar->As = fabs(mPar->As);
  if(mUseR)  mPar->Ry = mPar->Rx;
}

inline void BlastWave::useR(){
  mUseR=1;
}
inline void BlastWave::doNotUseR(){
  mUseR=0;
}
inline const BlastWaveResult& BlastWave::getBlastWaveResult(double aMass,
							   double aPt,
							   double aPhiP){
  if(mNewPar || mPrevMass!=aMass || mPrevPt!=aPt || mPrevPhiP!=aPhiP){
    mNewPar=0; mPrevMass=aMass; mPrevPt=aPt; mPrevPhiP=aPhiP;
    mSpace = 1;
    calc(aMass, aPt, aPhiP);
  }
  return mResult;
}
inline const BlastWaveResult& BlastWave::getBlastWaveResult(double aMass,
							   double aPt){
  if(mNewPar || mPrevMass!=aMass || mPrevPt!=aPt){
    mNewPar=0; mPrevMass=aMass; mPrevPt=aPt;
    mSpace = 1;
    calc(aMass, aPt);
  }
  return mResult;
}
inline const BlastWaveResult& BlastWave::getBlastWaveNoSpaceResult(double aMass, double aPt){
  if(mNewPar || mPrevMass!=aMass || mPrevPt!=aPt){
    mNewPar=0; mPrevMass=aMass; mPrevPt=aPt;
    mSpace = 0;
    calc(aMass, aPt);
  }
  return mResult;
}

inline BWParameters* BlastWave::par(){return mPar;}

inline double BlastWave::getEmissionProba(double aMass, double aX, double aY,
					  double aPt, double aPhiP, 
					  double aEta){
  double tR = sqrt(aX*aX/mPar->Rx/mPar->Rx+aY*aY/mPar->Ry/mPar->Ry);
  double tSpaceWeight = mPar->As==0.? (tR<1.) : 1./(1.+exp((tR-1.)/mPar->As));
  //cout << mPar->As << " " << tR << " " <<  tSpaceWeight << endl;
  if(tSpaceWeight>0.){
    double tPhiB = TMath::ATan2(aY*mPar->Rx*mPar->Rx/mPar->Ry/mPar->Ry,aX);
    double tRho = tR*(mPar->Rho0+mPar->Rho2*cos(2*tPhiB));
    double tCoshEta = TMath::CosH(aEta);
    double tMt = sqrt(aPt*aPt+aMass*aMass);
    return tMt*exp(-tMt/mPar->T*TMath::CosH(tRho)*tCoshEta)*
      tCoshEta*exp(aPt/mPar->T*TMath::SinH(tRho)*cos(aPhiP-tPhiB))*
      tSpaceWeight;    
  }
  else{
    return 0.;
  }
}

inline double BlastWave::getEmissionProba(double aMass, double aX, double aY,
					  double aPt, double aPhiP){
  double tR = sqrt(aX*aX/mPar->Rx/mPar->Rx+aY*aY/mPar->Ry/mPar->Ry);
  double tSpaceWeight = mPar->As==0.? (tR<1.) : 1./(1.+exp((tR-1.)/mPar->As));
  //cout << mPar->As << " " << tR << " " <<  tSpaceWeight << endl;
  if(tSpaceWeight>0.){ 
    double tPhiB = TMath::ATan2(aY*mPar->Rx*mPar->Rx/mPar->Ry/mPar->Ry,aX);
    double tRho = tR*(mPar->Rho0+mPar->Rho2*cos(2*tPhiB));
    double tMt = sqrt(aPt*aPt+aMass*aMass);
    return tMt*TMath::BesselK1(tMt/mPar->T*TMath::CosH(tRho))
      *exp(aPt/mPar->T*TMath::SinH(tRho)*cos(aPhiP-tPhiB))*
      tSpaceWeight;        
  }
  else{
    return 0.;
  }
}
inline void BlastWave::setBoseEinsteinStatistics(){mStat=1;}
inline void BlastWave::setFermiDiracStatistics(){mStat=-1;}
inline void BlastWave::setBoltzmanStatistics(){mStat=0;}
inline void BlastWave::setStat(int aStat){
  mStat= aStat;
}


#endif
