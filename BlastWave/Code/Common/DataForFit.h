#ifndef DataForFit_h
#define DataForFit_h
#include "BlastWave.h"

class TGraphErrors;
class TH1;

#include <iomanip>
#include <cmath>
using namespace std;

class DataForFit{
 public:
  DataForFit(TGraphErrors* aGraph, double aMass, 
	     double aPtMin=-9999., double aPtMax=9999., int aStat=0);
  DataForFit(TH1* aGraph, double aMass, 
	     double aPtMin=-9999., double aPtMax=9999., int aStat=0);
  void setBoseEinsteinStat(){mStat=1;}
  void setFermiDiracStat()  {mStat=-1;}


  virtual ~DataForFit();
  void addDataForFit(TGraphErrors* aGraph);

  double chi2With(BlastWave* aBW, int aAllPoint=0);
 
  int nDOF();
  virtual TH1D* histo(const char* aName, int NBin, BlastWave* aBW);
  virtual TH1D* histo(const char* aName, int NBin, double XMin, double XMax,
		      BlastWave* aBW);
  virtual const char* getName();
            
 protected:
  int mNPoints;
  int mCurrentPoint;
  double* mX;
  double* mY;
  double* mYErrSquare;
  char mName[50];

  double mMass;
  virtual double getBlastWaveY(BlastWave* aBW, double aPt) =0;

  int mStart;
  int mSpectra;

  int mStat;
};

inline int DataForFit::nDOF(){
  return mNPoints;
}
inline double DataForFit::chi2With(BlastWave* aBW, int aAllPoint){
  if(mCurrentPoint==mNPoints || aAllPoint){
    mCurrentPoint=0;
    mStart = 1;
    if(!aAllPoint) return -1.;
  }
  aBW->setStat(mStat);
  double tChi2(0.);
  double tVal;
  int tEndPoint = (mSpectra+aAllPoint)? mNPoints : mCurrentPoint+1;
  for(int ti=mCurrentPoint; ti<tEndPoint; ti++){
    tVal = mY[ti] - getBlastWaveY(aBW, mX[ti]);
    tVal *= tVal;
    tChi2 += (tVal/mYErrSquare[ti]);
  }
  mCurrentPoint= aAllPoint? 0 : tEndPoint;
  return tChi2;
}

inline const char* DataForFit::getName() { return mName;}

// ____________________________________________________________________

// ______________________________________________________________
// ------ HBT radii
// ------
class ROutSquareForFit : public DataForFit{
 public:
  ROutSquareForFit(TGraphErrors* aGraph, double aMass, double aPhiP, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat), mPhiAverage(0),mPhiP(aPhiP){}
  ROutSquareForFit(TGraphErrors* aGraph, double aMass)
    :DataForFit(aGraph, aMass), mPhiAverage(1){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  int mPhiAverage;
  double mPhiP;
};
inline double ROutSquareForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return (mPhiAverage? 
  		aBW->getBlastWaveResult(mMass, aPt).ROutSquare  :
  		aBW->getBlastWaveResult(mMass, aPt, mPhiP).ROutSquare) ;
}
// ------
class RSideSquareForFit : public DataForFit{
 public:
  RSideSquareForFit(TGraphErrors* aGraph, double aMass, double aPhiP, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat), mPhiAverage(0),mPhiP(aPhiP){}
  RSideSquareForFit(TGraphErrors* aGraph, double aMass)
    :DataForFit(aGraph, aMass), mPhiAverage(1){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  int mPhiAverage;
  double mPhiP;
};
inline double RSideSquareForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return (mPhiAverage? 
  		aBW->getBlastWaveResult(mMass, aPt).RSideSquare  :
  		aBW->getBlastWaveResult(mMass, aPt, mPhiP).RSideSquare) ;
}
// ------
class RLongSquareForFit : public DataForFit{
 public:
  RLongSquareForFit(TGraphErrors* aGraph, double aMass, double aPhiP, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat), mPhiAverage(0),mPhiP(aPhiP){}
  RLongSquareForFit(TGraphErrors* aGraph, double aMass)
    :DataForFit(aGraph, aMass), mPhiAverage(1){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  int mPhiAverage;
  double mPhiP;
};
inline double RLongSquareForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return (mPhiAverage? 
  		aBW->getBlastWaveResult(mMass, aPt).RLongSquare  :
  		aBW->getBlastWaveResult(mMass, aPt, mPhiP).RLongSquare) ;
}
// ------
class ROutSideSquareForFit : public DataForFit{
 public:
  ROutSideSquareForFit(TGraphErrors* aGraph, double aMass, double aPhiP, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat), mPhiAverage(0),mPhiP(aPhiP){}
  ROutSideSquareForFit(TGraphErrors* aGraph, double aMass)
    :DataForFit(aGraph, aMass), mPhiAverage(1){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  int mPhiAverage;
  double mPhiP;
};
inline double ROutSideSquareForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return (mPhiAverage? 
  		aBW->getBlastWaveResult(mMass, aPt).ROutSideSquare  :
  		aBW->getBlastWaveResult(mMass, aPt, mPhiP).ROutSideSquare) ;
}
// ------
class ROutSquareOscForFit : public DataForFit{
 public:
  ROutSquareOscForFit(TGraphErrors* aGraph, double aMass, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double ROutSquareOscForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).ROutSquareOsc;
}
// ------
class RSideSquareOscForFit : public DataForFit{
 public:
  RSideSquareOscForFit(TGraphErrors* aGraph, double aMass, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double RSideSquareOscForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).RSideSquareOsc;
}
// ------
class RLongSquareOscForFit : public DataForFit{
 public:
  RLongSquareOscForFit(TGraphErrors* aGraph, double aMass, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double RLongSquareOscForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).RLongSquareOsc;
}
// ------
class ROutSideSquareOscForFit : public DataForFit{
 public:
  ROutSideSquareOscForFit(TGraphErrors* aGraph, double aMass, int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat){}    
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double ROutSideSquareOscForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).ROutSideSquareOsc;
}
// ------
class ROutForFit : public DataForFit{
 public:
  ROutForFit(TGraphErrors* aGraph, double aMass, double aPhiP=0., int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat),mPhiP(aPhiP){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
};
inline double ROutForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return sqrt(aBW->getBlastWaveResult(mMass, aPt).ROutSquare);
}
// ------
class RSideForFit : public DataForFit{
 public:
  RSideForFit(TGraphErrors* aGraph, double aMass, double aPhiP=0., int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat),mPhiP(aPhiP){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
};
inline double RSideForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return  sqrt(aBW->getBlastWaveResult(mMass, aPt).RSideSquare);
}
// ------
class RLongForFit : public DataForFit{
 public:
  RLongForFit(TGraphErrors* aGraph, double aMass, double aPhiP=0., int aStat=0)
    :DataForFit(aGraph, aMass, -9999., 9999., aStat),mPhiP(aPhiP){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
};
inline double  RLongForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return  sqrt(aBW->getBlastWaveResult(mMass, aPt).RLongSquare);
}

// ______________________________________________________________
// ------ Non ID
// ------
class AvROutStarForFit : public DataForFit{
 public:
  AvROutStarForFit(TGraphErrors* aGraph, double aMass1, double aMass2, 
		   double aPhiP=0.)
    :DataForFit(aGraph,aMass1),mPhiP(aPhiP),mMass2(aMass2){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
  double mMass2;
};
inline double  AvROutStarForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).AvROutStar - 
    aBW->getBlastWaveResult(mMass2, aPt).AvROutStar;
}
// ------
class AvRSideForFit : public DataForFit{
 public:
  AvRSideForFit(TGraphErrors* aGraph, double aMass1, double aMass2, 
		   double aPhiP=0.)
    :DataForFit(aGraph,aMass1),mPhiP(aPhiP),mMass2(aMass2){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
  double mMass2;
};
inline double  AvRSideForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).AvRSide - 
    aBW->getBlastWaveResult(mMass2, aPt).AvRSide;
}
// ------
class AvRLongStarForFit : public DataForFit{
 public:
  AvRLongStarForFit(TGraphErrors* aGraph, double aMass1, double aMass2, 
		   double aPhiP=0.)
    :DataForFit(aGraph,aMass1),mPhiP(aPhiP),mMass2(aMass2){}
  double getBlastWaveY(BlastWave* aBW, double aPt);
 private:
  double mPhiP;
  double mMass2;
};
inline double  AvRLongStarForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveResult(mMass, aPt).AvRLongStar - 
    aBW->getBlastWaveResult(mMass2, aPt).AvRLongStar;
}


// ______________________________________________________________
// ------ Spectra
// ------
class SpectraForFit : public DataForFit{
 public:
  SpectraForFit(TGraphErrors* aGraph, double aMass, 
		double aPtMin=-9999., double aPtMax=9999., int aStat=0)
    : DataForFit(aGraph,aMass, aPtMin, aPtMax, aStat){mSpectra=1;};
  SpectraForFit(TH1* aHisto, double aMass, 
		double aPtMin=-9999., double aPtMax=9999., int aStat=0);
  double normBW();
  //double normData();
  //TH1D* histo(const char* aName, int NBin, BlastWave* aBW);

 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);

 private:
  double mNormBW;
  //double mNormData;
};
inline double SpectraForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  if(mStart){
    mStart = 0;
    double tBW;
    mNormBW =0.;
    double tNorm=0.;
    for(int ti=0; ti< mNPoints; ti++){
      tBW = aBW->getBlastWaveNoSpaceResult(mMass, mX[ti]).dNdptpt;
      tNorm += (tBW*tBW/mYErrSquare[ti]);
      mNormBW += (tBW*mY[ti]/mYErrSquare[ti]);
    }
    mNormBW/=tNorm;
  }
  return aBW->getBlastWaveNoSpaceResult(mMass, aPt).dNdptpt*mNormBW;
}
inline double SpectraForFit::normBW(){ return mNormBW;}
//inline double SpectraForFit::normData(){ return mNormData;}

// ______________________________________________________________
// ------ V2
// ------
class V2ForFit : public DataForFit{
 public:
  V2ForFit(TGraphErrors* aGraph, double aMass, 
	   double aMin=-9999., double aMax=9999., int aStat=0)
    :DataForFit(aGraph,aMass, aMin, aMax, aStat){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double V2ForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveNoSpaceResult(mMass, aPt).V2;
}
// ______________________________________________________________
// ------ V4
// ------
class V4ForFit : public DataForFit{
 public:
  V4ForFit(TGraphErrors* aGraph, double aMass, 
	   double aMin=-9999., double aMax=9999., int aStat=0)
    :DataForFit(aGraph,aMass, aMin, aMax, aStat){}
 protected:
  double getBlastWaveY(BlastWave* aBW, double aPt);
};
inline double V4ForFit::getBlastWaveY(BlastWave* aBW, double aPt){
  return aBW->getBlastWaveNoSpaceResult(mMass, aPt).V4;
}

#endif
