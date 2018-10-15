#ifndef BlastWaveFitter_h
#define BlastWaveFitter_h

class TGraph;
class TGraphErrors;
class TMinuit;
class TH2D;
class TH1;
#include "BlastWave.h"
#include "DataForFit.h"

#include <iostream>
using namespace std;

class BlastWaveFitter{
 public:
  // __________________________________________________________________
  // --- Init
  static BlastWaveFitter* instance();
  ~BlastWaveFitter();

  // __________________________________________________________________
  // --- Data input functions
  // X = Pt | Y prorporotional to dN/dpt/pt 
  SpectraForFit* addSpectra(TGraphErrors* aSpectra, 
			    double aMass, 
			    double aPtMin=-9999.,
			    double aPtMax=9999., int aStat=0); 
  // X = Pt | Y prorporotional to dN/dpt/pt 
  SpectraForFit* addSpectra(TH1* aHisto, 
			    double aMass, 
			    double aPtMin=-9999., 
			    double aPtMax=9999., int aStat=0); 
  // X = pt | Y = radii
  void addHbtRadii(TGraphErrors* aROut, 
		   TGraphErrors* aRSide, 
		   TGraphErrors* aRLong,
		   double aMass, int aStat=0); 
  // X = pt | Y = radii square | phip = angle wrt reaction plane
  void addHbtRadiiSquareWRTReactionPlane(TGraphErrors* aROutSquare, 
			      TGraphErrors* aRSideSquare, 
			      TGraphErrors* aROutSideSquare, 
			      TGraphErrors* aRLongSquare, 
			      double aMass,
			      double aPhiP, int aStat=0); 
  // X = pt | Y = radii square | phip = 45 degree
  void addHbtRadiiSquare(TGraphErrors* aROutSquare, 
			 TGraphErrors* aRSideSquare, 
			 TGraphErrors* aROutSideSquare,
			 TGraphErrors* aRLongSquare,  
			 double aMass, int aStat=0); 			      
  // X = pt | Y = second order Fourier coef of the HBT radii
  void addHbtRadiiSquareOsc(TGraphErrors* aROutSquareOsc,
			    TGraphErrors* aRSideSquareOsc, 
			    TGraphErrors* aROutSideSquareOsc,
			    TGraphErrors* aRLongSquareOsc, 
			    double aMass, int aStat=0);
  // X = pt | Y = v2 in absolute value (No %)
  void addV2(TGraphErrors* aV2, 
	     double aMass, 
	     double aPtMin=-9999., 
	     double aPtMax=9999., 
	     int aStat=0); 
  // X = pt | Y = v2 in absolute value (No %)
  void addV4(TGraphErrors* aV2, 
	     double aMass, 
	     double aPtMin=-9999., 
	     double aPtMax=9999., 
	     int aStat=0); 

  // __________________________________________________________________
  // --- Set Blast Wave parameters
  // Use this to fit spectra
  void setParameters(double aT,
		     double aRho0);
  // Use this to fit v2
  void setParameters(double aT,
		     double aRho0,
		     double aRho2,
		     double aRho4,
		     double aAspectRatio); // ROut-of-plane/RIn-plane
  // Use this to fit HBT intgrated in phi 
  void setParameters(double aT,
		     double aRho0,
		     double aR, double aAs,
		     double aTau, double aDeltat);
  // Use this to fit HBT wrt to reaction plane
  void setParameters(double aT,
		     double aRho0, 
		     double aRho2,
		     double aRho4,
		     double aRx, // In-plane 
		     double aRy, // Out-of-plane
		     double aAs,
		     double aTau, 
		     double aDeltat);
  void fixT();
  void fixRho0();
  void fixR(); // Fix Rx and Ry
  void useR(); // Equalize Rx and Ry   
  void fixTau();
  void fixDeltat();
  void fixRho2();
  void fixRho4();
  void fixRx();
  void fixRy();
  void fixAs();
  void freeAllParameters();
  void fixParameter(int aIndex);

  // __________________________________________________________________
  // --- Perform fit
  void minimize();
  void calculatePreciseErrors(); // call minos

  // __________________________________________________________________
  // --- Output
  TGraph* betaTContour(int aNSigma);
  TGraph* S2VsRho2Contour(int aNSigma);
  TGraph* TauVsDeltatContour(int aNSigma);
  // i1 and i2 are the parameters indexes
  TGraph* contour(char* aName, int aNSigma, int i1, int i2);
  TH2D* confLevelMap(char* aName, 
  	int ai1, int aNBinX, double aXMin, double aXMax, 
  	int ai2, int aNBinY, double aYMin, double aYMax);
  TH2D* confLevelMap(char* aName, int ai1, int aNBinX,
		     int ai2, int aNBinY);
  double chi2();
  int nDOF();
 
  void getParameter(Int_t parNo, 
		    Double_t& currentValue, 
		    Double_t& currentError);
  void coutFitResults(int ScaleErrorByChi2PerDof=0);
  BlastWave* blastWave(); // May be used to make access blast wave output histoes

  // __________________________________________________________________
  // --- Internal functions. Do not use.
  void clearData();
  int parameterIsFixed(int aIndex);

  // __________________________________________________________________
  // --- Private
 private: 
  BlastWaveFitter();
  static BlastWaveFitter* mInstance;
  
  int mDataArraySize;
  void adjustDataSize();
  int mNData;
  DataForFit** mData;
  int mSpace;

  int mNParameter;
  int mNParameterMax;
  int mNSpectra; // Add 1 parameters (scale) per spectra
  int mParameterIsFixed[10];
  BlastWave* mBlastWave;

  TMinuit* mMinuit;
  void setMinuitParameters();
  int* mContinueCalcChi2;
};

// ___________________________________________________________________
// ___________________________________________________________________
// --- Inline functions

inline int BlastWaveFitter::parameterIsFixed(int aIndex){
  return mParameterIsFixed[aIndex];
}

inline BlastWave* BlastWaveFitter::blastWave(){
  return mBlastWave;
}
inline double BlastWaveFitter::chi2(){
  double tChi2Sum(0.);
  int tContinue(1);
  double tChi2;
  while(tContinue){
    tContinue=0;
    for(int ti=0; ti<mNData; ti++){
      if(mContinueCalcChi2[ti]) {
	tChi2 = mData[ti]->chi2With(mBlastWave);
	if(tChi2!=-1.){
	  tChi2Sum += tChi2;
	  tContinue=1;
	}
	else{
	  mContinueCalcChi2[ti]=0;
	}
      }
      //cout << ti << " " << tChi2 << " " << tChi2Sum << " ";
    }
    //cout << endl;
  }
  for(int ti=0; ti<mNData; ti++){mContinueCalcChi2[ti]=1;}
  return tChi2Sum;
}
inline int BlastWaveFitter::nDOF(){
  int tNDOF=0;
  for(int ti=0; ti<mNData; ti++){
    tNDOF += mData[ti]->nDOF();
  }
  return tNDOF-mNParameter-mNSpectra;
}








#endif

