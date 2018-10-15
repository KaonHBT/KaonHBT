#include "BlastWaveFitter.h"
#include <TMinuit.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>

#include <iostream>
#include <cstdlib>
using namespace std;

BlastWaveFitter* BlastWaveFitter::mInstance =0;

BlastWaveFitter* gBlastWaveFitter;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  gBlastWaveFitter->blastWave()->setParameters(par);
  f = gBlastWaveFitter->chi2();
  //if(isnan(f)) {
  //cout << par[0] << " " << par[1] << " " << par[2] << " "
  // << par[3] << " " << par[4] << " " << par[5] << " " 
  // << par[6] << " " << par[7] << " " << par[8] << endl;
  //exit(0);
  //}
  //cout << f << " " << par[2] << " " << par[4] << " " << par[5] << endl;
}

BlastWaveFitter::BlastWaveFitter(){
  mBlastWave = new BlastWave(0.1,0.9,13.,0.,10.,2.);
  mNParameterMax=9;
  mNParameter = mNParameterMax; 
  mNSpectra = 0;
  mNData = 0;
  mDataArraySize = 10;
  mData = new DataForFit*[mDataArraySize];
  mContinueCalcChi2 = new int[mDataArraySize];
  for(int ti=0; ti<mDataArraySize; ti++)  mContinueCalcChi2[ti]=1;
  mMinuit = new TMinuit(mNParameter);
  mMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 0;
  mMinuit->mnexcm("SET STR", arglist ,1,ierflg);
  arglist[0] = 1.;
  mMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  setMinuitParameters();
  mSpace =0;
}

void BlastWaveFitter::setMinuitParameters(){
  Int_t ierflg = 0;
  BWParameters* tPar = mBlastWave->par();
  /*
  mMinuit->mnparm(0, "T", tPar->T, 0.01, 0.102274,0.1022741,ierflg);
  mMinuit->mnparm(1, "Rho0", tPar->Rho0, 0.005, 1.00313, 1.003131,ierflg);
  mMinuit->mnparm(2, "Rho2", tPar->Rho2, 0.001, -1., 1.,ierflg);
  mMinuit->mnparm(3, "Rho4", tPar->Rho4, 0.001, -1., 1.,ierflg);
  mMinuit->mnparm(4, "Rx", tPar->Rx, 0.1, 10.4827, 10.48271,ierflg);
  mMinuit->mnparm(5, "Ry", tPar->Ry, 0.1, 0., 100.,ierflg);
  mMinuit->mnparm(6, "As", tPar->As, 0.0005, 0., 0.5,ierflg);
  mMinuit->mnparm(7, "Tau", tPar->Tau, 0.1, 7.90111, 7.901111,ierflg);
  mMinuit->mnparm(8, "Deltat",tPar->Deltat, 0.01, 2.34736, 2.347361,ierflg);
  */
  ///*
  mMinuit->mnparm(0, "T", tPar->T, 0.03, 0.03, 0.4,ierflg);
  mMinuit->mnparm(1, "Rho0", tPar->Rho0, 0.005, 0., 3.,ierflg);
  mMinuit->mnparm(2, "Rho2", tPar->Rho2, 0.001, -1., 1.,ierflg);
  mMinuit->mnparm(3, "Rho4", tPar->Rho4, 0.001, -1., 1.,ierflg);
  mMinuit->mnparm(4, "Rx", tPar->Rx, 0.1, 0., 100.,ierflg);
  mMinuit->mnparm(5, "Ry", tPar->Ry, 0.1, 0., 100.,ierflg);
  mMinuit->mnparm(6, "As", tPar->As, 0.0005, 0., 0.5,ierflg);
  mMinuit->mnparm(7, "Tau", tPar->Tau, 0.1, 0., 100.,ierflg);
  mMinuit->mnparm(8, "Deltat",tPar->Deltat, 0.01, -20., 20.,ierflg);
  //*/
  for(int ti=0; ti<mNParameterMax; ti++){
    mParameterIsFixed[ti]=0;
  }
}

BlastWaveFitter* BlastWaveFitter::instance(){
  if(!mInstance) {
    mInstance = new BlastWaveFitter();
    gBlastWaveFitter = mInstance;
  }
  return mInstance;

}

BlastWaveFitter::~BlastWaveFitter(){
  delete mBlastWave;
  delete mMinuit;
  for(int ti=0; ti<mNData; ti++){
    delete mData[ti];
  }
  delete[] mData;
  delete[] mContinueCalcChi2;
}

void BlastWaveFitter::adjustDataSize(){
  if(mNData==mDataArraySize){
    cout << "WARNING Resize data array" << endl;
    DataForFit** tData = new DataForFit*[mNData];
    for(int ti=0; ti<mNData; ti++){
      tData[ti] = mData[ti];
    }
    delete[] mData;
    mDataArraySize*=2;
    mData = new DataForFit*[mDataArraySize];
    for(int ti=0; ti<mNData; ti++){
      mData[ti] = tData[ti];
    }
    delete[] tData;
    delete[] mContinueCalcChi2;
    mContinueCalcChi2 = new int[mDataArraySize];
    for(int ti=0; ti<mDataArraySize; ti++)  mContinueCalcChi2[ti]=1;
  }
}

void BlastWaveFitter::setParameters(double aT,
				    double aRho0){
  mNParameter = mNParameterMax;
  fixRho2();
  fixR();
  fixAs();
  fixTau();
  fixDeltat();
  mBlastWave->setParameters(aT,aRho0,10.,0.,10.,1.);
  setMinuitParameters();
}
void BlastWaveFitter::setParameters(double aT,
				    double aRho0,
				    double aRho2,
				    double aRho4,
				    double aAspectRatio){
  mNParameter = mNParameterMax;
  fixRx();
  fixAs();
  fixTau();
  fixDeltat();
  mBlastWave->setParameters(aT,aRho0,aRho2,aRho4,10.,aAspectRatio*10.,0.,10.,1.);
  setMinuitParameters();
}
void BlastWaveFitter::setParameters(double aT,
				    double aRho0,
				    double aR, 
				    double aAs,
				    double aTau, 
				    double aDeltat){
  mNParameter = mNParameterMax;
  fixRho2();
  useR();
  mBlastWave->setParameters(aT,aRho0,aR,aAs,aTau,aDeltat);
  setMinuitParameters();
}

void BlastWaveFitter::setParameters(double aT,
				    double aRho0, 
				    double aRho2,
				    double aRho4,
				    double aRx, 
				    double aRy, 
				    double aAs,
				    double aTau, 
				    double aDeltat){
  mNParameter = mNParameterMax;
  mBlastWave->setParameters(aT,aRho0,aRho2,aRho4,aRx,aRy,aAs,aTau,aDeltat);
  setMinuitParameters();
}


SpectraForFit* BlastWaveFitter::addSpectra(TGraphErrors* aSpectra, 
					   double aMass, 
					   double aPtMin, double aPtMax, int aStat){
  adjustDataSize();
  mData[mNData] = new SpectraForFit(aSpectra, aMass, aPtMin, aPtMax, aStat);
  mNData++;
  mNSpectra++;
  return (SpectraForFit*) mData[mNData-1];
}

SpectraForFit* BlastWaveFitter::addSpectra(TH1* aSpectra, 
					   double aMass, 
					   double aPtMin, double aPtMax, 
					   int aStat){
  adjustDataSize();
  mData[mNData] = new SpectraForFit(aSpectra, aMass, aPtMin, aPtMax, aStat);
  mNData++;
  return (SpectraForFit*) mData[mNData-1];
}

void BlastWaveFitter::addHbtRadii(TGraphErrors* aROut, 
				  TGraphErrors* aRSide, 
				  TGraphErrors* aRLong,
				  double aMass, int aStat){
  adjustDataSize();
  mData[mNData] = new ROutForFit(aROut, aMass, TMath::Pi()/4., aStat);
  mNData++;
  adjustDataSize();
  mData[mNData] = new RSideForFit(aRSide, aMass, TMath::Pi()/4., aStat);
  mNData++;
  adjustDataSize();
  mData[mNData] = new RLongForFit(aRLong, aMass, TMath::Pi()/4., aStat);
  mNData++;
  mBlastWave->setSpace(1);
  mSpace =1;
}

void BlastWaveFitter::addHbtRadiiSquareWRTReactionPlane(TGraphErrors* aROut, 
						  TGraphErrors* aRSide, 
						  TGraphErrors* aROutSide, 
						  TGraphErrors* aRLong, 
						  double aMass,
						  double aPhiP, int aStat ){
  if(aROut){
    adjustDataSize();
    mData[mNData] = new ROutSquareForFit(aROut, aMass, aPhiP, aStat);
    mNData++;
  }
  if(aRSide){
    adjustDataSize();
    mData[mNData] = new RSideSquareForFit(aRSide, aMass, aPhiP, aStat);
    mNData++;
  }
  if(aROutSide){
    adjustDataSize();
    mData[mNData] = new ROutSideSquareForFit(aROutSide, aMass, aPhiP, aStat);
    mNData++;
  }
  if(aRLong){
    adjustDataSize();
    mData[mNData] = new RLongSquareForFit(aRLong, aMass, aPhiP, aStat);
    mNData++;
  }
  mBlastWave->setSpace(1);
  mSpace =1;
}
void BlastWaveFitter::addHbtRadiiSquare(
  				  TGraphErrors* aROutSquare, 
			      TGraphErrors* aRSideSquare, 
			      TGraphErrors* aROutSideSquare,
			      TGraphErrors* aRLongSquare,  
			      double aMass, int aStat){
  adjustDataSize();
  mData[mNData] = new ROutSquareForFit(aROutSquare, aMass, aStat);
  mNData++;
  adjustDataSize();
  mData[mNData] = new RSideSquareForFit(aRSideSquare, aMass, aStat);
  mNData++;
  if(aROutSideSquare){
    adjustDataSize();
    mData[mNData] = new ROutSideSquareForFit(aROutSideSquare, aMass, aStat);
    mNData++;
  }
  adjustDataSize();
  mData[mNData] = new RLongSquareForFit(aRLongSquare, aMass, aStat);
  mNData++;
  mBlastWave->setSpace(1);
  mSpace =1;			      
}		      
void BlastWaveFitter::addHbtRadiiSquareOsc(
  				  TGraphErrors* aROutSquareOsc,
			      TGraphErrors* aRSideSquareOsc, 
			      TGraphErrors* aROutSideSquareOsc,
			      TGraphErrors* aRLongSquareOsc, 
			      double aMass, int aStat){
  adjustDataSize();
  mData[mNData] = new ROutSquareOscForFit(aROutSquareOsc, aMass, aStat);
  mNData++;
  adjustDataSize();
  mData[mNData] = new RSideSquareOscForFit(aRSideSquareOsc, aMass, aStat);
  mNData++;
  adjustDataSize();
  mData[mNData] = new ROutSideSquareOscForFit(aROutSideSquareOsc, aMass, aStat);
  mNData++;
  //adjustDataSize();
  //mData[mNData] = new RLongSquareOscForFit(aRLongSquareOsc, aMass, aStat);
  //mNData++;
  mBlastWave->setSpace(1);
  mSpace =1;					      
}


void BlastWaveFitter::addV2(TGraphErrors* aV2, 
			    double aMass, 
			    double aPtMin, double aPtMax, int aStat){
  adjustDataSize();
  mData[mNData] = new V2ForFit(aV2, aMass, aPtMin, aPtMax, aStat);
  mNData++;
}

void BlastWaveFitter::addV4(TGraphErrors* aV4, 
			    double aMass, 
			    double aPtMin, double aPtMax, int aStat){
  adjustDataSize();
  mData[mNData] = new V4ForFit(aV4, aMass, aPtMin, aPtMax, aStat);
  mNData++;
}


void BlastWaveFitter::fixT(){
  fixParameter(0);
}
void BlastWaveFitter::fixRho0(){
  fixParameter(1);
};
void BlastWaveFitter::fixRho2(){
  fixParameter(2);
};
void BlastWaveFitter::fixRho4(){
  fixParameter(3);
};
void BlastWaveFitter::fixRx(){
  fixParameter(4);
};
void BlastWaveFitter::fixRy(){
  fixParameter(5);
};
void BlastWaveFitter::fixAs(){
  fixParameter(6);
};
void BlastWaveFitter::fixTau(){
  fixParameter(7);
};
void BlastWaveFitter::fixDeltat(){
  fixParameter(8);
};
void BlastWaveFitter::useR(){
  fixRy();
  mBlastWave->useR();
}
void BlastWaveFitter::fixR(){
  fixRx();
  fixRy();
};


void BlastWaveFitter::fixParameter(int aIndex){
  if(! mParameterIsFixed[aIndex]){
    mNParameter --;
    mParameterIsFixed[aIndex] = 1;
    mMinuit->FixParameter(aIndex);
  }
}


void BlastWaveFitter::freeAllParameters(){
  mNParameter=mNParameterMax;
  for(int ti=0; ti<mNParameter; ti++){
    if(mParameterIsFixed[ti]==1){
      mMinuit->mnfree(-1*(ti+1));
      mParameterIsFixed[ti]=0;
    }
  }
  mBlastWave->doNotUseR();
}


void BlastWaveFitter::minimize(){
  double arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1000;
  arglist[1] = 10.;
  mMinuit->mnexcm("MINIMIZE", arglist ,2,ierflg);
  BWParameters* tPar = mBlastWave->par();
  mBlastWave->setSpace(mSpace);
  double tErrSig;
  mMinuit->GetParameter(0, tPar->T, tErrSig);
  mMinuit->GetParameter(1, tPar->Rho0, tErrSig);
  mMinuit->GetParameter(2, tPar->Rho2, tErrSig);
  mMinuit->GetParameter(3, tPar->Rho4, tErrSig);
  mMinuit->GetParameter(4, tPar->Rx, tErrSig);
  mMinuit->GetParameter(5, tPar->Ry, tErrSig);
  mMinuit->GetParameter(6, tPar->As, tErrSig);
  mMinuit->GetParameter(7, tPar->Tau, tErrSig);
  mMinuit->GetParameter(8, tPar->Deltat, tErrSig);
}

void BlastWaveFitter::calculatePreciseErrors(){
  double arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1000;
  arglist[1] = 10.;
  mMinuit->mnexcm("MINO", arglist ,1,ierflg);
}

TGraph* BlastWaveFitter::contour(char* aName, int aNSigma, int i1, int i2){
  freeAllParameters();
  for(int ti=0; ti<mNParameterMax; ti++){
    if(ti!=i1 && ti!=i2) fixParameter(ti);
  }
  double tNDOF = nDOF();
  double tChi2PerDof = tNDOF>0.? chi2()/tNDOF : 1.;
  if(tChi2PerDof<1.) tChi2PerDof = 1.;
  mMinuit->SetErrorDef(aNSigma*aNSigma*tChi2PerDof);
  TGraph* tG = (TGraph*) mMinuit->Contour(20,i1,i2);  
  mMinuit->SetErrorDef(1.);
  char tGName[50];
  sprintf(tGName,"%sNSig%i\0",aName,aNSigma);
  if(mMinuit->GetStatus()!=0){
   cerr << "Problem to find contour " << mMinuit->GetStatus() << endl;
   cerr << "Make dummy " << tGName << endl;
   double tX[1]; tX[0]=0.;
    double tY[1]; tY[0]=0.;  
    tG = new TGraph(1,tX,tY);
  }  
  tG->SetName(tGName);     	  
  return tG;
}

TGraph* BlastWaveFitter::betaTContour(int aNSigma){
  return contour("betaT",aNSigma,1,0);
}

TGraph* BlastWaveFitter::S2VsRho2Contour(int aNSigma){
  TGraph* tG = contour("S2Rhoa",aNSigma,2,4);
  int tN = tG->GetN();
  double* tY = tG->GetY();
  double tVal;
  double tErr;
  mMinuit->GetParameter(5,tVal,tErr);
  for(int ti=0; ti<tN; ti++){
     tY[ti] /= tVal;
  }  
  TGraph* ttG = new TGraph(tN, tY, tG->GetX());
  char tGName[20];
  strcpy(tGName,tG->GetName());
  tG->Delete();
  ttG->SetName(tGName);
  return ttG;
}

TGraph* BlastWaveFitter::TauVsDeltatContour(int aNSigma){
  return contour("TauDt",aNSigma,7,8);
}

void BlastWaveFitter::coutFitResults(int ScaleErrorByChi2PerDof){
  cout << "Fit Results" << endl;
  double tChi2=0.;
  double tNDOF=0;
  for(int ti=0; ti<mNData; ti++){
    double ttNDOF = mData[ti]->nDOF();
    double ttChi2 = mData[ti]->chi2With(mBlastWave,1);
    tNDOF += ttNDOF;
    tChi2 += ttChi2;
    cout << mData[ti]->getName() << " " << ttChi2 << " / " << ttNDOF 
           << " CL = " << TMath::Prob(ttChi2, (int) ttNDOF) <<  endl;
  }
  tNDOF -= mNParameter;
  cout << "TOTAL Chi2/dof = " << tChi2 << "/" << tNDOF 
         << " CL = " << TMath::Prob(tChi2, (int) tNDOF) << endl;
  double tErrScale = ScaleErrorByChi2PerDof? 0. : sqrt(tChi2/tNDOF);
  double tVal;
  double tErr;
  mMinuit->GetParameter(0,tVal,tErr);
  cout << "T (MeV) = " << tVal*1000. << " +- " << tErr*tErrScale*1000. << endl;
  mMinuit->GetParameter(1,tVal,tErr);
  cout << "Rho0 = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(2,tVal,tErr);
  cout << "Rho2 = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(3,tVal,tErr);
  cout << "Rho4 = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(4,tVal,tErr);
  cout << "Rx (fm) = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(5,tVal,tErr);
  cout << "Ry (fm) = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(6,tVal,tErr);
  cout << "As = " << tVal << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(7,tVal,tErr);
  cout << "Tau (fm/c) = " << fabs(tVal) << " +- " << tErr*tErrScale << endl;
  mMinuit->GetParameter(8,tVal,tErr);
  cout << "Deltat (fm/c) = " << tVal << " +- " << tErr*tErrScale << endl;

 //test chi2 map

cout<<" jdu delat chi2 mapu"<<endl;

TH2D *scan01 = new TH2D("scan01", "Chi2/ndf counter: Temperature vs #rho; #rho; Temperature", 60, 0.7,1.3,
                                                                                              60, 0.075,0.135);

TH2D *scan02 = new TH2D("scan02", "Chi2/ndf counter: Temperature vs R; R; Temperature", 60, 7.,13.,
                                                                                        60, 0.075,0.135);

TH2D *scan03 = new TH2D("scan03", "Chi2/ndf counter: Temperature vs #tau; #tau; Temperature", 60, 6.5,9.5,
                                                                                              60, 0.075,0.135);

TH2D *scan04 = new TH2D("scan04", "Chi2/ndf counter: Temperature vs #delta #tau; #delta #tau; Temperature", 60, 2.,3.,
                                                                                                        60,0.075,0.135);
/////
TH2D *scan05 = new TH2D("scan05", "Chi2/ndf counter: #rho vs R; R; #rho", 60, 7., 13.,
                                                                          60, 0.7,1.3);

TH2D *scan06 = new TH2D("scan06", "Chi2/ndf counter: #rho vs #tau; #tau; #rho", 60, 6.5, 9.5,
                                                                                60, 0.7,1.3);

TH2D *scan07 = new TH2D("scan07", "Chi2/ndf counter: #rho vs #delta #tau; #delta #tau; #rho", 60, 2., 3.,
                                                                                              60, 0.7,1.3);
///
TH2D *scan08 = new TH2D("scan08", "Chi2/ndf counter: R vs #tau; #tau; R", 60, 6.5, 9.5, 
                                                                          60, 7., 13.);

TH2D *scan09 = new TH2D("scan09", "Chi2/ndf counter: R vs #delta #tau; #delta #tau; R", 60, 2., 3., 
                                                                                        60, 7., 13.);
//
TH2D *scan10 = new TH2D("scan10", "Chi2/ndf counter: #tau vs #delta #tau; #delta #tau;#tau", 60, 2., 3., 
                                                                                             60, 6.5, 9.5);

Double_t par0, par1, par2, par3, par4, par5, par6, par7, par8;
Double_t par0Err, par1Err, par2Err, par3Err, par4Err, par5Err, par6Err, par7Err, par8Err;
mMinuit->GetParameter(0,par0,par0Err); // T
mMinuit->GetParameter(1,par1,par1Err); // Rho0
mMinuit->GetParameter(2,par2,par2Err); // Rho2
mMinuit->GetParameter(3,par3,par3Err); // Rho4
mMinuit->GetParameter(4,par4,par4Err); // Rx
mMinuit->GetParameter(5,par5,par5Err); // Ry
mMinuit->GetParameter(6,par6,par6Err); // As
mMinuit->GetParameter(7,par7,par7Err); // Tau
mMinuit->GetParameter(8,par8,par8Err); // Delta Tau

Double_t p0p1p2[9] = {0.,0., par2, par3, par4, par5, par6, par7,par8};
Double_t calcChi2 = 0;
cout<<"toto je teplota\t"<<par0<<endl;
/*
Double_t *parameters=NULL;
gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
Double_t ndf1 = 65;

p0p1p2[0] = 0.;
p0p1p2[1] = 0.;
p0p1p2[2] = par2;
p0p1p2[3] = par3;
p0p1p2[4] = par4;
p0p1p2[5] = par5;
p0p1p2[6] = par6;
p0p1p2[7] = par7;
p0p1p2[8] = par8;

for(Double_t p1 = 0.075; p1 <= 0.135; p1 += (0.06/60.)){
  p0p1p2[0] = p1;
  for(Double_t p2 = 0.7; p2 <= 1.3; p2 += (0.6/60.)){
    p0p1p2[1] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan01->Fill(p2, p1, calcChi2);
  }
}


p0p1p2[0] = 0.;
p0p1p2[1] = par1;
p0p1p2[2] = par2;
p0p1p2[3] = par3;
p0p1p2[4] = 0.;
p0p1p2[5] = par5;
p0p1p2[6] = par6;
p0p1p2[7] = par7;
p0p1p2[8] = par8;

for(Double_t p1 = 0.075; p1 <= 0.135; p1 += (0.06/60.)){ 
  p0p1p2[0] = p1;
  for(Double_t p2 = 7.; p2 <= 13.; p2 += (6./60.)){
    p0p1p2[4] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan02->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = 0.;
p0p1p2[1] = par1;
p0p1p2[2] = par2;
p0p1p2[3] = par3;
p0p1p2[4] = par4;
p0p1p2[5] = par5;
p0p1p2[6] = par6;
p0p1p2[7] = 0.;
p0p1p2[8] = par8;

for(Double_t p1 = 0.075; p1 <= 0.135; p1 += (0.06/60.)){ 
  p0p1p2[0] = p1;
  for(Double_t p2 = 6.5; p2 <= 9.5; p2 += (3./60.)){
    p0p1p2[7] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan03->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = 0.;
p0p1p2[1] = par1;
p0p1p2[2] = par2;
p0p1p2[3] = par3;
p0p1p2[4] = par4;
p0p1p2[5] = par5;
p0p1p2[6] = par6;
p0p1p2[7] = par7;
p0p1p2[8] =  0.;

for(Double_t p1 = 0.075; p1 <= 0.135; p1 += (0.06/60.)){ 
  p0p1p2[0] = p1;
  for(Double_t p2 = 2.; p2 <= 3.; p2 += (1./60.)){
    p0p1p2[8] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan04->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  0.;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  0.;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6;
p0p1p2[7] =  par7;
p0p1p2[8] = par8;
for(Double_t p1 = 0.7; p1 <= 1.3; p1 += (0.6/60.)){ 
  p0p1p2[1] = p1;
  for(Double_t p2 = 7.; p2 <= 13.; p2 += (6./60.)){
    p0p1p2[4] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan05->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  0.;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  par4;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6;
p0p1p2[7] =  0.;
p0p1p2[8] = par8;
for(Double_t p1 = 0.7; p1 <= 1.3; p1 += (0.6/60.)){ 
  p0p1p2[1] = p1;
  for(Double_t p2 = 6.5; p2 <= 9.5; p2 += (3./60.)){
    p0p1p2[7] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan06->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  0.;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  par4;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6;
p0p1p2[7] =  par7;
p0p1p2[8] =  0.;
for(Double_t p1 = 0.7; p1 <= 1.3; p1 += (0.6/60.)){ 
  p0p1p2[1] = p1;
  for(Double_t p2 = 2.; p2 <= 3.; p2 += (1./60.)){
    p0p1p2[8] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan07->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  par1;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  0.;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6;
p0p1p2[7] =  0.;
p0p1p2[8] =  par8;
for(Double_t p1 = 7.; p1 <= 13.; p1 += (6./60.)){ 
  p0p1p2[4] = p1;
  for(Double_t p2 = 6.5; p2 <= 9.5; p2 += (3./60.)){
    p0p1p2[7] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan08->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  par1;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  0.;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6;
p0p1p2[7] =  par7;
p0p1p2[8] =  0.;
for(Double_t p1 = 7.; p1 <= 13.; p1 += (6./60.)){ 
  p0p1p2[4] = p1;
  for(Double_t p2 = 2.; p2 <= 3.; p2 += (1./60.)){
    p0p1p2[8] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan09->Fill(p2, p1, calcChi2);
  }
}

p0p1p2[0] = par0;
p0p1p2[1] =  par1;
p0p1p2[2] =  par2;
p0p1p2[3] =  par3;
p0p1p2[4] =  par4;
p0p1p2[5] =  par5;
p0p1p2[6] =  par6; 
p0p1p2[7] = 0.;
p0p1p2[8] =  0.;


for(Double_t p1 = 6.5; p1 <= 9.5; p1 += (3./60.)){ 
  p0p1p2[7] = p1;
  for(Double_t p2 = 2.; p2 <= 3.; p2 += (1./60.)){
    p0p1p2[8] = p2;
    gMinuit->Eval(0, parameters, calcChi2, p0p1p2, 3);
    calcChi2=calcChi2/ndf1;
    scan10->Fill(p2, p1, calcChi2);
  }
}
*/
  /*
gMinuit->SetErrorDef(1);
   TGraph *gr1 = (TGraph*)gMinuit->Contour(80,1,0);
   gr1->SetFillColor(38);

gMinuit->SetErrorDef(2);
   TGraph *gr2 = (TGraph*)gMinuit->Contour(80,1,0);
   gr2->SetFillColor(38);

gMinuit->SetErrorDef(3);
   TGraph *gr3 = (TGraph*)gMinuit->Contour(80,1,0);
   gr3->SetFillColor(38);


gMinuit->SetErrorDef(4);
   TGraph *gr4 = (TGraph*)gMinuit->Contour(80,1,0);
   gr4->SetFillColor(38);


gMinuit->SetErrorDef(5);
   TGraph *gr5 = (TGraph*)gMinuit->Contour(80,1,0);
   gr5->SetFillColor(38);
*/


TString outputscan;
outputscan="scan.root";

   TFile *oFile=new TFile(outputscan.Data(), "RECREATE"); 

   scan01->Write();
   scan02->Write();
   scan03->Write();
   scan04->Write();
   scan05->Write();
   scan06->Write();
   scan07->Write();
   scan08->Write();
   scan09->Write();
   scan10->Write();
  /* gr1->Write();
   gr2->Write();
   gr3->Write();
   gr4->Write();
   gr5->Write();
*/


   oFile->Close();

cout<<"udelal jsem chi2 mapu !!"<<endl;
  //test chi2 map
}

TH2D* BlastWaveFitter::confLevelMap(char* aName, 
				    int ai1, int aNBinX,
				    int ai2, int aNBinY){
  double tVal1, tErr1, tVal2, tErr2;
  mMinuit->GetParameter(ai1,tVal1,tErr1);
  mMinuit->GetParameter(ai2,tVal2,tErr2);
  double tSqrtChi2Perdof = sqrt(chi2()/nDOF());
  tErr1 *= tSqrtChi2Perdof;
  tErr2 *= tSqrtChi2Perdof;
  return confLevelMap(aName,
		      ai1,aNBinX,tVal1-tErr1*3.5,tVal1+tErr1*3.5,
		      ai2,aNBinY,tVal2-tErr2*3.5,tVal2+tErr2*3.5);

}

TH2D* BlastWaveFitter::confLevelMap(char* aName, 
				int ai1, int aNBinX, double aXMin, double aXMax, 
  	int ai2, int aNBinY, double aYMin, double aYMax){
	  	
  double* tPrevPar = (double*) mBlastWave->par();	  	
  double* tPar = new double[mNParameterMax];
  for(int ti=0; ti<mNParameterMax; ti++) tPar[ti]=tPrevPar[ti];

  double tMinimum = chi2()/nDOF();
  TH2D* tH = new TH2D(aName,aName,aNBinX,aXMin,aXMax,aNBinY,aYMin,aYMax);
  for(int ti=1; ti<=aNBinX;ti++){
    for(int tj=1; tj<=aNBinY; tj++){
      tPar[ai1] = tH->GetXaxis()->GetBinCenter(ti);
      tPar[ai2] = tH->GetYaxis()->GetBinCenter(tj);
      mBlastWave->setParameters(tPar);
      tH->SetCellContent(ti,tj,chi2()/nDOF()-tMinimum);
    }
  }
  mBlastWave->setParameters(tPrevPar);
  delete[] tPar;	 
  return tH; 	
}

void BlastWaveFitter::clearData(){
  for(int ti=0; ti<mNData; ti++){
    delete mData[ti];
  }
  mNData = 0;
  mNSpectra =0;
}

void BlastWaveFitter::getParameter(Int_t parNo, Double_t& currentValue, Double_t& currentError){
  mMinuit->GetParameter(parNo, currentValue, currentError);
}


