#include "IOHandler.h"
#include "BlastWaveFitter.h"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TH2.h"

#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

double IOHandler::MPi = 0.1396; 
double IOHandler::MK = 0.4937;
double IOHandler::MKStar = 0.892;
double IOHandler::MP = 0.9383;
double IOHandler::MLa = 1.1157;
double IOHandler::MXi = 1.32132;
double IOHandler::MOm = 1.67245;
double IOHandler::MPhi = 1.019413;

IOHandler::IOHandler(BlastWaveFitter* aBWFitter, 
		     const char* aInputFileName)
  : mBWFitter(aBWFitter), 
    mOutputFile(0) {
  mCent=-1;
  mNCent = 0;

  mHBT=0;
  mV2=0;
  mSpectra=0;
  mNParticles=0;
  mFixedAs = 1;

  if(aInputFileName){
    mInputFile = new TFile(aInputFileName);
  }
  else{
    mInputFile = 0;
  }


}

void IOHandler::init(){
  initLoad();

  mValidParticle = new int[mNParticles];
  mSpectras = new SpectraForFit*[mNParticles];
  mYield = new double*[mNParticles];
  for(int ti=0; ti<mNParticles; ti++){
    mSpectras[ti]=0;
    mValidParticle[ti] = 1;
    mYield[ti] = new double[mNCent];
  }

  mChi2 = new double[mNCent];
  mConfLevel = new double[mNCent];

  mT = new double[mNCent];
  mTErr = new double[mNCent];
  mRho0 = new double[mNCent];
  mRho0Err = new double[mNCent];
  mRho2 = new double[mNCent];
  mRho2Err = new double[mNCent];
  mRho4 = new double[mNCent];
  mRho4Err = new double[mNCent];
  mRx = new double[mNCent];
  mRxErr = new double[mNCent];
  mRy = new double[mNCent];
  mRyErr = new double[mNCent];
  mAs = new double[mNCent];
  mAsErr = new double[mNCent];
  mTau = new double[mNCent];
  mTauErr = new double[mNCent];
  mDeltat = new double[mNCent];
  mDeltatErr = new double[mNCent];
}

IOHandler::~IOHandler(){
  if(mSpectra){
    for(int ti=0; ti<mNParticles; ti++){
      delete[] mYield[ti];
    }
    delete[] mYield;
    delete[] mSpectras;
  }
  delete[] mValidParticle;
}

void IOHandler::close(){
  if(mInputFile) mInputFile->Close();
  mOutputFile->Close();
}

void IOHandler::setOutputFileName(const char* aFileName){
  setOutputFile(aFileName);
}
void IOHandler::setOutputFile(const char* aFileName){
  if(mOutputFile) mOutputFile->Close();
  cout << "Fit results will be save in " << aFileName << endl;
  mOutputFile = new TFile(aFileName,"RECREATE");
}
void IOHandler::setOutputFile(TFile* aFile){
  if(mOutputFile) mOutputFile->Delete();
  mOutputFile = aFile;
}
  
void IOHandler::nextCent(){
  mCent++;
}
int IOHandler::nextFit(){
  mCent++;
  if(mCent<mNCent){
    mBWFitter->clearData();
    fixParameters();
    if(mSpectra){
      for(int ti=0; ti<mNParticles; ti++){
	mSpectras[ti] = 0;
      }
    }
    load();
    return 1;
  }
  else{
    return 0;
  }
}



void IOHandler::storeParameters(int ScaleErrorByChi2PerDof){
  mChi2[mCent] = mBWFitter->chi2();
  int tNDof = mBWFitter->nDOF();
  mConfLevel[mCent] = TMath::Prob(mChi2[mCent],tNDof-1);
  mChi2[mCent]/=tNDof;
  double tErrScale = ScaleErrorByChi2PerDof? sqrt(mChi2[mCent]) : 1.;
  mBWFitter->getParameter(0,mT[mCent],mTErr[mCent]);
  mTErr[mCent]*=tErrScale;
  mBWFitter->getParameter(1,mRho0[mCent],mRho0Err[mCent]);
  mRho0Err[mCent]*=tErrScale;
  mBWFitter->getParameter(2,mRho2[mCent],mRho2Err[mCent]);
  mRho2Err[mCent]*=tErrScale;
  mBWFitter->getParameter(3,mRho4[mCent],mRho4Err[mCent]);
  mRho4Err[mCent]*=tErrScale;
  mBWFitter->getParameter(4,mRx[mCent],mRxErr[mCent]);
  mRxErr[mCent]*=tErrScale;
  mBWFitter->getParameter(5,mRy[mCent],mRyErr[mCent]);
  mRyErr[mCent]*=tErrScale;
  mBWFitter->getParameter(6,mAs[mCent],mAsErr[mCent]);
  mAsErr[mCent]*=tErrScale;
  mBWFitter->getParameter(7,mTau[mCent],mTauErr[mCent]);
  mTauErr[mCent]*=tErrScale;
  mBWFitter->getParameter(8,mDeltat[mCent],mDeltatErr[mCent]);
  mDeltatErr[mCent]*=tErrScale;
}

void IOHandler::writeSpectra(){
  mOutputFile->cd();
  if(mSpectra){
    for(int ti=0; ti<mNParticles; ti++){
      if(mSpectras[ti]){
	char tName[50];
	sprintf(tName,"Spectra%s%i",mParticleName[ti],mCent);
	mSpectras[ti]->histo(tName, 100, mBWFitter->blastWave())->Write();
	char ttName[50];
	sprintf(ttName,"SpectraW%s%i",mParticleName[ti],mCent);
	mSpectras[ti]->histo(ttName, 100, 0.05, 3.05, 
			     mBWFitter->blastWave())->Write();
	mYield[ti][mCent] = mSpectras[ti]->normBW();
      }      
    }
    //char ttName[50];
    //sprintf(ttName,"TVsRho0Cont%i",mCent);
    //mBWFitter->confLevelMap(ttName,0,20,1,20)->Write();
  }
}

void IOHandler::writeHBT(){
  if(mHBT){
    mOutputFile->cd();
    for(int ti=0; ti<mNParticles; ti++){
      if(mValidParticle[ti]){
	char tName[50];
	sprintf(tName,"%s%i",mParticleName[ti],mCent);
	mBWFitter->blastWave()->rOut(tName,mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rSide(tName,mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rLong(tName,mParticleMass[ti])->Write();

	mBWFitter->blastWave()->rOutSquare(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rSideSquare(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rLongSquare(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rOutSideSquare(tName, mParticleMass[ti])->Write();

	mBWFitter->blastWave()->rOutSquareOsc(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rSideSquareOsc(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rLongSquareOsc(tName, mParticleMass[ti])->Write();
	mBWFitter->blastWave()->rOutSideSquareOsc(tName, mParticleMass[ti])->Write();
      }
    }
    if(mV2){
      double tPt[] = {0.2,0.4,0.6};
      for(int ti=0; ti<3; ti++){
	char tNameOut[50];
	sprintf(tNameOut,"PiROut2Pt%iC%i",(int) tPt[ti]*100, mCent);
	mBWFitter->blastWave()->rOutSquareVsPhiP(tNameOut,MPi,tPt[ti]);
	char tNameSide[50];
	sprintf(tNameSide,"PiRSide2Pt%iC%i",(int) tPt*100, mCent);
	mBWFitter->blastWave()->rSideSquareVsPhiP(tNameSide,MPi,tPt[ti]);
	char tNameLong[50];
	sprintf(tNameLong,"PiRLong2Pt%iC%i",(int) tPt*100, mCent);
	mBWFitter->blastWave()->rLongSquareVsPhiP(tNameLong,MPi,tPt[ti]);
	char tNameOutSide[50];
	sprintf(tNameOutSide,"PiROutSide2Pt%iC%i",(int) tPt*100, mCent);
	mBWFitter->blastWave()->rOutSideSquareVsPhiP(tNameOutSide,MPi,tPt[ti]);
      }
    }
  }
}

void IOHandler::writeEllipticFlow(){
  if(mV2){
    mOutputFile->cd();
    for(int ti=0; ti<mNParticles; ti++){
      if(mValidParticle[ti]){
	char tName[50];
	sprintf(tName,"%s%i",mParticleName[ti],mCent);
	mBWFitter->blastWave()->v2VsPt(tName,mParticleMass[ti])->Write();
	mBWFitter->blastWave()->v4VsPt(tName,mParticleMass[ti])->Write();
	mBWFitter->blastWave()->v6VsPt(tName,mParticleMass[ti])->Write();
      }
    }
    if(mSpectra){
      char ttName[50];
      sprintf(ttName,"TVsRxCont%i",mCent);
      //mBWFitter->confLevelMap(ttName,0,20,3,20)->Write();
    }
  }
}

void IOHandler::writeContours(){
  cout << "Calculate contours" << endl;
  if(mSpectra){
    char tName11[50];
    sprintf(tName11,"TRho0Cont%i",mCent);
    mBWFitter->contour(tName11,1,0,1)->Write();
    char tName12[50];
    sprintf(tName12,"TRho0Cont%i",mCent);
    mBWFitter->contour(tName12,2,0,1)->Write();
    if(mV2){
      char tName21[50];
      sprintf(tName21,"TRxCont%i",mCent);
      mBWFitter->contour(tName21,1,0,4)->Write();
      char tName22[50];
      sprintf(tName22,"TRxCont%i",mCent);
      mBWFitter->contour(tName22,2,0,4)->Write();
    }
  }
  fixParameters();
}

void IOHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  if(mFixedAs) mBWFitter->fixAs(); // need an extra switch
  mBWFitter->fixRho4();
  if(!mHBT){

    mBWFitter->fixRy();
    mBWFitter->fixTau();
    mBWFitter->fixDeltat();
  }
  if(!mV2){
    mBWFitter->useR();
    mBWFitter->fixRho2();
 //   mBWFitter->fixT();
 //   mBWFitter->fixRho0();
 //   mBWFitter->fixRx();
 //   mBWFitter->fixTau();
 //   mBWFitter->fixDeltat();
  }
  if(!(mV2 || mHBT)){
    mBWFitter->fixRx();
  }
}


void IOHandler::writeParameters(){
  mOutputFile->cd();
  TGraph* tGChi2 = new TGraph(mNCent,mNParticipant,mChi2);
  tGChi2->SetName("Chi2");
  tGChi2->Write();
  TGraph* tGConfLevel = new TGraph(mNCent,mNParticipant,mConfLevel);
  tGConfLevel->SetName("ConfLevel");
  tGConfLevel->Write();
  if(!mBWFitter->parameterIsFixed(0)){
    TGraphErrors* tGT = new TGraphErrors(mNCent,mNParticipant,mT,0,mTErr);
    tGT->SetName("T");
    tGT->Write();
  }
  if(!mBWFitter->parameterIsFixed(1)){
    TGraphErrors* tGRho0 = new TGraphErrors(mNCent,mNParticipant,mRho0,
					  0,mRho0Err);
    tGRho0->SetName("Rho0");
    tGRho0->Write();
    double* tMeanBetaT = new double[mNCent];
    double* tMeanBetaTErr = new double[mNCent];
    for(int ti=0; ti<mNCent; ti++){
      Rho0ToMeanBetaT(mRho0[ti],mRho0Err[ti],
		      tMeanBetaT[ti],tMeanBetaTErr[ti]);
    }
    TGraphErrors* tGBetaT = new TGraphErrors(mNCent,mNParticipant,tMeanBetaT,
					     0,tMeanBetaTErr);
    tGBetaT->SetName("BetaT");
    tGBetaT->Write();
    if(!mBWFitter->parameterIsFixed(2)){
      TGraphErrors* tGRho2 = new TGraphErrors(mNCent,mNParticipant,mRho2,
					      0,mRho2Err);
      tGRho2->SetName("Rho2");
      tGRho2->Write();

      /*
      double* tMeanBetaTOsc = new double[mNCent];
      double* tMeanBetaTOscErr = new double[mNCent];
      for(int ti=0; ti<mNCent; ti++){
	double tErr = sqrt(mRho0Err[ti]*mRho0Err[ti]+
			   mRho2Err[ti]*mRho2Err[ti]);
	Rho0ToMeanBetaT(mRho0[ti]+mRhoA[ti],tErr,
			tMeanBetaTOsc[ti],tMeanBetaTOscErr[ti]);
	tMeanBetaTOsc[ti]-=tMeanBetaT[ti];
	tErr = tMeanBetaTOscErr[ti];
	tMeanBetaTOscErr[ti] = sqrt(tErr*tErr+
				    mRho0Err[ti]*mRho0Err[ti]);
      }
      TGraphErrors* tGBetaTOsc = new TGraphErrors(mNCent,mNParticipant,
						  tMeanBetaTOsc,
						  0,tMeanBetaTOscErr);
      tGBetaTOsc->SetName("BetaTOsc");
      tGBetaTOsc->Write();
      delete[] tMeanBetaTOsc;
      delete[] tMeanBetaTOscErr;
      */
    }
    delete[] tMeanBetaT;
    delete[] tMeanBetaTErr;
  }  
  if(!mBWFitter->parameterIsFixed(3)){
    TGraphErrors* tGRho4 = new TGraphErrors(mNCent,mNParticipant,mRho4,
					  0,mRho4Err);
    tGRho4->SetName("Rho4");
    tGRho4->Write();
  }
  if(!mBWFitter->parameterIsFixed(4)){
    TGraphErrors* tGRx = new TGraphErrors(mNCent,mNParticipant,mRx,
					    0,mRxErr);
    tGRx->SetName("Rx");
    tGRx->Write();
  }
  if(!mBWFitter->parameterIsFixed(5)){
    TGraphErrors* tGRy = new TGraphErrors(mNCent,mNParticipant,mRy,
					    0,mRyErr);
    tGRy->SetName("Ry");
    tGRy->Write();
  }
  if(!mBWFitter->parameterIsFixed(6)){
    TGraphErrors* tGAs = new TGraphErrors(mNCent,mNParticipant,mAs,
					    0,mAsErr);
    tGAs->SetName("As");
    tGAs->Write();
  }
  if(!mBWFitter->parameterIsFixed(7)){
    TGraphErrors* tGTau = new TGraphErrors(mNCent,mNParticipant,mTau,
					    0,mTauErr);
    tGTau->SetName("Tau");
    tGTau->Write();
  }
  if(!mBWFitter->parameterIsFixed(8)){
    TGraphErrors* tGDeltat = new TGraphErrors(mNCent,mNParticipant,mDeltat,
					    0,mDeltatErr);
    tGDeltat->SetName("Deltat");
    tGDeltat->Write();
  }
  if(mSpectra){
    for(int ti=0; ti<mNParticles; ti++){
      if(mSpectras[ti]){
	char tName[50];
	sprintf(tName,"Yield%s",mParticleName[ti]);
	TGraph* tGYield = new TGraph(mNCent,mNParticipant,mYield[ti]);
	tGYield->SetName(tName);
	tGYield->Write();
      } 
    }
  }
}


void IOHandler::Rho0ToMeanBetaT(double aRho0, double aRho0Err, 
				double& aBetaT, double& aBetaTErr){
  aBetaT=0.;
  aBetaTErr=0.;
  double tR;
  double tTanH;
  for(int ti=0; ti<100; ti++){
    tR= (ti+0.5)/100.;
    tTanH = TMath::TanH(aRho0*tR);
    aBetaT+= tR*tTanH;
    aBetaTErr+= tR*(1.-tTanH*tTanH);
  }
  aBetaT/= 50.; // 2/NBin
  aBetaTErr*= (aRho0*aRho0Err/50.);
}

void IOHandler::writeEStruct(){
  if(mV2){ // makes no sense otherwise
    mOutputFile->cd();
    for(int ti=0; ti<mNParticles; ti++){
      if(mValidParticle[ti]){
	char tName[50];
	sprintf(tName,"%s%i",mParticleName[ti],mCent);
	mBWFitter->blastWave()->meanPtVsPhi(tName,mParticleMass[ti])->Write();
	mBWFitter->blastWave()->autoCorrVsPhi(tName,mParticleMass[ti])->Write();
      }
    }
  }
}
