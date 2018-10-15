#include "HbtPRCDataHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
#include <iostream>
using namespace std;
HbtPRCDataHandler::HbtPRCDataHandler(BlastWaveFitter* aBWFitter, 
				     const char* aInputFileName, 
				     int aAsHbt, // 0=no, 1=Fourier, 2=Broken 
				     int aV2,
				     int aSpectraFit) // 0=no fit, 1= fit, 2 = fit only
  //:STAR200asHBTSpectraHandler(aBWFitter,aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectraFit = aSpectraFit;
  mSpectra= mSpectraFit!=0? 1 : 0;
  mAsHbt = aAsHbt;
  if(mSpectraFit!=2){
    mV2= (aAsHbt && aV2);
    mHBT=1;
  }
}

void HbtPRCDataHandler::initLoad(){
  switch(mSpectraFit){
  case 0:
    switch(mAsHbt){
    case 0:
      {
	setOutputFile("data/HbtPRC/IntHbt.fit.root\0");
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/IntHbtErrScaled.root\0");      
	TFile fPar("data/HbtPRC/SpectraC6.fit.root");
	mFixedT = ((TGraphErrors*) fPar.Get("T"))->GetY();
	mFixedRho0 = ((TGraphErrors*) fPar.Get("Rho0"))->GetY();
	fPar.Close();
      }
      break;
    case 1:
      {
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
	setOutputFile("data/HbtPRC/AsHbtFourier0.fit.root\0");
	TFile fPar("data/HbtPRC/SpectraC5.fit.root");
	mFixedT = ((TGraphErrors*) fPar.Get("T"))->GetY();
	mFixedRho0 = ((TGraphErrors*) fPar.Get("Rho0"))->GetY();
	fPar.Close();
      }
      break;
    case 2:
      {
	//if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/PRLHbtVsPhiPErrScaled.root\0");
	if(mV2){
	  setOutputFile("data/HbtPRC/AsHbtPhiV2.fit.root\0");
	}
	else{
	  setOutputFile("data/HbtPRC/AsHbtPhi.fit.root\0");
	}
	TFile fPar("data/HbtPRC/SpectraC5.fit.root");
	mFixedT = ((TGraphErrors*) fPar.Get("T"))->GetY();
	mFixedRho0 = ((TGraphErrors*) fPar.Get("Rho0"))->GetY();
	fPar.Close();
      }
      break;
    case 3:
      {
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
	if(mV2){
	  setOutputFile("data/HbtPRC/AsHbtFourierV2ErrScaled.fit.root\0");
	}
	else{
	  setOutputFile("data/HbtPRC/AsHbtFourierErrScaled.fit.root\0");
	}
	TFile fPar("data/HbtPRC/SpectraC5.fit.root");
	mFixedT = ((TGraphErrors*) fPar.Get("T"))->GetY();
	mFixedRho0 = ((TGraphErrors*) fPar.Get("Rho0"))->GetY();
	fPar.Close();
      }
      break;
    case 4:
      {
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
	if(mV2){
	  setOutputFile("data/HbtPRC/AsHbtRcPhiV2.fit.root\0");
	}
	else{
	  setOutputFile("data/HbtPRC/AsHbtRcPhi.fit.root\0");
	}
	TFile fPar("data/HbtPRC/SpectraC5.fit.root");
	mFixedT = ((TGraphErrors*) fPar.Get("T"))->GetY();
	mFixedRho0 = ((TGraphErrors*) fPar.Get("Rho0"))->GetY();
	fPar.Close();
      }
      break;
    }
    break;
  case 1:
    switch(mAsHbt){
    case 0:
      setOutputFile("data/HbtPRC/IntHbtSpectra.fit.root\0");
      if(!mInputFile) mInputFile = new TFile("data/HbtPRC/IntHbtErrScaled.root\0");
      break;
    case 1:
      if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
      setOutputFile("data/HbtPRC/AsHbtFourier0Spectra.fit.root\0");
      break;
    case 2:
      //if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
	if(!mInputFile) mInputFile = new TFile("data/HbtPRC/PRLHbtVsPhiPErrScaled.root\0");
      if(mV2){
	setOutputFile("data/HbtPRC/AsHbtPhiSpectraV2.fit.root\0");
      }
      else{
	setOutputFile("data/HbtPRC/AsHbtPhiSpectra.fit.root\0");
      }
      break;
    case 3:
      if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
      if(mV2){
	setOutputFile("data/HbtPRC/AsHbtFourierSpectraV2.fit.root\0");
      }
      else{
	setOutputFile("data/HbtPRC/AsHbtFourierSpectra.fit.root\0");
      }
      break;
    case 4:
      if(!mInputFile) mInputFile = new TFile("data/HbtPRC/AsHbtPRLFourierErrScaled.root\0");
      if(mV2){
	setOutputFile("data/HbtPRC/AsHbtRcPhiSpectraV2.fit.root\0");
      }
      else{
	setOutputFile("data/HbtPRC/AsHbtRcPhiSpectra.fit.root\0");
      }
      break;
    }
    break;
  case 2:
    mInputFile=0;
    switch(mAsHbt){
    case 0:
      setOutputFile("data/HbtPRC/SpectraC6.fit.root\0");
      break;
    case 1:
      setOutputFile("data/HbtPRC/SpectraC5.fit.root\0");
      break;
    case 2:
      setOutputFile("data/HbtPRC/SpectraC5.fit.root\0");
      break;
    case 3:
      setOutputFile("data/HbtPRC/SpectraC5.fit.root\0");
      break;
    }
    break;
  }

  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/HbtPRC/STARSpectraPRL.root");
  }
  if((mV2 || mSpectraFit==2) & !mV2InputFile){
    mV2InputFile = new TFile("data/HbtPRC/STARV2PRC.root");
  }
  // Set centralities
  mNCent = mAsHbt? 5 : 6;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 351.5; // 0-5 C0
  mNParticipant[1] = 299.0;  // 5-10 C1
  mNParticipant[2] = 234.5;  // 10-20 C2
  mNParticipant[3] = 167.2; // 20-30 C3
  if(mAsHbt){
    mNParticipant[4] = 56.7;  // 30-80 C4-5-7-8
  }
  else{
    mNParticipant[4] = 96.5;  // 30-50 
    mNParticipant[5] = 30.2;  // 50-80 
  }

  // Set particle infos
  mNParticles = 9;
  mParticleName = new char*[mNParticles];
  char tParticleName[9][6] = {"Pip","Pim","Kp","Km","Pp","Pm","Pi","K","P"};
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
  mParticleMass[1] = MPi;
  mParticleMass[2] = MK;
  mParticleMass[3] = MK;
  mParticleMass[4] = MP;
  mParticleMass[5] = MP;
  mParticleMass[6] = MPi;
  mParticleMass[7] = MK;
  mParticleMass[8] = MP;
}

void HbtPRCDataHandler::loadSpectra(int iPart){
  int tCentMergedAsHbt[] = {0,1,2,3,4,9};
  int tCentMergedHbt[]   = {0,1,2,3,4,6,9};
  int* tCentMerged = mAsHbt? tCentMergedAsHbt : tCentMergedHbt;
  if(tCentMerged[mCent] == tCentMerged[mCent+1]-1){ // no merging
    char tGraphName[100];
    sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[iPart],mCent);
    mSpectras[iPart] = 
      mBWFitter->addSpectra((TGraphErrors*) mSpectraInputFile->Get(tGraphName),
			    mParticleMass[iPart],0.5,1.2);
  }
  else{
    int tCent = tCentMerged[mCent];
    int tN(0);
    double* tX(0);
    double* tEX(0);
    double* tY(0);
    double* tEY(0);
    while(tCent<tCentMerged[mCent+1]){
      char tGraphName[100];
      sprintf(tGraphName,"STARdEdx%sAuAu%i",mParticleName[iPart],tCent);
      TGraphErrors* tG = (TGraphErrors*) mSpectraInputFile->Get(tGraphName);
      if(!tN){
	tN = tG->GetN();
	tX = tG->GetX();
	tEX = tG->GetEX();
	tY = tG->GetY();
	tEY = tG->GetEY();
	for(int ti=0; ti<tN; ti++){
	  tEY[ti] = tEY[ti]*tEY[ti];
	}
      }
      else{
	double* ttY = tG->GetY();
	double* ttEY = tG->GetEY();
	for(int ti=0; ti<tN; ti++){
	  tY[ti]+= ttY[ti];
	  tEY[ti] += (ttEY[ti]*ttEY[ti]);
	}
      }
      tCent++;
    }
    int tNCent = tCentMerged[mCent+1]-tCentMerged[mCent];
    for(int ti=0; ti<tN; ti++){
      tY[ti] /= tNCent;
      tEY[ti] = sqrt(tEY[ti])/tNCent;
    }
    mSpectras[iPart] = 
      mBWFitter->addSpectra(new TGraphErrors(tN,tX,tY,tEX,tEY),
			    mParticleMass[iPart],0.5,1.2);   
  }
}

void HbtPRCDataHandler::load(){
  //  mCent=3;
  if(mSpectraFit!=0){
    for(int ti=0; ti<6; ti++){
      loadSpectra(ti);
      mValidParticle[ti]=0;
    }
  }
  if(mSpectraFit!=2){
    switch(mAsHbt){
    case 0: // integrated
      {
	char tNameo[100];
	sprintf(tNameo,"RoutC%i",mCent+1);
	//sprintf(tNameo,"pimROutC%i",mCent);
	char tNames[100];
	sprintf(tNames,"RsideC%i",mCent+1);
	//sprintf(tNames,"pimRSideC%i",mCent);
	char tNamel[100];
	sprintf(tNamel,"RlongC%i",mCent+1);
	//sprintf(tNamel,"pimRLongC%i",mCent);
	mBWFitter->addHbtRadii((TGraphErrors*) mInputFile->Get(tNameo),
			       (TGraphErrors*) mInputFile->Get(tNames),
			       (TGraphErrors*) mInputFile->Get(tNamel),
			       MPi);
      }
      break;
    case 1:
      {
	char tNameo[100];
	sprintf(tNameo,"Ro2F0C%i",mCent);
	char tNames[100];
	sprintf(tNames,"Rs2F0C%i",mCent);
	char tNamel[100];
	sprintf(tNamel,"Rl2F0C%i",mCent);
	mBWFitter->addHbtRadiiSquare((TGraphErrors*) mInputFile->Get(tNameo),
				     (TGraphErrors*) mInputFile->Get(tNames),
				     0,
				     (TGraphErrors*) mInputFile->Get(tNamel),
				     MPi);
      }
      break;
    case 2:
      {
	//	char tTestName[50];
	//sprintf(tTestName,"Ro2Phi0C%i",mCent);      
	//TGraphErrors* tTestG = (TGraphErrors*) mInputFile->Get(tTestName);
	//while(mCent<mNCent && !tTestG){
	//mCent++;
	//char ttTestName[50];
	//sprintf(ttTestName,"Ro2Phi0C%i",mCent);      
	//tTestG = (TGraphErrors*) mInputFile->Get(ttTestName);
	//}
	//if(mCent<mNCent){
	//fixParameters();
	int tNPhi=4;
	double tStepPhi = TMath::Pi()/4.;
	for(int tIPhi=0; tIPhi<tNPhi; tIPhi++){
	  char tNameO[50];
	  sprintf(tNameO,"Ro2Phi%iC%i",tIPhi,mCent);
	  char tNameS[50];
	  sprintf(tNameS,"Rs2Phi%iC%i",tIPhi,mCent);
	  char tNameOS[50];
	  sprintf(tNameOS,"Ros2Phi%iC%i",tIPhi,mCent);
	  char tNameL[50];
	  sprintf(tNameL,"Rl2Phi%iC%i",tIPhi,mCent);
	  mBWFitter->addHbtRadiiSquareWRTReactionPlane((TGraphErrors*) mInputFile->Get(tNameO),
						       (TGraphErrors*) mInputFile->Get(tNameS),
						       (TGraphErrors*) mInputFile->Get(tNameOS),
						       (TGraphErrors*) mInputFile->Get(tNameL),
						       MPi,tIPhi*tStepPhi);
	  }
	  //}
	  //else{
	  //nextFit();
	  //}
      }
      break;
    case 3:
      {
	char tNameo[100];
	sprintf(tNameo,"Ro2F0C%i",mCent);
	char tNames[100];
	sprintf(tNames,"Rs2F0C%i",mCent);
	char tNamel[100];
	sprintf(tNamel,"Rl2F0C%i",mCent);
	char tNameo2[100];
	sprintf(tNameo2,"Ro2F2C%i",mCent);
	char tNames2[100];
	sprintf(tNames2,"Rs2F2C%i",mCent);
	char tNameos2[100];
	sprintf(tNameos2,"Ros2F2C%i",mCent);
	char tNamel2[100];
	sprintf(tNamel2,"Rl2F2C%i",mCent);
	mBWFitter->addHbtRadiiSquare((TGraphErrors*) mInputFile->Get(tNameo),
				     (TGraphErrors*) mInputFile->Get(tNames),
				     0,
				     (TGraphErrors*) mInputFile->Get(tNamel),
				     MPi);
	mBWFitter->addHbtRadiiSquareOsc((TGraphErrors*) mInputFile->Get(tNameo2),
					(TGraphErrors*) mInputFile->Get(tNames2),
					(TGraphErrors*) mInputFile->Get(tNameos2),
					(TGraphErrors*) mInputFile->Get(tNamel2),
					MPi);
      }
      break;
    case 4:
      {
	int tNPhi=2;
	double tStepPhi = TMath::Pi()/2.;
	for(int tIPhi=0; tIPhi<tNPhi; tIPhi++){
	  char tNameO[50];
	  sprintf(tNameO,"Ro2RcPhi%iC%i",tIPhi,mCent);
	  char tNameS[50];
	  sprintf(tNameS,"Rs2RcPhi%iC%i",tIPhi,mCent);
	  char tNameL[50];
	  sprintf(tNameL,"Rl2RcPhi%iC%i",tIPhi,mCent);
	  mBWFitter->addHbtRadiiSquareWRTReactionPlane((TGraphErrors*) mInputFile->Get(tNameO),
						       (TGraphErrors*) mInputFile->Get(tNameS),
						       0,
						       (TGraphErrors*) mInputFile->Get(tNameL),
						       MPi,tIPhi*tStepPhi);
	}
	char tNameOS[50];
	sprintf(tNameOS,"Ros2RcPhi3C%i",mCent);
	mBWFitter->addHbtRadiiSquareWRTReactionPlane(0,0,(TGraphErrors*) mInputFile->Get(tNameOS),0,MPi,TMath::Pi()/4.);
						     
      }
      break;
    }
  }
  if(mV2 || mSpectraFit==2){
    double tLowLimPt[3] = {0.4,0.,0.3};
    double tHighLimPt[3] = {0.85, 0.5, 1.1};
    for(int ti=6; ti<mNParticles; ti++){
      mValidParticle[ti]=1;
      char tGraphName[100];
      sprintf(tGraphName,"%sC%i",mParticleName[ti],mCent);
      mBWFitter->addV2((TGraphErrors*) mV2InputFile->Get(tGraphName),
		       mParticleMass[ti],
		       tLowLimPt[ti-6], tHighLimPt[ti-6]);
    }
  }
  else{
    mValidParticle[6]=1;
  }
  
}
			  

void HbtPRCDataHandler::close(){
  //STAR200asHBTSpectraHandler::close();
  IOHandler::close();
  mSpectraInputFile->Close();
  if(mV2InputFile) mV2InputFile->Close();
}

void HbtPRCDataHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  if(mSpectraFit==0){
    mBWFitter->setParameters(mFixedT[mCent], mFixedRho0[mCent], 0., 0.,
			     9.,9.,0.,7.,2.);
    mBWFitter->fixT();
    mBWFitter->fixRho0();
  }
  if(mSpectraFit==2){
    //mBWFitter->fixRho2();
    mBWFitter->fixRx();
    //mBWFitter->fixRy();
    mBWFitter->fixTau();
    mBWFitter->fixDeltat();
  }
  else{
    if(mAsHbt==0 ||  mAsHbt==1){
      mBWFitter->fixRho2();
      mBWFitter->useR();
    }
  }
  mBWFitter->fixRho4();
  mBWFitter->fixAs(); // need an extra switch
}


