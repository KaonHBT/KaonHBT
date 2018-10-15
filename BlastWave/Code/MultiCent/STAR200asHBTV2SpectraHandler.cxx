#include "STAR200asHBTV2SpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
using namespace std;
STAR200asHBTV2SpectraHandler::STAR200asHBTV2SpectraHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName, int aFourier, int aV2)
  //:STAR200asHBTSpectraHandler(aBWFitter,aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mV2=aV2;
  mFourier = aFourier; 
  mHBT=1;
  mSpectra=1;
}

void STAR200asHBTV2SpectraHandler::initLoad(){
  // set Files
  if(mV2){ 
    if(mFourier){
      setOutputFile("data/STAR200asHBTV2Spectra.fit.root\0");
    }
    else{
      setOutputFile("data/STAR200asHBTVsPhiV2Spectra.fit.root\0");
    }
  }
  else{
    if(mFourier){
      setOutputFile("data/STAR200asHBTSpectra.fit.root\0");
    }
    else{
      setOutputFile("data/STAR200asHBTVsPhiSpectra.fit.root\0");
    }
  }
  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/Spectra/STARdEdx200Spectra.root");
  }
  if(mV2 & !mV2InputFile){
    mV2InputFile = new TFile("data/v2/STAR200V2PRC.root");//STARPidAuAu200.root");
  }
  if(!mInputFile){
    if(mFourier){
      mInputFile = new TFile("data/HBT/HbtRP200GeV.root");
    }
    else{
      mInputFile = new TFile("data/HBT/STAR200asHBT.root");
    }
  }
  // Set centralities
  mNCent = 5;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 351.5; // 0-5 C0
  mNParticipant[1] = 299.0;  // 5-10 C1
  mNParticipant[2] = 234.5;  // 10-20 C2
  mNParticipant[3] = 167.2; // 20-30 C3
  mNParticipant[4] = 56.7;  // 30-80 C4-5-7-8

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

void STAR200asHBTV2SpectraHandler::loadSpectra(int iPart){
  int tCentMerged[] = {0,1,2,3,4,9};
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

void STAR200asHBTV2SpectraHandler::load(){
  for(int ti=0; ti<6; ti++){
    loadSpectra(ti);
    mValidParticle[ti]=0;
  }
  if(mFourier){
    if(mFourier==1){
    char tNameo[100];
    sprintf(tNameo,"Ro2_0C%i",mCent);
    char tNames[100];
    sprintf(tNames,"Rs2_0C%i",mCent);
    char tNamel[100];
    sprintf(tNamel,"Rl2_0C%i",mCent);
    char tNameo2[100];
    sprintf(tNameo2,"Ro2_2C%i",mCent);
    char tNames2[100];
    sprintf(tNames2,"Rs2_2C%i",mCent);
    char tNameos2[100];
    sprintf(tNameos2,"Ros2_2C%i",mCent);
    char tNamel2[100];
    sprintf(tNamel2,"Rl2_2C%i",mCent);
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
  }
  else{
    int tNPhi =4;
    double tStepPhi = TMath::Pi()/4.;
    for(int tIPhi=0; tIPhi<tNPhi; tIPhi++){
      char tNameO[50];
      sprintf(tNameO,"ROSqPhi%iC%i",tIPhi,mCent);
      char tNameS[50];
      sprintf(tNameS,"RSSqPhi%iC%i",tIPhi,mCent);
      char tNameOS[50];
      sprintf(tNameOS,"ROSSqPhi%iC%i",tIPhi,mCent);
      char tNameL[50];
      sprintf(tNameL,"RLSqPhi%iC%i",tIPhi,mCent);
      mBWFitter->addHbtRadiiSquareWRTReactionPlane(
	      (TGraphErrors*) mInputFile->Get(tNameO),
	      (TGraphErrors*) mInputFile->Get(tNameS),
    	      (TGraphErrors*) mInputFile->Get(tNameOS),
	      (TGraphErrors*) mInputFile->Get(tNameL),
	      MPi,tIPhi*tStepPhi);
    }
  }
  if(mV2){
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
			  

void STAR200asHBTV2SpectraHandler::close(){
  //STAR200asHBTSpectraHandler::close();
  IOHandler::close();
  mSpectraInputFile->Close();
  if(mV2InputFile) mV2InputFile->Close();
}

void STAR200asHBTV2SpectraHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  mBWFitter->fixRho4();
  mBWFitter->fixAs(); // need an extra switch
}


