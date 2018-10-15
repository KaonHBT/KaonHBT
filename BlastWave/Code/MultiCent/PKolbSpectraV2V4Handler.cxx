#include "PKolbSpectraV2V4Handler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <iostream>
#include <cstdlib>
using namespace std;

PKolbSpectraV2V4Handler::PKolbSpectraV2V4Handler(BlastWaveFitter* aBWFitter,
				       int aV4,
				       const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName),mV4(aV4)
{
  mSpectra =1;
  mHBT=0;
  mV2=1;
}

void PKolbSpectraV2V4Handler::initLoad(){
  // set Files
  if(mV4){
    setOutputFile("data/PKolbSpectraV2V4.fit.root\0");
  }
  else{
    setOutputFile("data/PKolbSpectraV2.fit.root\0");
  }
  if(!mInputFile){
    mInputFile = new TFile("data/PKolbHydro.root");
  }
  // Set centralities
  mNCent = 1;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 351.459;
  // Set particle infos
  mNParticles = 1;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Pi","Pim","Kp","Km","Pp","Pm"};
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
}

void PKolbSpectraV2V4Handler::load(){
  mSpectras[0] = 
    mBWFitter->addSpectra((TGraphErrors*) mInputFile->Get("Spectra"),
			  MPi,0.2,1.2,0);
  mBWFitter->addV2((TGraphErrors*) mInputFile->Get("v2"),MPi,0.2,1.2,0);
  if(mV4) mBWFitter->addV4((TGraphErrors*) mInputFile->Get("v4"),MPi,0.2,1.2,0);
}


void PKolbSpectraV2V4Handler::fixParameters(){
  mBWFitter->freeAllParameters();
  if(!mV4) mBWFitter->fixRho4();
  mBWFitter->fixAs(); // need an extra switch
  mBWFitter->fixRy();
  mBWFitter->fixTau();
  mBWFitter->fixDeltat();
}
			  
