#include "STAR200asHBTSpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
using namespace std;
STAR200asHBTSpectraHandler::STAR200asHBTSpectraHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =1;
  mHBT=1;
  mV2=1;
}

void STAR200asHBTSpectraHandler::initLoad(){
  // set Files
  setOutputFile("data/STAR200asHBTSpectra.fit.root\0");
  if(!mSpectraInputFile){
    mSpectraInputFile = new TFile("data/Spectra/STARdEdx200Spectra.root");
  }
  if(!mInputFile){
    mInputFile = new TFile("data/HBT/HbtRP200GeV.root");
  }
  // Set centralities
  mNCent = 5;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 300.; // 0-5 C0
  mNParticipant[1] = 141.635;  // 5-10 C1
  mNParticipant[2] = 62.4953;  // 10-20 C2
  mNParticipant[3] = 30.; // 20-30 C3
  mNParticipant[4] = 10.;  // 30-80 C4-5-7-8

  // Set particle infos
  mNParticles = 6;
  mParticleName = new char*[mNParticles];
  char tParticleName[10][6] = {"Pip","Pim","Kp","Km","Pp","Pm"};
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

  // Set blast wave fitter parameters
  mBWFitter->fixAs();

}


void STAR200asHBTSpectraHandler::loadSpectra(int iPart){
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


void STAR200asHBTSpectraHandler::load(){
  for(int ti=0; ti<mNParticles; ti++){
    loadSpectra(ti);
  }

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

  for(int ti=1; ti<mNParticles; ti++){
    mValidParticle[ti]=0;
  }

}
			  

void STAR200asHBTSpectraHandler::close(){
  IOHandler::close();
  mSpectraInputFile->Close();
}

void STAR200asHBTSpectraHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  mBWFitter->fixAs(); // need an extra switch
}
