#include "STAR200asHBTHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <cstdlib>
using namespace std;
STAR200asHBTHandler::STAR200asHBTHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :IOHandler(aBWFitter,aInputFileName)
{
  mSpectra =0;
  mHBT=0;
  mV2=1;
}

void STAR200asHBTHandler::initLoad(){
  // set Files
  setOutputFile("data/HBT/STAR200asHBT.fit.root\0");
  if(!mInputFile){
    mInputFile = new TFile("data/HBT/HbtRP200GeV.root");
  }
  // Set centralities
  mNCent = 5;
  mNParticipant = new double[mNCent];
  mNParticipant[0] = 300.; // 0-10
  mNParticipant[1] = 141.635; // 20-40
  mNParticipant[2] = 62.4953;  // 40-60
  mNParticipant[3] = 30.; // 20-40
  mNParticipant[4] = 10.;  // 40-60

   
  // Set particle infos
  mNParticles = 1;
  mParticleName = new char*[mNParticles];
  char tParticleName[1][6] = {"Pi"};
  mParticleMass = new double[mNParticles];
  mParticleMass[0] = MPi;
  for(int ti=0; ti<mNParticles; ti++){
    mParticleName[ti] = new char[10];
    strcpy(mParticleName[ti],tParticleName[ti]);
  }
  // Set blast wave fitter parameters
}

void STAR200asHBTHandler::load(){

  for(int ti=0; ti<mNParticles; ti++){
    
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
				 mParticleMass[ti]);
    mBWFitter->addHbtRadiiSquareOsc((TGraphErrors*) mInputFile->Get(tNameo2),
				    (TGraphErrors*) mInputFile->Get(tNames2),
				    (TGraphErrors*) mInputFile->Get(tNameos2),
				    (TGraphErrors*) mInputFile->Get(tNamel2),
				    mParticleMass[ti]);
  }
}

void STAR200asHBTHandler::fixParameters(){
  mBWFitter->freeAllParameters();
  mBWFitter->fixAs(); // need an extra switch
  mBWFitter->fixT();
  mBWFitter->fixRho0();
}
			  
