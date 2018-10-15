#include "STAR200V2SpectraHandler.h"
#include "STAR200V2Handler.h"
#include "STARdEdx200SpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <stdlib.h>

STAR200V2SpectraHandler::STAR200V2SpectraHandler(BlastWaveFitter* aBWFitter, const char* aInputFileName)
  :STAR200V2Handler(aBWFitter,aInputFileName)
{
  mV2=1;
  mSpectra=1;
  mHBT=0;
}

void STAR200V2SpectraHandler::initLoad(){
  STAR200V2Handler::initLoad();
  mSpectraHandler = new STARdEdx200SpectraHandler(mBWFitter);
  mSpectraHandler->init();
  setOutputFileName("data/STAR200V2Spectra.root");
  mSpectraHandler->setOutputFile(mOutputFile);
  // Set blast wave fitter parameters
  mBWFitter->freeAllParameters();
  mBWFitter->fixAs();
  mBWFitter->fixRy();
  mBWFitter->fixTau();
  mBWFitter->fixDeltat();
}

void STAR200V2SpectraHandler::load(){
  STAR200V2Handler::load();
  mSpectraHandler->nextCent();
  mSpectraHandler->load();
}
			  
void STAR200V2SpectraHandler::writeSpectra(){
  mSpectraHandler->writeSpectra();
}
