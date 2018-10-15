#include "STAR200V2SpectraHandler.h"
#include "STAR200V2Handler.h"
#include "STARdEdx200SpectraHandler.h"
#include "BlastWaveFitter.h"
#include "DataForFit.h"
#include "TFile.h"
#include <cstdlib>
using namespace std;
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
  setOutputFileName("data/STAR200V2Spectra.fit.root");
  mSpectraHandler->setOutputFile(mOutputFile);
  // Set blast wave fitter parameters
}

void STAR200V2SpectraHandler::load(){
  STAR200V2Handler::load();
  mSpectraHandler->nextCent();
  mSpectraHandler->load();
}
			  
void STAR200V2SpectraHandler::writeSpectra(){
  mSpectraHandler->writeSpectra();
}
