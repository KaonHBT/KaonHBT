#include "BlastWaveFitter.h"
#include "PHENIXAuAu200SpectraHandler.h"
#include "STARdEdx200SpectraHandler.h"
#include "STAR200V2Handler.h"
#include "STAR200V2V4Handler.h"
#include "STAR200V2SpectraHandler.h"
#include "STAR200XiHandler.h"
#include "PHENIX200SpectraV2Handler.h"
#include "PHENIX200SpectraV2HBTHandler.h"
#include "STAR200OmSpectraHandler.h"
#include "STAR200asHBTHandler.h"
#include "NA49SpectraHandler.h"
#include "STAR200asHBTSpectraHandler.h"
#include "STAR200asHBTV2SpectraHandler.h"
#include "STAR200HBTSpectraHandler.h"
#include "PHENIX200SpectraHBTHandler.h"
#include "STAR200HBTHandler.h"
#include "DataTestHandler.h"
#include "PKolbSpectraV2V4Handler.h"
#include "HbtPRCDataHandler.h"
#include "STARKplusKplusHBTSpectraHandler.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv){

  TROOT ("myROOT","ttt");
  int aWhatKind = argc>1 ? atoi(argv[1]) : 0; 
  int aWhatCalc = argc>2 ? atoi(argv[2]) : 0; // 0 = migrad 1 = migrad+minos 2 = migrad+contours, 3 = migrad+minos+contours

  // _____________________________________________________________________
  // --- Instantiate Blast wave
  BlastWaveFitter* BWFitter =  BlastWaveFitter::instance();
  double T;      // Temperature GeV
  double Rho0;    // Flow rapidity
  double Rho2;    // Flow modulation (phi dependence)
  double Rho4;
  double Rx;      // in-plane radius (fm) 
  double Ry;      // out-of-plane radius (fm)
  double As;      // Smoothness of the source, 0= box -> 0.3 = smooth
  double Tau;      // System life time (fm/c)
  double Deltat; // Emission duration
  BWFitter->setParameters(0.1,1.,0.,0.,13.,13.,0.,9.,2.);
  //BWFitter->setParameters(0.111,0.936,0.0255,0.0058,11.52,13.,0.,9.,2.); // speed up for Art
  //BWFitter->setParameters(0.098,1.0,0.0098,0.,12.5,13.,0.,9.,2.);
  //BWFitter->fixAs();

  // _____________________________________________________________________
  // --- Instantiate IO handler
  IOHandler* Handler;
  switch(aWhatKind){
  case 0:
    Handler = new PHENIXAuAu200SpectraHandler(BWFitter);
    break;
  case 1:
    Handler = new STARdEdx200SpectraHandler(BWFitter);
    break;
  case 2:
    Handler = new STAR200V2V4Handler(BWFitter,0);
    //Handler = new STAR200V2Handler(BWFitter);
    break;
  case 3:
    Handler = new STAR200V2V4Handler(BWFitter);
    break;
  case 4:
    Handler = new STAR200V2SpectraHandler(BWFitter);
    break;
  case 5:
    Handler = new STAR200XiHandler(BWFitter);
    break;
  case 6:
    Handler = new PHENIX200SpectraV2Handler(BWFitter);
    break;
  case 7:
    Handler = new STAR200OmSpectraHandler(BWFitter);
    break;
  case 8:
    Handler = new STAR200asHBTHandler(BWFitter);
    break;
  case 9:
    Handler = new NA49SpectraHandler(BWFitter);
    break;
  case 10:
    Handler = new STAR200asHBTV2SpectraHandler(BWFitter,0,1,0);
    break;
  case 11:
    Handler = new STAR200asHBTV2SpectraHandler(BWFitter,0,-1,1);
    break;
  case 12:
    Handler = new PHENIX200SpectraV2HBTHandler(BWFitter);
    break;
  case 13:
    Handler = new STAR200HBTSpectraHandler(BWFitter);
    break;
  case 14:
    Handler = new PHENIX200SpectraHBTHandler(BWFitter);
    break;
  case 15:
    Handler = new STAR200HBTHandler(BWFitter);
    break;
  case 16:
    Handler = new PKolbSpectraV2V4Handler(BWFitter,0);
    break;
  case 17:
    Handler = new PKolbSpectraV2V4Handler(BWFitter);
    break;
  case 20: // phi integrated HBT
    Handler = new HbtPRCDataHandler(BWFitter);
    break;
  case 21: // asHBT coef 0
    Handler = new HbtPRCDataHandler(BWFitter,0,1);
    break;
  case 22: // asHBT 
    Handler = new HbtPRCDataHandler(BWFitter,0,3,0);
    break;
  case 23: // asHBT + v2
    Handler = new HbtPRCDataHandler(BWFitter,0,3,1);
    break;
  case 24:
    Handler = new HbtPRCDataHandler(BWFitter,0,0,0,2);
    break;
  case 25:
    Handler = new HbtPRCDataHandler(BWFitter,0,3,0,2);
    break;
  case 26:
    Handler = new HbtPRCDataHandler(BWFitter,0,0,0,0);
    break;
  case 27:
    Handler = new HbtPRCDataHandler(BWFitter,0,3,1,0);
    break;
  case 28:
    Handler = new HbtPRCDataHandler(BWFitter,0,1,0,0);
    break;
  case 29:
    Handler = new HbtPRCDataHandler(BWFitter,0,3,0,0);
    break;
  case 30:
    Handler = new HbtPRCDataHandler(BWFitter,0,2,0,0);
    break;
  case 31:
    Handler = new HbtPRCDataHandler(BWFitter,0,2,1,0);
    break;
  case 32:
    Handler = new HbtPRCDataHandler(BWFitter,0,4,0,0);
    break;
  case 33:
    Handler = new HbtPRCDataHandler(BWFitter,0,4,1,0);
    break;
  case 50: // phi integrated HBT
    Handler = new HbtPRCDataHandler(BWFitter,0,0,0,1);
    break;
  case 51: // phi integrated HBT, independent fit to spectra 
    Handler = new HbtPRCDataHandler(BWFitter,0,0,0,0);
    break;
  case 52: // HBT vs phi 
    Handler = new HbtPRCDataHandler(BWFitter,0,2,0,1);
    break;
  case 53: // HBT vs phi, independent fit to spectra 
    Handler = new HbtPRCDataHandler(BWFitter,0,2,0,0);
    break;
  case 54: // HBT vs phi and v2
    Handler = new HbtPRCDataHandler(BWFitter,0,2,1,1);
    break;
  case 55: // HBT vs phi and v2, independent fit to spectra 
    Handler = new HbtPRCDataHandler(BWFitter,0,2,1,0);
    break;
  case 101:
    Handler = new STARdEdx200SpectraHandler(BWFitter);
    Handler->freeAs();
    break;
  case 113:
    Handler = new STAR200HBTSpectraHandler(BWFitter);
    Handler->freeAs();
    break;
  case 200:
    Handler = new STARKplusKplusHBTSpectraHandler(BWFitter);
    break;
 
  case 1000:
    Handler = new DataTestHandler(BWFitter);
    break;
  }
  Handler->init();

  // _____________________________________________________________________
  // --- Start fitting
  while(Handler->nextFit()){
    if(aWhatCalc>=0) BWFitter->minimize();
    if(aWhatCalc==1 || aWhatCalc ==3) BWFitter->calculatePreciseErrors();
    Handler->storeParameters();
    Handler->writeSpectra(); // do nothing if no spectra
    Handler->writeHBT(); 
    Handler->writeEllipticFlow();
    Handler->writeEStruct();
    if(aWhatCalc>1) Handler->writeContours();
    BWFitter->coutFitResults();
  }
  cout << "End fit" <<endl;
  // _____________________________________________________________________
  // --- End
  Handler->writeParameters();
  if(aWhatKind==6){ // calculate deuteron v2 (min bias)
    BWFitter->blastWave()->v2VsPt("d3",1.876)->Write();
  }
  Handler->close();

}
