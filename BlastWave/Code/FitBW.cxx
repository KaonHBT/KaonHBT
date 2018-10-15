#include "BlastWaveFitter.h"

#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2.h"

#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;

int main(int argc, char **argv){
  double MPi = 0.1396; 
  double MK = 0.4937;
  double MKStar = 0.892;
  double MP = 0.9383;
  double MLa = 1.1157;
  double MXi = 1.32132;
  double MOm = 1.67245;
  double  MPhi = 1.019413;
  if(argc!=3){
    cout << "Need 2 arguments" << endl;
    return 1;
  }
  int iDataType = atoi(argv[1]); // 0: central130, 1: mid-perif130, 2: perif130, 3: cent200, 4: mid200, 5: perif 200
  int iMinimize = atoi(argv[2]); // 0: none, 1: migrad, 2: migrad+minos, 3: migrad spectra, 4: 

  TROOT ("myROOT","ttt");
  // _____________________________________________________________________
  // --- Instantiate Blast wave
  BlastWaveFitter* BWFitter =  BlastWaveFitter::instance();
  SpectraForFit* SpecPim =0;
  SpectraForFit* SpecKm =0;
  SpectraForFit* SpecP =0;
  SpectraForFit* SpecLa =0;
  SpectraForFit* SpecPip =0;
  SpectraForFit* SpecKp =0;
  SpectraForFit* SpecPB =0;
  SpectraForFit* SpecLaB =0;
  SpectraForFit* SpecK0s =0;
  SpectraForFit* SpecXi =0;
  SpectraForFit* SpecXiB =0;
  SpectraForFit* SpecOm =0;
   SpectraForFit* SpecPhi =0; 
   SpectraForFit* SpecKStar =0; 
  // _____________________________________________________________________
  // --- Get data and plug them in
  TFile fIn("data/DataForFullBW.root");
  
  double T;      // Temperature GeV
  double Rho0;    // Flow rapidity
  double Rho2;    // Flow modulation (phi dependence)
  double Rho4(0);
  double Rx;      // in-plane radius (fm) 
  double Ry;      // out-of-plane radius (fm)
  double As;      // Smoothness of the source, 0= box -> 0.3 = smooth
  double Tau;      // System life time (fm/c)
  double Deltat; // Emission duration
  char fOutName[500];

  if(iDataType<100){
    switch(iDataType){
    case 0:
      SpecPim  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent0to5"), 
			     MPi, 0.4, 2.);
      SpecPip  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent0to5"), 
			     MPi, 0.4, 2.);
      SpecKm   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent0to5"), 
			     MK, 0.3, 1.5);
      SpecKp   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent0to5"), 
			     MK, 0.3, 1.5);
      SpecPB = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pbarCent0to5"), 
			     MP, 0.3, 2.5);
      SpecP    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent0to5"), 
			     MP, 0.3, 2.5);
      SpecLa    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent0to5"), 
			     MLa, 0.3, 2.5);
      SpecLaB    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent0to5"), 
			     MLa, 0.3, 2.5);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPimC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RSPimC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RLPimC0to12"), 
			    MPi);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPipC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RSPipC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RLPipC0to12"), 
			    MPi);
      /*BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Phenix130ROPim"),
	(TGraphErrors*) fIn.Get("Phenix130RSPim"),
	(TGraphErrors*) fIn.Get("Phenix130RLPim"), MPi);
	BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Phenix130ROPip"),
			  (TGraphErrors*) fIn.Get("Phenix130RSPip"),
			  (TGraphErrors*) fIn.Get("Phenix130RLPip"), MPi);*/
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PiCent0to11"), MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PCent0to11"), MP);
      T = 0.106;      
      Rho0 = 0.89;   
      Rho2 = 0.060;   
      Rx = 13.2;     
      Ry = 13.0;    
      As = 0.;
      Tau = 9.2;    
      Deltat = 0.003; 
      strcpy(fOutName,"cBWAuAu130Cent.root" );
      break;
    case 1:
      SpecPim  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent15to30"), 
			     MPi, 0.4, 1.2);
      SpecPip  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent15to30"), 
			     MPi, 0.4, 1.2);
      SpecKm   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent15to30"), 
			     MK, 0.3, 1.5);
      SpecKp   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent15to30"), 
			     MK, 0.3, 1.5);
      SpecPB = 
	BWFitter->addSpectra((TGraphErrors*)fIn.Get("Phenix130pbarCent15to30"), 
			     MP, 0.3, 2.);
      SpecP = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent15to30"), 
			     MP, 0.3, 2.);
      SpecLa = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent20to35"), 
			     MLa, 0.3, 2.);
      SpecLaB = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent20to35"), 
			     MLa, 0.3, 2.);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPimC12to32"),
			    (TGraphErrors*) fIn.Get("Star130RSPimC12to32"),
			    (TGraphErrors*) fIn.Get("Star130RLPimC12to32"), 
			    MPi);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPipC12to32"),
			    (TGraphErrors*) fIn.Get("Star130RSPipC12to32"),
			    (TGraphErrors*) fIn.Get("Star130RLPipC12to32"), 
			    MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PiCent11to45"), MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PCent11to45"), MP);
      T = 0.107;// 0.111;
      Rho0 = 0.85; //0.818;
      Rho2 = 0.050; //0.058;
      Rx = 10.4;//8.3;
      Ry = 11.8;//9.6;
      As = 0.;
      Tau = 7.7;//5.9;
      Deltat = 0.06;//0.0002;
      strcpy(fOutName,"cBWAuAu130MidPerif.root" );
      break;
    case 2:
      SpecPim  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent60to92"), 
			     MPi, 0.4, 1.2);
      SpecPip  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent60to92"), 
			     MPi, 0.4, 1.2);
      SpecKm   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent60to92"), 
			     MK, 0.3, 1.5);
      SpecKp   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent60to92"), 
			     MK, 0.3, 1.5);
      SpecPB =
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pbarCent60to92"), 
			     MP, 0.3, 2.);
      SpecP = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent60to92"), 
			     MP, 0.3, 2.);
      SpecLa = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent35to75"), 
			     MLa, 0.3, 2.);
      SpecLaB = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent35to75"), 
			     MLa, 0.3, 2.);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPimC32to72"),
			    (TGraphErrors*) fIn.Get("Star130RSPimC32to72"),
			    (TGraphErrors*) fIn.Get("Star130RLPimC32to72"), 
			    MPi);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPipC32to72"),
			    (TGraphErrors*) fIn.Get("Star130RSPipC32to72"),
			    (TGraphErrors*) fIn.Get("Star130RLPipC32to72"), 
			    MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PiCent45to85"), MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("V2PCent45to85"), MP);
      T = 0.1;// 0.111;
      Rho0 = 0.79; //0.818;
      Rho2 = 0.05; //0.058;
      Rx = 8.0;//8.3;
      Ry = 10.1;//9.6;      
      As = 0.;
      Tau = 6.5;//5.9;
      Deltat = 0.6;//0.0002;
      strcpy(fOutName,"cBWAuAu130Perif.root" );
      break;
    case 3: // preliminary central AuAu 200
      SpecK0s = BWFitter->addSpectra((TGraphErrors*) fIn.Get("PaulsSpecK0s0to5")
				     , MK, 0.5, 2.);
      SpecLa = BWFitter->addSpectra((TGraphErrors*) fIn.Get("HuisSpecLa0to5"),
				    MLa, 0., 2.5);
      SpecPim = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC0"), 
				     MPi);
      SpecKm = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC0"), MK);
      SpecP = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC0"), MP);
      SpecPB = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC0"), MP);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("MercRO0to10"),
			    (TGraphErrors*) fIn.Get("MercRS0to10"),
			    (TGraphErrors*) fIn.Get("MercRL0to10"),
			    MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2K0s0to5"), MK, 0.5, 2.);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2LaB0to5"), 
		      MLa, 0.5, 2.5);
      T = 0.0916;      
      Rho0 = 0.993;   
      Rho2 = 0.0161;   
      Rx = 13.86;     
      Ry = 14.25;   
      As = 0.;   
      Tau = 10.34;    
      Deltat = 1.89; 
      strcpy(fOutName,"BWAuAu200Cent.root" );
      break;
    case 4: // preliminary MidPerif AuAu 200  
      SpecPim = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC2"), MPi);
      SpecKm = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC2"), MK);
      SpecP = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC2"), MP);
      SpecPB = BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC2"), MP);  
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("MercRO10to30"), 
			    (TGraphErrors*) fIn.Get("MercRS10to30"), 
			    (TGraphErrors*) fIn.Get("MercRL10to30"), MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2K0s5to30"), 
		      MK, 0.5, 2.);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2LaB5to30"), 
		      MLa, 0.5, 2.5);
      T = 0.092;// 0.111;
      Rho0 = 0.870; //0.818;
      Rho2 = 0.038; //0.058;
      Rx = 8.67;//8.3;
      Ry = 10.5;//9.6;  
      As = 0.;
      Tau = 7.4;//5.9;
      Deltat = 1.5;//0.0002;
      strcpy(fOutName,"BWAuAu200MidPerif.root" );
      break;
    case 5: // preliminary Perif AuAu 200    
        SpecPim = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisPimC6"), MPi);
      SpecKm = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaisKmC6"), MK);
      SpecP = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispC6"), MP);
      SpecPB = 
      BWFitter->addSpectra((TGraphErrors*) fIn.Get("KaispBC6"), MP);
      SpecK0s = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("PaulsSpecK0s40to60"), 
			     MK, 0.5, 2.);
      SpecLa = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("HuisSpecLa40to60"), 
			     MLa, 0., 2.5);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("MercRO30to70"), 
			    (TGraphErrors*) fIn.Get("MercRS30to70"), 
			    (TGraphErrors*) fIn.Get("MercRL30to70"), MPi);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2K0s30to70"), 
		      MK, 0.5, 2.);
      BWFitter->addV2((TGraphErrors*) fIn.Get("PaulsV2LaB30to70"), 
		      MLa, 0.5, 2.5);
      T = 0.092;// 0.111;
      Rho0 = 0.870; //0.818;
      Rho2 = 0.038; //0.058;
      Rx = 8.67;//8.3;
      Ry = 10.5;//9.6;  
      As = 0.;
      Tau = 7.4;//5.9;
      Deltat = 1.5;//0.0002;
      strcpy(fOutName,"BWAuAu200Perif.root" );
      break;
    case 6:
      iMinimize = 3;
      SpecPim = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasPim"), MPi,0.3,2.);
      SpecKm = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasKm"), MK);
      SpecKp = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasKp"), MK);
      //SpecP = 
        // BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasP"), MP);
      //SpecPB = 
        // BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasPB"), MP);
      SpecLa = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasLa"), MLa);
      SpecLaB = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("MasLaB"), MLa);         
      T = 0.1;
      Rho0 = 0.9;
      Rho2 = 0.;
      Rx = 10.;
      Ry = 10.;  
      As = 0.;
      Tau = 10.;
      Deltat = 1.;
      strcpy(fOutName,"Masashi130.root" );    
      break;    
    case 7:
      iMinimize = 3;
      SpecXi = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Xi130C0to10"), MXi);      
      SpecXiB = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("AXi130C0to10"), MXi);                   
      SpecOm = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Om130C0to10"), MOm);      
      SpecPhi = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phi130Cent0to11"), MPhi);              
            T = 0.1;
      Rho0 = 0.9;
      Rho2 = 0.;
      Rx = 10.;
      Ry = 10.;  
      As = 0.;
      Tau = 10.;
      Deltat = 1.;
      strcpy(fOutName,"MultiStrange130.root" );    
      break;    
   case 8:  
      iMinimize = 3;
      SpecPhi = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phi130Cent0to11"), MPhi);              
            T = 0.1;
      Rho0 = 0.9;
      Rho2 = 0.;
      Rx = 10.;
      Ry = 10.;
      Tau = 10.;
      Deltat = 1.;
      strcpy(fOutName,"Phi130.root" );    
      break;          
    case 9: // Test for Mike		         
      T = 0.1;
      Rho0 = 0.9;
      Rho2 = 0.05;
      Rx = 16.5;
      Ry = 19.5;  
      As = 0.;
      Tau = 9.;
      Deltat = 2.;
      //strcpy(fOutName,"TestForMike.root" );    
      BWFitter->setParameters(T,Rho0,Rho2,Rho4,Rx,Ry,As,Tau,Deltat);
      BWFitter->blastWave()->testForMike("HomeMade/BlastWave/pions-vs-phi.BW");
      return 0;
      break;        
      case 10: // Hbt only with T and rho fix to value from Xi fir
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPimC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RSPimC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RLPimC0to12"), 
			    MPi);
      BWFitter->addHbtRadii((TGraphErrors*) fIn.Get("Star130ROPipC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RSPipC0to12"),
			    (TGraphErrors*) fIn.Get("Star130RLPipC0to12"), 
			    MPi);
		T = 0.155;
		Rho0 = 0.72;
		Rho2 = 0.;
		Rx = 13.;
		Ry = 13.;  
      As = 0.;
		Tau = 9.;
		Deltat = 2.;
		strcpy(fOutName,"HbtHighTemp.root" );
		break;			            
    }
  }
  else{
    int iCent = iDataType- ((int) floor(iDataType/100.))*100;
    cout << iCent << endl;
    iMinimize = 3;
    char tPimName[50];
    sprintf(tPimName,"KaisPimC",iCent);
    SpecPim = 
      BWFitter->addSpectra((TGraphErrors*) fIn.Get(tPimName), MPi);
    char tKmName[50];
    sprintf(tKmName,"KaisKmC",iCent );
    SpecKm = 
      BWFitter->addSpectra((TGraphErrors*) fIn.Get(tKmName), MK);
    char tpName[50];
    sprintf(tpName,"KaispC",iCent );
    SpecP = 
      BWFitter->addSpectra((TGraphErrors*) fIn.Get(tpName), MP);
    char tpBName[50];
    sprintf(tpBName,"KaispBC",iCent );
    SpecPB = 
    BWFitter->addSpectra((TGraphErrors*) fIn.Get(tpBName), MP);
    T = 0.1;
    Rho0 = 0.9;
    Rho2 = 0.;
    Rx = 10.;
    Ry = 10.;
    Tau = 10.;
    Deltat = 1.;
    sprintf(fOutName,"BWAuAu200C%i.root",iCent);
  }

  //fIn.Close();




  // _____________________________________________________________________
  // --- Fitting



  BWFitter->setParameters(T,Rho0,Rho2,Rho4,Rx,Ry,As,Tau,Deltat);
  //BWFitter->fixAs();
  //BWFitter->fixT();
  //BWFitter->fixRho0();
  //BWFitter->fixRho2();
  //BWFitter->useR(); // merge Rx and Ry
  //BWFitter->fixRx();
  //BWFitter->fixRy();
  //BWFitter->fixTau();
  //BWFitter->fixDeltat();
  cout << "Ready to minimize. NDOF = " << BWFitter->nDOF() << endl;
  switch(iMinimize){
  case 0: // nothing
    break;
  case 1: // minimize
    BWFitter->minimize();
    break;
  case 2: // minimize and minos
    BWFitter->minimize();
    BWFitter->calculatePreciseErrors();
    break;
  case 3: // mininize only spectra
    BWFitter->fixRho2();
    BWFitter->fixRx();
    BWFitter->fixRy();
    BWFitter->fixTau();
    BWFitter->fixDeltat();
    BWFitter->minimize();
    break;
   case 4: // minimize only time
   BWFitter->fixT();
   BWFitter->fixRho0();
    BWFitter->fixRho2();
    BWFitter->fixRx();
    BWFitter->fixRy();   
    BWFitter->minimize();    
    BWFitter->calculatePreciseErrors();    
    break;
    case 5: // minimize only Hbt
   BWFitter->fixT();
   BWFitter->fixRho0();
    BWFitter->fixRho2();    
	BWFitter->useR();
    BWFitter->minimize();        
    break;
  }
  BWFitter->coutFitResults();
    
  // _____________________________________________________________________
  // --- Output 
  TFile* fOut = new TFile(fOutName,"RECREATE");  
  if(iMinimize!=0 && iMinimize!=5){
	   BWFitter->betaTContour(1)->Write(); 
	   BWFitter->betaTContour(2)->Write();   
	   BWFitter->betaTContour(3)->Write(); 
	   //BWFitter->confLevelMap("BetaT",1,19,0.4,1.2,0,19,0.05,0.25)->Write();
	   if(iMinimize==1 || iMinimize==2){
		   BWFitter->S2VsRho2Contour(1)->Write();
		   BWFitter->S2VsRho2Contour(2)->Write();
		   BWFitter->S2VsRho2Contour(3)->Write();
		   //BWFitter->confLevelMap("TauDt",5,19,7.,11.,6,9,0.,1.)->Write();
		   //BWFitter->TauVsDeltatContour(1)->Write();
		   //BWFitter->TauVsDeltatContour(2)->Write();	  
		   //BWFitter->TauVsDeltatContour(3)->Write();	  
	  }
	  
  }
  // _____________________________________________________________________
  // --- Get data againwithout including them in fit
  //fIn.Open("DataForFullBW.root");  
    if(iDataType<100){
    switch(iDataType){
	case 0:
	  	SpecKStar = 
	  		BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130KStarC0to14"), MKStar);
      SpecXi = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Xi130C0to10"), MXi);      
      SpecXiB = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("AXi130C0to10"), MXi);                   
      SpecOm = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Om130C0to10"), MOm);      
      SpecPhi = 
         BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phi130Cent0to11"), MPhi);     	
		break;
    case 7:
	  	SpecKStar = 
	  		BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130KStarC0to14"), MKStar);    
      SpecPim  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PimCent0to5"), 
			     MPi, 0.4, 1.2);
      SpecPip  = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130PipCent0to5"), 
			     MPi, 0.4, 1.2);
      SpecKm   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KmCent0to5"), 
			     MK, 0.3, 1.5);
      SpecKp   = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130KpCent0to5"), 
			     MK, 0.3, 1.5);
      SpecPB = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pbarCent0to5"), 
			     MP, 0.3, 2.);
      SpecP    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Phenix130pCent0to5"), 
			     MP, 0.3, 2.);
      SpecLa    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaCent0to5"), 
			     MLa, 0.3, 2.);
      SpecLaB    = 
	BWFitter->addSpectra((TGraphErrors*) fIn.Get("Star130LaBCent0to5"), 
			     MLa, 0.3, 2.);
  }
  
  }
  fIn.Close();
  

	fOut->cd();
  // Do it on the spectra to get the scaling (dN/dY) right
  if(SpecPim) 
    SpecPim->histo("BestPimPt", 100, BWFitter->blastWave())->Write();
  if(SpecPip) 
    SpecPip->histo("BestPipPt", 100, BWFitter->blastWave())->Write();
  if(SpecKm) SpecKm->histo("BestKmPt", 100, BWFitter->blastWave())->Write();
  if(SpecKp) SpecKp->histo("BestKpPt", 100, BWFitter->blastWave())->Write();
  if(SpecK0s) 
    SpecK0s->histo("BestK0sPt", 100, BWFitter->blastWave())->Write();
  if(SpecP) SpecP->histo("BestPPt", 100, BWFitter->blastWave())->Write();
  if(SpecPB) SpecPB->histo("BestPBPt", 100, BWFitter->blastWave())->Write();
  if(SpecLa) SpecLa->histo("BestLaPt", 100, BWFitter->blastWave())->Write();
  if(SpecLaB) 
    SpecLaB->histo("BestLaBPt", 100, BWFitter->blastWave())->Write();
  if(SpecXi) 
    SpecXi->histo("BestXiPt", 100, BWFitter->blastWave())->Write();
  if(SpecXiB) 
    SpecXiB->histo("BestXiBPt", 100, BWFitter->blastWave())->Write();
  if(SpecOm) 
    SpecOm->histo("BestOmPt", 100, BWFitter->blastWave())->Write();
  if(SpecPhi) 
    SpecPhi->histo("BestPhiPt", 100, BWFitter->blastWave())->Write();
  if(SpecKStar) 
    SpecKStar->histo("BestKStarPt", 100, BWFitter->blastWave())->Write();
    
    
  BWFitter->blastWave()->rOut("Pi",MPi)->Write();
  BWFitter->blastWave()->rSide("Pi",MPi)->Write();
  BWFitter->blastWave()->rLong("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongVsMt("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSideSquare("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rSideSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rLongSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->rOutSideSquareOsc("Pi",MPi)->Write();
  BWFitter->blastWave()->v2VsPt("Pi",MPi)->Write();
  BWFitter->blastWave()->v2VsPt("P",MP)->Write();
  BWFitter->blastWave()->v2VsPt("La",MLa)->Write();
  


  BWFitter->blastWave()->avSepOut("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOutSpace("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOutTime("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepSide("PiK",MPi,MK)->Write();  
  BWFitter->blastWave()->avSepLong("PiK",MPi,MK)->Write();
  BWFitter->blastWave()->avSepOut("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOutSpace("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOutTime("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepSide("PiP",MPi,MP)->Write();  
  BWFitter->blastWave()->avSepLong("PiP",MPi,MP)->Write();
  BWFitter->blastWave()->avSepOut("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepOutSpace("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepOutTime("KP",MK,MP)->Write();
  BWFitter->blastWave()->avSepSide("KP",MK,MP)->Write();  
  BWFitter->blastWave()->avSepLong("KP",MK,MP)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi200",MPi,0.2)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi200",MPi,0.2)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi400",MPi,0.4)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi400",MPi,0.4)->Write();

  BWFitter->blastWave()->rOutSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rSideSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rLongSquareVsPhiP("Pi600",MPi,0.6)->Write();
  BWFitter->blastWave()->rOutSideSquareVsPhiP("Pi600",MPi,0.6)->Write();

  BWFitter->blastWave()->v2VsPt("K",MK)->Write();
  BWFitter->blastWave()->rOut("K",MK)->Write();
  BWFitter->blastWave()->rSide("K",MK)->Write();
  BWFitter->blastWave()->rLong("K",MK)->Write();
  BWFitter->blastWave()->rOutVsMt("K",MK)->Write();
  BWFitter->blastWave()->rSideVsMt("K",MK)->Write();
  BWFitter->blastWave()->rLongVsMt("K",MK)->Write();



  fOut->Close();

  BWFitter->coutFitResults();
  
  delete BWFitter;


}
