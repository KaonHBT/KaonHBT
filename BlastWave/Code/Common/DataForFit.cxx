#include "DataForFit.h"

#include "TH1.h"
#include "TGraphErrors.h"

#include <iostream>
#include <cstdlib>
using namespace std;

DataForFit::DataForFit(TGraphErrors* aGraph, 
		       double aMass, 
		       double aPtMin, 
		       double aPtMax, 
		       int aStat)
  : mCurrentPoint(0),mMass(aMass),mStart(1),mSpectra(0),mStat(aStat)
{
  if(!aGraph) {
    cerr << "!!!!!!!!!!!!!!!!!!!!!! " << endl 
	 << "No input graph " << endl << "!!!!!!!!!!!!!!!!!!!!!! " << endl;
    mNPoints =0;
  }
  else{
    strcpy(mName,aGraph->GetName());
    int tN = aGraph->GetN();
    mX = new double[tN];
    mY = new double[tN];
    mYErrSquare = new double[tN];
    mNPoints=0;
    double tX;
    double tYErr;
    for(int ti=0; ti< tN; ti++){
      tX = aGraph->GetX()[ti];
      tYErr = aGraph->GetEY()[ti];
      if((tX>aPtMin) && (tX<aPtMax) && (tYErr>0.)){       
	mX[mNPoints] = tX;
	mY[mNPoints] = aGraph->GetY()[ti];	
	mYErrSquare[mNPoints] = tYErr*tYErr;
	mNPoints++;
      }
    }
  }
  cout << "Add " << aGraph->GetName() << " " << mNPoints
       << " " << mX[0] << " < pt < " << mX[mNPoints-1] << " " << endl;
//cout<<"check1"<<endl;
}

DataForFit::DataForFit(TH1* aHisto, 
		       double aMass, 
		       double aPtMin, 
		       double aPtMax,
		       int aStat)
  : mCurrentPoint(0),mMass(aMass),mStart(1),mSpectra(0),mStat(aStat)
{
  if(!aHisto) {
    cerr << "!!!!!!!!!!!!!!!!!!!!!! " << endl 
	 << "No input graph " << endl << "!!!!!!!!!!!!!!!!!!!!!! " << endl;
    mNPoints =0;
  }
  else{
    strcpy(mName,aHisto->GetName());
    cout << "Add " << aHisto->GetName() << endl;
    int tN = aHisto->GetNbinsX();
    mX = new double[tN];
    mY = new double[tN];
    mYErrSquare = new double[tN];
    mNPoints=0;
    double tX;
    double tYErr;
    for(int ti=1; ti<= tN; ti++){      
      tX = aHisto->GetBinCenter(ti);
      tYErr = aHisto->GetBinError(ti);
      if(tX> aPtMin && tX < aPtMax && tYErr>0.){       
	mX[mNPoints] = tX;
	mY[mNPoints] = aHisto->GetBinContent(ti);
	mYErrSquare[mNPoints] = tYErr*tYErr;
	mNPoints++;
      }
    }
  }
}


DataForFit::~DataForFit(){
  delete[] mX;
  delete[] mY;
  delete[] mYErrSquare;
}

void DataForFit::addDataForFit(TGraphErrors* aGraph){
  cerr << "ERROR! DataForFit::AddDataForFit Not implemented yet" << endl;
}

TH1D* DataForFit::histo(const char* aName, int NBin, BlastWave* aBW){
  TH1D* tH = new TH1D(aName, aName, NBin, mX[0], mX[mNPoints-1]);
  mStart=1;
  aBW->setStat(mStat);
  for(int ti=1; ti<=tH->GetNbinsX(); ti++){
    tH->SetBinContent(ti,getBlastWaveY(aBW, tH->GetBinCenter(ti)));
  }
  return tH;
}

TH1D* DataForFit::histo(const char* aName, int NBin, double XMin, double XMax,
			BlastWave* aBW){
  TH1D* tH = new TH1D(aName, aName, NBin, XMin, XMax);
  mStart=1;
  aBW->setStat(mStat);
  for(int ti=1; ti<=tH->GetNbinsX(); ti++){
    tH->SetBinContent(ti,getBlastWaveY(aBW, tH->GetBinCenter(ti)));
  }
  return tH;
}
 
SpectraForFit::SpectraForFit(TH1* aHisto, double aMass, 
			     double aPtMin, double aPtMax, int aStat)
    : DataForFit(aHisto,aMass, aPtMin, aPtMax, aStat){
  for(int ti=0; ti<mNPoints; ti++){
    mY[ti] /= mX[ti];
    mYErrSquare[ti] /= (mX[ti]*mX[ti]);
  }
  mSpectra=1;
}



//TH1D* SpectraForFit::histo(const char* aName, int NBin, BlastWave* aBW){
//TH1D* tH = DataForFit::histo(aName, NBin, aBW);
//tH->Scale(mNormData/mNormBW);
//return tH;
//}
