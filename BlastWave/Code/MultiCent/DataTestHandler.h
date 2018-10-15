#ifndef DataTestHandler_h
#define DataTestHandler_h
#include "IOHandler.h"
class TFile;

class DataTestHandler : public IOHandler{
 public:
  DataTestHandler(BlastWaveFitter* aBWFitter, 
		   const char* aInputFileName=0);
  ~DataTestHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  virtual void fixParameters();
  virtual void loadSpectra(int iPart, int iStat=0);
  virtual void close();
  TFile* mSpectraInputFile;
};

#endif
