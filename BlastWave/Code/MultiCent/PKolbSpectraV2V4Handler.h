#ifndef PKolbSpectraV2V4Handler_h
#define PKolbSpectraV2V4Handler_h
#include "IOHandler.h"
class TFile;

class PKolbSpectraV2V4Handler : public IOHandler{
 public:
  PKolbSpectraV2V4Handler(BlastWaveFitter* aBWFitter,
		     int aV4=1,
		     const char* aInputFileName=0);
  ~PKolbSpectraV2V4Handler(){}
 protected:
  TFile* mSpectraInputFile;
  virtual void fixParameters();
  virtual void initLoad();
  virtual void load();
  int mV4;
};

#endif
