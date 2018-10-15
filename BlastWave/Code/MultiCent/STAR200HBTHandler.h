#ifndef STAR200HBTHandler_h
#define STAR200HBTHandler_h
//#include "STAR200asHBTHandler.h"
#include "IOHandler.h"
class TFile;

class STAR200HBTHandler : public IOHandler{//STAR200asHBTHandler{
 public:
  STAR200HBTHandler(BlastWaveFitter* aBWFitter, 
			       const char* aInputFileName=0);
  ~STAR200HBTHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
};

#endif
