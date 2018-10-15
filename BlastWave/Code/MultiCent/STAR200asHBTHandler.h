#ifndef STAR200asHBTHandler_h
#define STAR200asHBTHandler_h
#include "IOHandler.h"

class STAR200asHBTHandler : public IOHandler{
 public:
  STAR200asHBTHandler(BlastWaveFitter* aBWFitter, 
		   const char* aInputFileName=0);
  ~STAR200asHBTHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
  virtual void fixParameters();
};

#endif
