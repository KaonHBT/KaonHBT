#ifndef STAR200XiHandler_h
#define STAR200XiHandler_h
#include "IOHandler.h"

class STAR200XiHandler : public IOHandler{
 public:
  STAR200XiHandler(BlastWaveFitter* aBWFitter, 
		   const char* aInputFileName=0);
  ~STAR200XiHandler(){}
 protected:
  virtual void initLoad();
  virtual void load();
};

#endif
