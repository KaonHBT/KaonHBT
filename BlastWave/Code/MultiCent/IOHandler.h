#ifndef IOHandler_h
#define IOHandler_h

class BlastWaveFitter;
class TFile;
class SpectraForFit;


class IOHandler{
 public:
  IOHandler(BlastWaveFitter* aBWFitter, 
	    const char* aInputFileName=0);
  virtual void setOutputFileName(const char* aFileName);
  virtual void setOutputFile(TFile* aFile);

  virtual ~IOHandler();

  void init();
  void nextCent();
  int nextFit();
  virtual void storeParameters(int ScaleErrorByChi2PerDof=0);
  virtual void writeSpectra();
  virtual void writeHBT();
  virtual void writeEllipticFlow();
  virtual void writeEStruct();
  virtual void writeParameters();
  virtual void writeContours();
  virtual void close();

  virtual void freeAs() {mFixedAs = 0;}

  //void cdToOutputFile(){mOutputFile->cd();}

 protected:
  virtual void fixParameters();
  virtual void initLoad()=0;
  virtual void load()=0;

  BlastWaveFitter* mBWFitter;

  int mNCent;
  int mCent;
  char mFileName[500];
  TFile* mInputFile;
  TFile* mOutputFile;
  void setOutputFile(const char* aFileName);
  
  int mHBT;
  int mV2;
  int mSpectra;
  int mFixedAs;
  int mNParticles;
  char** mParticleName; 
  double* mParticleMass;
  int* mValidParticle;
  SpectraForFit** mSpectras;
  double** mYield;

  // Parameter graph
  double* mNParticipant;
  double* mChi2;
  double* mConfLevel;
  double* mT;
  double* mTErr;
  double* mRho0;
  double* mRho0Err;
  double* mRho2;
  double* mRho2Err;
  double* mRho4;
  double* mRho4Err;
  double* mRx;
  double* mRxErr;
  double* mRy;
  double* mRyErr;
  double* mAs;
  double* mAsErr;
  double* mTau;
  double* mTauErr;
  double* mDeltat;
  double* mDeltatErr;
  void Rho0ToMeanBetaT(double aRho0, double aRho0Err, 
		       double& aBetaT, double& aBetaTErr);

  static double MPi;
  static double MK;
  static double MKStar;
  static double MP;
  static double MLa;
  static double MXi;
  static double MOm;
  static double MPhi;

};
#endif
