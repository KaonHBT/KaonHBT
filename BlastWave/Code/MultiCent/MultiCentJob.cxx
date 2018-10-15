#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstdio>
using namespace std;

int main(int argc, char **argv){
  int aWhatKind = atoi(argv[1]);
  char RunDir[] = "/auto/rnc2/fretiere/Models/BlastWave"; 
  char JobName[100];
  sprintf(JobName,"BW%i",aWhatKind);
  char JobScript[200];
  sprintf(JobScript,"%s/Job/%s.script",RunDir,JobName);
  ofstream fJob(JobScript);
  fJob << "#!/bin/csh" << endl;
  fJob << "#BSUB -J " << JobName << endl; 
  fJob << "#BSUB -N" << endl;
  char JobLog[200];
  sprintf(JobLog,"%s/Job/%s.log",RunDir,JobName);
  char RmJobLog[200];
  sprintf(RmJobLog,"rm -f %s",JobLog);
  system(RmJobLog);
  fJob << "#BSUB -e " << JobLog << endl;
  fJob << "#BSUB -o " << JobLog << endl;
  fJob << "#BSUB -q medium" << endl;
  fJob << "setenv NODEBUG 1" << endl;
  fJob << "cd /auto/rnc2/fretiere/Models/BlastWave" << endl;
  fJob << "BIN/MultiCentFitBW.exe " << aWhatKind << " 3" << endl;
  fJob.close();
  char Submit[200];
  sprintf(Submit,"bsub < %s",JobScript);
  //cout << Submit << endl;
  system(Submit);
}
