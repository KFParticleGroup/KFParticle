
/*
  do it in the interpreter before running this macro
  gSystem->AddIncludePath("-I$KFPARTICLE_DIR/install/usr/local/include");
  gSystem->Load("$KFPARTICLE_DIR/install/usr/local/lib/libKFParticle.so");
*/

#include "ConfigConstants.h"

void RunKFUnitTestMacro()
{
#ifndef ALIPHYSICS
  gSystem->AddIncludePath(" -I$KFPARTICLE_DIR/install/usr/local/include");
  //gSystem->Load("$KFPARTICLE_DIR/install/usr/local/lib/libKFParticle.dylib");
  gSystem->Load("$KFPARTICLE_DIR/install/usr/local/lib/libKFParticle.so");
  gROOT->ProcessLine("KFParticle start_script_part3324151532");
#endif

  gROOT->ProcessLine(".x MakeKFParticleTrees.cxx");
  gROOT->ProcessLine(".x MakeKFUnitTestHistos.cxx");
}
