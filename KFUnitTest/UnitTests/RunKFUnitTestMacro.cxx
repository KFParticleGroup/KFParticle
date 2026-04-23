#include "ConfigConstants.h"

void RunKFUnitTestMacro() {
    #ifndef ALIPHYSICS
    gSystem->AddIncludePath(" -I$KFPARTICLE_DIR/install/include");
    //gSystem->Load("$KFPARTICLE_DIR/install/lib/libKFParticle.dylib"); 
    gSystem->Load("$KFPARTICLE_DIR/install/lib/libKFParticle.so");
    gROOT->ProcessLine("KFParticle start_script_part3324151532");
    #endif

    gROOT->ProcessLine(".x MakeKFParticleTrees.cxx");
    gROOT->ProcessLine(".x MakeKFUnitTestHistos.cxx");
}
