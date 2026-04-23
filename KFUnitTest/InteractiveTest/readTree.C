/*
 .L $NOTE/JIRA/PWGPP-724/covarianceParam.C
 */
//#include "AliAnalysisTaskFilteredTree.h"
#include <iostream>

#include "../UnitTests/ConfigConstants.h"


TFile fileMC(MCInputFile);

TString outputFileName = "RootInteractive_" + KFTreeOutputFile;
TFile fileReco(outputFileName);

TTree *treeMC = 0;
TTree * treeKFBR=0;
TTree * treeKFAR=0;

int readTree(){
   
    if (fileMC.IsZombie()) {
        std::cerr << "fileMC: Error opening file" << endl;
        //exit(-1);
        return -1;
    }
    if (fileReco.IsZombie()) {
        std::cerr << "fileReco: Error opening file" << endl;
        //exit(-1);
        return -1;
    }
        
    treeMC = (TTree*)fileMC.Get("Events");
    treeKFAR = (TTree*)fileReco.Get("KFAfterReco");
        
    std::cout << "treeMC:   # events = " << treeMC->GetEntries() << std::endl;
    std::cout << "treeKFAR: # events = " << treeKFAR->GetEntries() << std::endl;

    return 0;   
}