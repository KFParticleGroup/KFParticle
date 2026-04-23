
void CreateInteractiveTest() {
    gROOT->ProcessLine(".x MakeKFParticleTreesForRootInteractive.cxx");
    gROOT->ProcessLine("TPython::Exec( \"import interactiveTest\" )");
   
}