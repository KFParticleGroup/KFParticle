# KFUnitTest
[KFUnitTest](https://github.com/vsandul/KFUnitTest) is an instrument to check correctness of [KFParticle package](https://github.com/Ivan-Kisel/KFParticle) performance. 

To use this package, [KFParticle package](https://github.com/Ivan-Kisel/KFParticle) and [ROOT](https://root.cern/) are neccessary.

*Despite it is called **KFUnitTest**, actually the package is more of an integration test.*

# Instructions to install
## ROOT installation
The present version of the **KFUnitTest** package has been developed and tested within **ROOT v.6.26.04**. Proper work for earlier versions of **ROOT** is not guaranteed.
To install the **ROOT** follow the instruction in [this link](https://root.cern/install/).

## KFParticle installation
1. Download the neccessary branch of **KFParticle** from one of the repositories to folder `<path-to-your-KF-folder/your-KF-folder>`
1. Go to `<your-KF-folder>` and build the package via 
        
        $ cd your-KF-folder
        $ mkdir build
        $ cd build
        $ cmake ..
        $ make

1. Install KFParticle into <install> directory
        
        $ cd ..
        $ mkdir install
        $ make DESTDIR=<path-to-your-KF-folder>/install install

1. Point the install directory as the directory containing libraries and headers of KFParticle. You have to do this in every used terminal, otherwise just put this line into `~/.bashrc` file.
        
        $ export KFPARTICLE_DIR=<path-to-your-KF-folder>/install
        
## KFUnitTest installation
1. Download **KFnitTest** from [the repository](https://github.com/vsandul/KFUnitTest)
1. * Download the MC dataset you are interested in: for this go to subfolder `<MCData/name-of-dataset>` and execute bashscript `download-data-name-of-dataset.sh`;
   * Or download GEANT4 programs and create your own dataset: go to `<GeantPrograms/>` subfolder and execute bashscript `download-program-name-of-dataset.sh` ([GEANT4](https://geant4.web.cern.ch/) is required); after build and execute the GEANT4 program (see `README.md` inside of `GeantPrograms` folder).

# Configuring KFUnitTest
The main configuration file to configure the workflow of the unit test is the `ConfigConstants.h` in the `<UnitTests>` folder. Here you can define version of KFParticle you use, options of reconstruction (modes and constraints), PDG checks and names of input and output files.

# Running KFUnitTest
Before running KFUnitTest, be sure you have downloaded or generated Monte-Carlo datasets (see **KFUnitTest installation** paragraph) and set right option in the `ConfigConstants.h` file (see **Configuring KFUnitTest** paragraph).

The main unit test macro is in the `<UnitTests>` and calls `RunKFUnitTestMacro`. To perform the unit test, run this macro with **ROOT**:
```bash
$ root RunKFUnitTestMacro.cxx
```
The result of this macro work are 2 `.root` files: `outputKFTrees_<data>.root` and `outputKFHistos_<data>.root`. The first one is the file with trees of `KFParticle` structures - before and after reconstruction (KFBR and KFAR). The second file is the file with QA plots and histograms, i.e. pulls and residuals of quanitites of reconstructed particles. This plots and histograms have to be analyzed by user to get a judgement about the quality of the **KFParticle** package performance.

# Interactive KFUnitTest
The other way to make a quality analisys is to use Interactive KFUnitTest. To use this ability, [RootInteractive](https://pypi.org/project/RootInteractive/) has to be installed.
To perform Interactive KFUnitTest, one may install KFUnitTest, download or generate Monte-Carlo data and configure the test with `ConfigConstants.h` file, then go to folder `<InteractiveTest>` and run `CreateInteractiveTest.cxx` script with **ROOT**:
```bash
$ root CreateInteractiveTest.cxx
```
The result of the script work is 2 `.html` files: `InterTest_pulls.html` and `InterTest_residuals.html` with pulls and residuals plots of particles respectively. This files could be opened with any web-browser.

# Extra KFUnitTest settings
## Additional histograms and plots
To provide user with the ability to add his/her own histograms, `OutputObjInterface` was written.
Shortly, user has to write his/her own class inherited from `OutputObjInterface` class (in the `<UnitTests/OutputObjInterface>` folder); then put the instance of the created class into `OutputMap` structure (`OutputMap.h` file).
More detailed for more detailed example you can find `ExampleOutput.h` file in the `<UnitTests/OutputObjInterface>` folder. Don't foreget to add your class object into `OutputMap.h`! (as it is done for example).

## Covariance matrix generation

# KFUnitTest workflow


