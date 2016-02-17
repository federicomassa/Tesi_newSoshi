#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"

#include "ParticleAnalysis/ZmumuAnalysis.h"

#include <iostream>
#include <string>

// ************* usage :
//  $ twoTrkEffRun name_of_output_directory /path/to/input/dataset/
//
// NB: you must have all your data in /path/to/input/dataset/ with a separate subdirectory for each sample
//  
// TIDEAmbi
//     /mel_data/mmilesi/SimpleSamplesTIDE/SingleRho_RDO/Rel19_RDOs/TIDEAmbi_v*/
// DefaultAmbi
//     /mel_data/mmilesi/SimpleSamplesTIDE/SingleRho_RDO/Rel19_RDOs/DefaultAmbi/"
// DefaultAmbi (no splitting)
//     /mel_data/mmilesi/SimpleSamplesTIDE/SingleRho_RDO/Rel19_RDOs/DefaultAmbi_NoSplitting/

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir  = argv[1];
  std::string inDir      = argv[2];
  char* pileupWeightFile = argv[3];
  int targetLayer        = atoi(argv[4]);

  std::cout << "Input directory    : " << inDir << std::endl;
  std::cout << "Output directory   : " << submitDir << std::endl;
  std::cout << "Pileup weight file : " << pileupWeightFile << std::endl;
  std::cout << "Target Layer       : " << targetLayer << std::endl;

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  SH::SampleHandler sh;
  
  // this takes all the files in directory
  SH::scanDir(sh, inDir); 

  // this takes one single file - FOR TESTING
  //SH::DiskListLocal list(inDir.c_str()); // need this if you want to run on a single file!
  //SH::scanDir (sh, list, "group.det-indet.4246665.ID1._000056.pool.root");

  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString ("nc_tree", "CollectionTree");

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler(sh);

  // Add our analysis to the job:
  ZmumuAnalysis::SetPileupWeightFile(pileupWeightFile);
  ZmumuAnalysis *alg = new ZmumuAnalysis;
  alg -> SetTargetLayer(targetLayer);
  job.algsAdd(alg);

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  // we can use other drivers to run things on the Grid, with PROOF, etc.

  // process the job using the driver
  driver.submit(job, submitDir);

  return 0;
}

