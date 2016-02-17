#ifndef HistManager_ProcessFeaturesInterface_H
#define HistManager_ProcessFeaturesInterface_H

// c++ include(s):
#include <iostream>
#include <vector>

#include "TLorentzVector.h"

class ProcessFeatures {

  public:
    unsigned int	nLeptons;
    std::vector<int> LeptonID;
    std::vector<TLorentzVector> Lepton;

    unsigned int	nJets;
    std::vector<TLorentzVector> Jet;

    ProcessFeatures() {
      nLeptons = -1;
//      LeptonID.clear();
//      Lepton.clear();
      nJets = -1;
//      Jet.clear();
    };

};
#endif //ParticleAnalysis_ProcessFeaturesInterface_H
