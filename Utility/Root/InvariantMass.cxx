#include <Utility/InvariantMass.h>
#include <TLorentzVector.h>


//computes invariant mass of vector of tracks, which every track associated to same-mass particles
double InvariantMass::TrackInvariantMass(const std::vector<xAOD::TrackParticleContainer::const_iterator>& vec, const double elementMass) {

  TLorentzVector total(0.0,0.0,0.0,0.0);

  for (auto itr = vec.begin();
       itr != vec.end();
       itr++) {
    TLorentzVector tmp;
    tmp.SetPtEtaPhiM((**itr)->pt(), (**itr)->eta(), (**itr)->phi(), elementMass);

    total += tmp;
  }

  return total.M();

} 

double InvariantMass::TruthInvariantMass(const std::vector<xAOD::TruthParticleContainer::const_iterator>& vec) {

  TLorentzVector total(0.0,0.0,0.0,0.0);

  for (auto itr = vec.begin();
       itr != vec.end();
       itr++) {
    total += (**itr)->p4();
  }

  return total.M();

}

