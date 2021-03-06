#ifndef INVARIANT_MASS_H
#define INVARIANT_MASS_H

#include <vector>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODTruth/TruthParticleContainer.h>

namespace InvariantMass {
  double TrackInvariantMass(const std::vector<xAOD::TrackParticleContainer::const_iterator>& vec, const double elementMass);
  double TruthInvariantMass(const std::vector<xAOD::TruthParticleContainer::const_iterator>& vec);
}

#endif //INVARIANT_MASS
