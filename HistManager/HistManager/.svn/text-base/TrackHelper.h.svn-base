#ifndef XAODTRACKING_TRACKHELPER_H
#define XAODTRACKING_TRACKHELPER_H

///  
/// Simple set of helpers to make access
/// to complement the existing functionality
/// with the goal of 
///
/// @author Gabriel Facini <Gabriel.Facini@cern.ch>
/// @nosubgrouping
//
/// Wed Sep 17 10:51:35 CEST 2014
///
///
 
// c++ include(s):
#include <iostream>
#include <map>
// EDM include(s):
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackMeasurementValidation.h"
#include "xAODTruth/TruthVertex.h"

#include "AthLinks/ElementLink.h"

namespace xAOD {

  namespace TrackHelper {
    ///@brief number of Si hits
    int numberOfSiHits(const xAOD::TrackParticle *trkPart);
    ///@brief number of dead Si sensors
    int numberOfSiDeadSensors(const xAOD::TrackParticle *trkPart);
    ///@brief chi2 per degree of freedom for the overall track fit
    float chiSqPerDoF(const xAOD::TrackParticle *trkPart);
    ///@brief truth match "probability" ... really a score
    float truthMatchProb(const xAOD::TrackParticle *trkPart);
    ///@brief the truth particle if a link exists
    const xAOD::TruthParticle* truthParticle(const xAOD::TrackParticle *trkPart);
    ///@brief the hits on track if a linke exists
    const std::vector<ElementLink< xAOD::TrackStateValidationContainer > > measureOnTrack(const xAOD::TrackParticle *trkPart);
    ///@brief count incorrectly assigned and missing hits on track
    void countMissingAndWrongHitsOnTrack(const xAOD::TrackParticle *trk, const xAOD::TrackParticle *refTrk );
    ///@brief check if two MSOS correspond to sibling PRDs ( from the same split cluster )
    bool areSiblings( const xAOD::TrackStateValidation* iMSOS, const xAOD::TrackStateValidation* jMSOS );
    ///@brief check if two PRDs are sibling ( from the same split cluster )
    bool areSiblings( const xAOD::TrackMeasurementValidation* iPRD, const xAOD::TrackMeasurementValidation* jPRD);
    ///@brief build a multimap of the cluster siblings
    std::multimap< xAOD::TrackMeasurementValidation*, xAOD::TrackMeasurementValidation* > getSiblingMap( const xAOD::TrackMeasurementValidationContainer* clusters );

    ///@brief Check if the truth particle is a primary charged particle within detector acceptance
    bool isPrimary( const xAOD::TruthParticle*   );
  
    ///@brief  Check if the truth particle is a secondary charged particle within detector acceptance
    bool isSecondary( const xAOD::TruthParticle* );
	
    ///@brief  Check if the truth particle is a combinatorial fake or has no truth associated to it
    bool isFake( const xAOD::TrackParticle*      );
  
    ///@brief Check if the track particle originated from is a primary charged particle within detector acceptance
    bool isPrimary( const xAOD::TrackParticle*   );

    ///@brief  Check if the track particle originated from is a secondary charged particle within detector acceptance
    bool isSecondary( const xAOD::TrackParticle* );

  } // close namespace TrackHelper

} // close namespace xAOD


#endif // XAODTRACKING_TRACKHELPER_H
