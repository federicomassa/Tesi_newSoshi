#ifndef PHYSICS_EVENT_H
#define PHYSICS_EVENT_H

#include <xAODTracking/TrackParticleContainer.h>
#include <xAODTruth/TruthParticleContainer.h>

#include <vector>
#include <utility>
#include <Utility/assert.h>

class TLorentzVector;

class PhysicsEvent {
 private:
  std::vector<xAOD::TruthParticleContainer::const_iterator> m_truthItr;
  std::vector<xAOD::TrackParticleContainer::const_iterator> m_matchedItr;
  std::vector<xAOD::TrackParticleContainer::const_iterator> m_recoItr;
  std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> > m_truthMatchedPair;
  std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> > m_recoMatchedPair;
  std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> > m_recoIsolation;
  unsigned int m_nTarget;  
  double m_etaMax;
  double m_decayMass;
  double m_parentMass;

  
  static const double m_ZMass;

  bool m_isTruthSet, m_isMatchedSet, m_isRecoSet, m_isPairSet, m_isRecoPairSet;
  void MatchReco();

 public:
  static const double m_largeNumber;
  bool isTruthVisible;
  int truthCutCode, recoCutCode, recoCutCode27, recoCutCode32, recoCutCode4;

 PhysicsEvent(const int& nTarget, const double& maxEta, const double& decayMass, const double& parentMass) : m_nTarget(nTarget), m_etaMax(maxEta), m_decayMass(decayMass), m_parentMass(parentMass) {

    Assert("in PhysicsEvent::PhysicsEvent(...): size should be >= 2", 
	   m_nTarget >= 2);

    m_isTruthSet = false; 
    m_isMatchedSet = false;
    m_isRecoSet = false;
    m_isPairSet = false;
    m_isRecoPairSet = false;
    isTruthVisible = false;
    recoCutCode = -1;
    recoCutCode27 = -1;
    recoCutCode32 = -1;
    recoCutCode4 = -1;
  }

 PhysicsEvent(const PhysicsEvent& physicsEvent) {/* : m_nTarget(physicsEvent.GetNTarget()), */
    /* m_etaMax(physicsEvent.GetEtaMax()), */
    /* m_decayMass(physicsEvent.GetDecayMass()), */
    /* m_parentMass(physicsEvent.GetParentMass()), */
    /* m_largeNumber(physicsEvent.GetLargeNumber()) { */

    /* if (physicsEvent.isTruthSet()) */
    /*   m_truthItr = physicsEvent.GetTruth(); */
    /* if (physicsEvent.isMatchedSet()) */
    /*   m_matchedItr = physicsEvent.GetMatched(); */
    /* if (physicsEvent.isRecoSet()) */
    /*   m_recoItr = physicsEvent.GetReco(); */
    /* if (physicsEvent.isPairSet()) */
    /*   m_truthMatchedPair = physicsEvent.GetPair(); */

    /* Assert("in PhysicsEvent::PhysicsEvent(...): size should be >= 2",  */
    /* 	   m_nTarget >= 2); */

    /* m_isTruthSet = physicsEvent.isTruthSet();  */
    /* m_isMatchedSet = physicsEvent.isMatchedSet(); */
    /* m_isRecoSet = physicsEvent.isRecoSet(); */
    /* m_isPairSet = physicsEvent.isPairSet(); */

    *this = physicsEvent;
  }

  PhysicsEvent& operator=(const PhysicsEvent& physicsEvent) {
    m_nTarget = physicsEvent.GetNTarget();
    m_etaMax = physicsEvent.GetEtaMax();
    m_decayMass = physicsEvent.GetDecayMass();
    m_parentMass = physicsEvent.GetParentMass();

    isTruthVisible = physicsEvent.isTruthVisible;
    truthCutCode = physicsEvent.truthCutCode;
    recoCutCode = physicsEvent.recoCutCode;
    recoCutCode27 = physicsEvent.recoCutCode27;
    recoCutCode32 = physicsEvent.recoCutCode32;
    recoCutCode4 = physicsEvent.recoCutCode4;

    m_isTruthSet = physicsEvent.isTruthSet();
    m_isMatchedSet = physicsEvent.isMatchedSet();
    m_isRecoSet = physicsEvent.isRecoSet();
    m_isPairSet = physicsEvent.isPairSet();

    if (physicsEvent.isTruthSet())
      m_truthItr = physicsEvent.GetTruth();
    if (physicsEvent.isMatchedSet())
      m_matchedItr = physicsEvent.GetMatched();
    if (physicsEvent.isRecoSet()) {
      m_recoItr = physicsEvent.GetReco();
      m_recoIsolation = physicsEvent.GetRecoIsolation();
      MatchReco();
    }
    if (physicsEvent.isPairSet())
      m_truthMatchedPair = physicsEvent.GetPair();

    Assert("in PhysicsEvent::PhysicsEvent(...): size should be >= 2", 
	   m_nTarget >= 2);

    return *this;
  }
  
  unsigned int GetNTarget() const {return m_nTarget;}
  double GetEtaMax() const {return m_etaMax;}
  double GetDecayMass() const {return m_decayMass;}
  double GetParentMass() const {return m_parentMass;}
  double GetLargeNumber() const {return m_largeNumber;}

  void SetTruth(const std::vector<xAOD::TruthParticleContainer::const_iterator>& truth) {
    Assert("In PhysicsEvent::SetTruth(): size not corresponding", truth.size() == m_nTarget);
    Assert("In PhysicsEvent::SetTruth(): truth already set", !m_isTruthSet);

    m_truthItr = truth;
    m_isTruthSet = true;

  }
  
  void SetMatched(const std::vector<xAOD::TrackParticleContainer::const_iterator>& matched) {
    Assert("In PhysicsEvent::SetMatched(): size not corresponding", matched.size() == m_nTarget);
    Assert("In PhysicsEvent::SetMatched(): matched already set", !m_isMatchedSet);    
    Assert("In PhysicsEvent::SetMatched(): truth not set yet", m_isTruthSet);

    m_matchedItr = matched;
    m_isMatchedSet = true;
  }
  
  void SetReco(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco, const std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> >& isolation) {
    Assert("In PhysicsEvent::SetReco(): size not corresponding", reco.size() == m_nTarget && isolation.size() == m_nTarget);
    Assert("In PhysicsEvent::SetReco(): reco already set", !m_isRecoSet);    
    Assert("In PhysicsEvent::SetReco(): truth not set yet", m_isTruthSet);
    m_recoItr = reco;
    m_recoIsolation = isolation;
    m_isRecoSet = true;
    MatchReco();
  }

  void SetPair(const std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >& pair) {
    Assert("In PhysicsEvent::SetPair(): size not corresponding", pair.size() == m_nTarget);
    Assert("In PhysicsEvent::SetPair(): pair already set", !m_isPairSet);

    m_truthMatchedPair = pair;
    m_isPairSet = true;

  }

  const std::vector<xAOD::TruthParticleContainer::const_iterator>& GetTruth() const {
    if (!m_isTruthSet) {
      std::cout << "in PhysicsEvent::GetTruth(): truth is not set" << std::endl;
      exit(1);
    }

    return m_truthItr;
  }

  const std::vector<xAOD::TrackParticleContainer::const_iterator>& GetMatched() const {
    if (!m_isMatchedSet) {
      std::cout << "in PhysicsEvent::GetMatched(): matched is not set" << std::endl;
      exit(1);
    }

    return m_matchedItr;
  } 

  const std::vector<xAOD::TrackParticleContainer::const_iterator>& GetReco() const {
    if (!m_isRecoSet) {
      std::cout << "in PhysicsEvent::GetReco(): reco is not set" << std::endl;
      exit(1);
    }

    return m_recoItr;
  }

  const std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> >& GetRecoIsolation() const {
    Assert("In PhysicsEvent::GetRecoIsolation(): reco not set", m_isRecoSet);

    return m_recoIsolation;
  }


  const std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >& GetPair() const {
    if (!m_isPairSet) {
      std::cout << "in PhysicsEvent::GetPair(): pair is not set" << std::endl;
      exit(1);
    }

    return m_truthMatchedPair;
  }

  const std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >& GetRecoPair() const {
    if (!m_isRecoPairSet) {
      std::cout << "in PhysicsEvent::GetRecoPair(): reco pair is not set" << std::endl;
      exit(1);
    }

    return m_recoMatchedPair;
  }

  std::vector<xAOD::TrackParticleContainer::const_iterator> GetOnShellReco() const;
  std::vector<xAOD::TrackParticleContainer::const_iterator> GetOffShellReco() const;

  static std::vector<xAOD::TrackParticleContainer::const_iterator> GetOnShellReco(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v, 
										  const double onShellMass,
										  const double elementMass);  

  static std::vector<xAOD::TrackParticleContainer::const_iterator> GetOffShellReco(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v,
										   const double onShellMass,
										   const double elementMass);

  static int GetNCenteredOffShellMuons(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v,
				       const double onShellMass,
				       const double elementMass);

  int GetNCenteredOffShellMuons();

  const bool& isTruthSet() const {return m_isTruthSet;}
  const bool& isMatchedSet() const {return m_isMatchedSet;}
  const bool& isRecoSet() const {return m_isRecoSet;}
  const bool& isPairSet() const {return m_isPairSet;}
  const bool& isRecoPairSet() const {return m_isRecoPairSet;}

  const std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator>* GetPairWithMaxEta();
  const xAOD::TruthParticle* GetTruthWithMaxEta();
  const xAOD::TrackParticle* GetMatchedWithMaxEta();
  std::vector<const xAOD::TrackParticle*> GetMatchedFake();
  const xAOD::TrackParticle* GetRecoWithMaxEta(); //WARNING: here the max eta refers to the reco track
  
  double GetTruthMass();
  TLorentzVector GetTruthP4();
  double GetTruthCharge();
  double GetTruthSmallestDR();
  double GetTruthMaxPt();
  double GetTruthMinPt();
  double GetTruthMaxEtaPt();
  double GetTruthMaxEta();
  double GetTruthMinEta();
  const xAOD::TruthParticle* GetTruthMaxZ();
  const xAOD::TruthParticle* GetTruthMinZ();

  unsigned int GetTruthOutsideDetector();

  double GetMatchedMass();
  TLorentzVector GetMatchedP4();
  double GetMatchedCharge();
  double GetMatchedSmallestDR();
  double GetMatchedMaxPt();
  double GetMatchedMinPt();
  double GetMatchedMaxEtaPt();
  double GetMatchedMaxEta();
  double GetMatchedMinEta();
  const xAOD::TrackParticle* GetMatchedMaxZ();
  const xAOD::TrackParticle* GetMatchedMinZ();

  double GetRecoMass();
  TLorentzVector GetRecoP4();
  TLorentzVector GetRecoTruthP4();
  double GetRecoCharge();
  double GetRecoTruthCharge();
  double GetRecoSmallestDR();
  double GetRecoMaxPt();
  double GetRecoMinPt();
  double GetRecoMaxEta();
  double GetRecoMinEta();
  double GetRecoMatchedMaxDR();
  const xAOD::TrackParticle* GetRecoMaxZ();
  const xAOD::TrackParticle* GetRecoMinZ();
};


#endif
