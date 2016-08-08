#include <Utility/PhysicsEvent.h>
#include <Utility/combinations.h>
#include <Utility/DeltaFunctions.h>
#include <Utility/TrackHelper.h>

#include <xAODTruth/TruthVertex.h>

#include <TLorentzVector.h>
#include <TMath.h>
#include <utility>
#include <iostream>
#include <set>

const double PhysicsEvent::m_largeNumber = 10000000;
const double PhysicsEvent::m_ZMass       = 91187.6;

std::vector<xAOD::TrackParticleContainer::const_iterator> PhysicsEvent::GetOnShellReco() const {
  if (!m_isRecoSet) {
    std::cout << "ERROR\t in GetOnShellReco(): reco not set" << std::endl;
    exit(1);
  }
  if (m_recoItr.size() != 4) {
    std::cout << "ERROR\t in GetOnShellReco(): wrong target number" << std::endl;
    exit(1);
  }
    
  /*  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > comb = 
    combinations<xAOD::TrackParticleContainer::const_iterator>(m_recoItr, 2);

  double minResidual = m_largeNumber;
  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator bestComb;

  for (auto combItr = comb.begin();
       combItr != comb.end();
       combItr++) {
    TLorentzVector combVec(0.0,0.0,0.0,0.0);
    double charge = 0.0;

    for (auto itr = (*combItr).begin();
	 itr != (*combItr).end();
	 itr++) {
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM((**itr)->pt(), (**itr)->eta(), (**itr)->phi(), m_decayMass);
      charge += (**itr)->charge();
      combVec += tmp;
    }
    double combMass = combVec.M();
    if (TMath::Abs(combMass - m_ZMass) < minResidual &&
	TMath::Abs(charge) < 1E-6) {
      bestComb = combItr;
      minResidual = TMath::Abs(combMass - m_ZMass);
    }
  }

  return *bestComb;

  */

  std::vector<xAOD::TrackParticleContainer::const_iterator> onShellTracks;
  onShellTracks.push_back(m_recoItr[0]);
  onShellTracks.push_back(m_recoItr[1]);

  return onShellTracks;

}

std::vector<xAOD::TrackParticleContainer::const_iterator> PhysicsEvent::GetOnShellReco(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v,
										       const double intermediateMass,
										       const double elementMass) {

  if (reco_v.size() != 4) {
    std::cout << "ERROR\t in GetOnShellReco(): wrong target number" << std::endl;
    exit(1);
  }

  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > comb =
    combinations<xAOD::TrackParticleContainer::const_iterator>(reco_v, 2);

  double minResidual = m_largeNumber;
  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator bestComb;

  for (auto combItr = comb.begin();
       combItr != comb.end();
       combItr++) {
    TLorentzVector combVec(0.0,0.0,0.0,0.0);
    double charge = 0.0;

    for (auto itr = (*combItr).begin();
         itr != (*combItr).end();
         itr++) {
      TLorentzVector tmp;
      tmp.SetPtEtaPhiM((**itr)->pt(), (**itr)->eta(), (**itr)->phi(), elementMass);
      charge += (**itr)->charge();
      combVec += tmp;
    }
    double combMass = combVec.M();
    if (TMath::Abs(combMass - intermediateMass) < minResidual &&
        TMath::Abs(charge) < 1E-6) {
      bestComb = combItr;
      minResidual = TMath::Abs(combMass - intermediateMass);
    }
  }

  return *bestComb;

}


std::vector<xAOD::TrackParticleContainer::const_iterator> PhysicsEvent::GetOffShellReco() const {

  /*  std::vector<xAOD::TrackParticleContainer::const_iterator> onShellVec = GetOnShellReco();
  std::vector<xAOD::TrackParticleContainer::const_iterator> offShellVec;

  for (auto itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {

    if (std::find(onShellVec.begin(), onShellVec.end(), *itr) == onShellVec.end())
      offShellVec.push_back(*itr);

  }
  
  //debug
  if (offShellVec.size() != 2) {
    std::cout << "ERROR\t in PhysicsEvent::GetOffShellReco(): offShell size is wrong" << std::endl;
    exit(1);
  }

  return offShellVec;
  */

  std::vector<xAOD::TrackParticleContainer::const_iterator> offShellTracks;
  offShellTracks.push_back(m_recoItr[2]);
  offShellTracks.push_back(m_recoItr[3]);

  return offShellTracks;

}

std::vector<xAOD::TrackParticleContainer::const_iterator> PhysicsEvent::GetOffShellReco(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v,
											const double intermediateMass,
											const double elementMass) {

  std::vector<xAOD::TrackParticleContainer::const_iterator> onShellVec = GetOnShellReco(reco_v, intermediateMass, elementMass);
  std::vector<xAOD::TrackParticleContainer::const_iterator> offShellVec;

  for (auto itr = reco_v.begin();
       itr != reco_v.end();
       itr++) {

    if (std::find(onShellVec.begin(), onShellVec.end(), *itr) == onShellVec.end())
      offShellVec.push_back(*itr);

  }

  //debug
  if (offShellVec.size() != 2) {
    std::cout << "ERROR\t in PhysicsEvent::GetOffShellReco(): offShell size is wrong" << std::endl;
    exit(1);
  }

  return offShellVec;

}

int PhysicsEvent::GetNCenteredOffShellMuons(const std::vector<xAOD::TrackParticleContainer::const_iterator>& reco_v,
					    const double intermediateMass,
					    const double elementMass) {

  auto vec = GetOffShellReco(reco_v, intermediateMass, elementMass);
  int nCenteredMuons = 0;
  for (auto itr = vec.begin();
       itr != vec.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < 2.7)
      nCenteredMuons++;
  }

  return nCenteredMuons;

}


int PhysicsEvent::GetNCenteredOffShellMuons() {

  auto vec = GetOffShellReco();
  int nCenteredMuons = 0;
  for (auto itr = vec.begin();
       itr != vec.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < 2.7)
      nCenteredMuons++;
  }

  return nCenteredMuons;

}

/*
void PhysicsEvent::MatchReco() {
  Assert("In MatchReco(): truth not set", m_isTruthSet);
  Assert("In MatchReco(): reco not set", m_isRecoSet);
  Assert("In MatchReco(): reco matched pair already set", !m_isRecoPairSet);

  std::set<xAOD::TrackParticleContainer::const_iterator> recoMatchedSet;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator truthItr = m_truthItr.begin();
       truthItr != m_truthItr.end();
       truthItr++) {
    double minDR = m_largeNumber;
    std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator bestTrk_itr;
    for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator recoItr = m_recoItr.begin();
	 recoItr != m_recoItr.end();
	 recoItr++) {

      if (recoMatchedSet.find((*recoItr)) != recoMatchedSet.end())
	continue;

      double dR = deltaR((**truthItr)->phi(), (**recoItr)->phi(), (**truthItr)->eta(), (**recoItr)->eta());
      if (dR < minDR) {
	minDR = dR;
	bestTrk_itr = recoItr;
      }
    }
    m_recoMatchedPair.push_back(std::make_pair((*truthItr), (*bestTrk_itr)));
    recoMatchedSet.insert((*bestTrk_itr));
  }

  m_isRecoPairSet = true;

}

*/

void MatchRecoIterate(std::vector<xAOD::TruthParticleContainer::const_iterator>& truth,
		      std::vector<xAOD::TrackParticleContainer::const_iterator>& tracks,
		      std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >& pair) {
  
  Assert("In MatchRecoIterate(): sizes not corresponding", truth.size() == tracks.size());
  if (truth.size() == 0)
    return;

  //take every combination, match the one with the minimum distance. Erase from vectors and iterate until vector size is 0
  double minDR = PhysicsEvent::m_largeNumber;
  std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator bestTruth;
  std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator bestTrack;

  for (auto truth_itr = truth.begin();
       truth_itr != truth.end();
       truth_itr++) {
    for (auto track_itr = tracks.begin();
	 track_itr != tracks.end();
	 track_itr++) {
      double dR = deltaR((**truth_itr)->phi(), (**track_itr)->phi(), (**truth_itr)->eta(), (**track_itr)->eta());
      if (dR < minDR) {
	minDR = dR;
	bestTruth = truth_itr;
	bestTrack = track_itr;
      }
    }
  }

  pair.push_back(std::make_pair(*bestTruth, *bestTrack));
  truth.erase(bestTruth);
  tracks.erase(bestTrack);

  MatchRecoIterate(truth, tracks, pair);

}

void PhysicsEvent::MatchReco() {
  Assert("In MatchReco(): truth not set", m_isTruthSet);
  Assert("In MatchReco(): reco not set", m_isRecoSet);
  Assert("In MatchReco(): reco matched pair already set", !m_isRecoPairSet);

  auto tracks = m_recoItr;
  auto truth = m_truthItr;

  MatchRecoIterate(truth, tracks, m_recoMatchedPair);

  Assert("In MatchReco(): reco size wrong", m_recoMatchedPair.size() == 4);

  m_isRecoPairSet = true;

}


double PhysicsEvent::GetTruthMass() {
  Assert("In GetTruthMass(): Truth not set", m_isTruthSet);

  return GetTruthP4().M();
}

/*
double PhysicsEvent::GetTruthMaxDeltaZ() {
  Assert("In GetTruthMaxDeltaZ(): Truth not set", m_isTruthSet);

  double maxDZ = -1;
  std::vector<const xAOD::TruthVertex*> vtx;

  for (auto itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {

    const xAOD::TruthVertex* truthV = (**itr)->prodVtx();
    if (truthV == 0)
      return -1;
    else 
      vtx.push_back(truthV);
  }

  ///////////////

  for (auto itr1 = vtx.begin(); itr1 != vtx.end(); itr1++) {
    for (auto itr2 = vtx.begin(); (itr2 - vtx.begin()) > (itr1 - vtx.begin()) && itr2 != vtx.end(); itr2++) {
      if (TMath::Abs((*itr1)->z() - (*itr2)->z()) > maxDZ)
	maxDZ = TMath::Abs((*itr1)->z() - (*itr2)->z());
    }
  } 

  return maxDZ;

}
*/


double PhysicsEvent::GetMatchedMass() {
  Assert("In GetMatchedMass(): Matched not set", m_isMatchedSet);

  return GetMatchedP4().M();
}

double PhysicsEvent::GetRecoMass()  {
  Assert("In GetRecoMass(): Reco not set", m_isRecoSet);

  return GetRecoP4().M();
}

TLorentzVector PhysicsEvent::GetTruthP4()  {
  Assert("In GetTruthP4(): Truth not set", m_isTruthSet);
  TLorentzVector p(0.0,0.0,0.0,0.0);

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    double px = (**itr)->pt()*TMath::Cos((**itr)->phi());
    double py = (**itr)->pt()*TMath::Sin((**itr)->phi());
    double pz = (**itr)->pt()*TMath::SinH((**itr)->eta());
    TLorentzVector tmp(px,py,pz,TMath::Sqrt(TMath::Power(px,2)+
					    TMath::Power(py,2)+
					    TMath::Power(pz,2)+
					    TMath::Power(m_decayMass,2)));
    p += tmp;
  }

  return p;
}

TLorentzVector PhysicsEvent::GetTruthP4WithPhotons()  {
  Assert("In GetTruthP4WithGamma(): Truth not set", m_isTruthSet);
  Assert("In GetTruthP4WithGamma(): Photons not set", m_arePhotonsSet);
  TLorentzVector p(0.0,0.0,0.0,0.0);

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    double px = (**itr)->pt()*TMath::Cos((**itr)->phi());
    double py = (**itr)->pt()*TMath::Sin((**itr)->phi());
    double pz = (**itr)->pt()*TMath::SinH((**itr)->eta());
    TLorentzVector tmp(px,py,pz,TMath::Sqrt(TMath::Power(px,2)+
					    TMath::Power(py,2)+
					    TMath::Power(pz,2)+
					    TMath::Power(m_decayMass,2)));
    p += tmp;
  }

  for (std::vector<const xAOD::TruthParticle*>::iterator itr = m_primaryPhotons.begin();
       itr != m_primaryPhotons.end();
       itr++) {
    double px = (*itr)->pt()*TMath::Cos((*itr)->phi());
    double py = (*itr)->pt()*TMath::Sin((*itr)->phi());
    double pz = (*itr)->pt()*TMath::SinH((*itr)->eta());
    TLorentzVector tmp(px,py,pz,TMath::Sqrt(TMath::Power(px,2)+
					    TMath::Power(py,2)+
					    TMath::Power(pz,2)));
    p += tmp;
    
  }
  
  return p;
}



TLorentzVector PhysicsEvent::GetMatchedP4()  {
  Assert("In GetMatchedP4(): Matched not set", m_isMatchedSet);
  TLorentzVector p(0.0,0.0,0.0,0.0);

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    double px = (**itr)->pt()*TMath::Cos((**itr)->phi());
    double py = (**itr)->pt()*TMath::Sin((**itr)->phi());
    double pz = (**itr)->pt()*TMath::SinH((**itr)->eta());
    TLorentzVector tmp(px,py,pz,TMath::Sqrt(TMath::Power(px,2)+
					    TMath::Power(py,2)+
					    TMath::Power(pz,2)+
					    TMath::Power(m_decayMass,2)));
    p += tmp;
  }

  return p;
}

TLorentzVector PhysicsEvent::GetRecoP4()  {
  Assert("In GetRecoP4(): Reco not set", m_isRecoSet);
  TLorentzVector p(0.0,0.0,0.0,0.0);

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    double px = (**itr)->pt()*TMath::Cos((**itr)->phi());
    double py = (**itr)->pt()*TMath::Sin((**itr)->phi());
    double pz = (**itr)->pt()*TMath::SinH((**itr)->eta());
    TLorentzVector tmp(px,py,pz,TMath::Sqrt(TMath::Power(px,2)+
					    TMath::Power(py,2)+
					    TMath::Power(pz,2)+
					    TMath::Power(m_decayMass,2)));
    p += tmp;
  }

  return p;
}

TLorentzVector PhysicsEvent::GetRecoTruthP4()  {
  Assert("In GetRecoP4(): Reco not set", m_isRecoSet);
  TLorentzVector p(0.0,0.0,0.0,0.0);

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    const xAOD::TruthParticle* truth = xAOD::TrackHelper::truthParticle(**itr);
    p += truth->p4();
  }

  return p;
}

double PhysicsEvent::GetTruthMaxPt()  {

  double maxPt = -1.0;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if ((**itr)->pt() > maxPt)
      maxPt = (**itr)->pt();
  }

  return maxPt;
}

double PhysicsEvent::GetTruthMinPt()  {

  double minPt = m_largeNumber;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if ((**itr)->pt() < minPt)
      minPt = (**itr)->pt();
  }

  return minPt;
}

const std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator>* PhysicsEvent::GetPairWithMaxEta()  {
  
  if (!m_isPairSet)
    return NULL;
  
  double maxEta = -1.0;
  std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >::iterator maxEtaPairItr;

  for (std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >::iterator itr = m_truthMatchedPair.begin();
       itr != m_truthMatchedPair.end();
       itr++) {
    double eta = TMath::Abs((*((*itr).first))->eta());
    if (eta > maxEta) {
      maxEta = eta;
      maxEtaPairItr = itr;
    }
  }

  return &(*maxEtaPairItr);
  
}

const xAOD::TruthParticle* PhysicsEvent::GetTruthWithMaxEta()  {
  const std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator>* maxEtaPair = GetPairWithMaxEta();
  
  if (maxEtaPair == 0) {
    if (m_isTruthSet) {
      double maxEta = -1.0;
      xAOD::TruthParticleContainer::const_iterator maxEtaTruthItr;
      
      for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
	   itr != m_truthItr.end();
	   itr++) {
	double eta = TMath::Abs((**itr)->eta());
	if (eta > maxEta) {
	  maxEta = eta;
	  maxEtaTruthItr = *itr;
	}
      }
      
      return *maxEtaTruthItr;

    }
  }
  
  else if (maxEtaPair)
    return *(maxEtaPair->first);
  

  return NULL;
}


const xAOD::TrackParticle* PhysicsEvent::GetMatchedWithMaxEta()  {
  const std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator>* maxEtaPair = GetPairWithMaxEta();
  
  if (maxEtaPair)
    return *(maxEtaPair->second);

  return NULL;
}

std::vector<const xAOD::TrackParticle*> PhysicsEvent::GetMatchedFake() {

  std::vector<const xAOD::TrackParticle*> tmp;

  for (auto itr = m_truthMatchedPair.begin();
       itr != m_truthMatchedPair.end();
       itr++) {

    double bias = (*((*itr).first))->charge() - (*((*itr).second))->charge();
    if (TMath::Abs(bias) > 1E-6)
      tmp.push_back(*((*itr).second));

  }


  return tmp;

}


const xAOD::TrackParticle* PhysicsEvent::GetRecoWithMaxEta()  {
  
  double maxEta = -1.0;
  xAOD::TrackParticleContainer::const_iterator maxEtaRecoItr;
  
  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    double eta = TMath::Abs((**itr)->eta());
    if (eta > maxEta) {
      maxEta = eta;
      maxEtaRecoItr = *itr;
    }
  }
  
  if (TMath::Abs(maxEta + 1.0) > 1E-6)
    return *maxEtaRecoItr;
  else
    return NULL;

}

double PhysicsEvent::GetTruthMaxEtaPt()  {
  double maxEta = GetTruthMaxEta();
  double residual = m_largeNumber;
  xAOD::TruthParticleContainer::const_iterator maxEtaItr;
  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if (TMath::Abs(TMath::Abs((**itr)->eta()) - maxEta) < residual) {
      residual = TMath::Abs(TMath::Abs((**itr)->eta()) - maxEta);
      maxEtaItr = *itr;
    }
  }

  return (*maxEtaItr)->pt();
}


double PhysicsEvent::GetMatchedMaxPt()  {

  double maxPt = -1.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if ((**itr)->pt() > maxPt)
      maxPt = (**itr)->pt();
  }

  return maxPt;
}

double PhysicsEvent::GetMatchedMinPt()  {

  double minPt = m_largeNumber;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if ((**itr)->pt() < minPt)
      minPt = (**itr)->pt();
  }

  return minPt;
}


double PhysicsEvent::GetMatchedMaxEtaPt()  {
  double maxEta = GetTruthMaxEta();
  double residual = m_largeNumber;
  xAOD::TrackParticleContainer::const_iterator maxEtaItr;
  for (std::vector<std::pair<xAOD::TruthParticleContainer::const_iterator, xAOD::TrackParticleContainer::const_iterator> >::iterator itr = m_truthMatchedPair.begin();
       itr != m_truthMatchedPair.end();
       itr++) {
    if (TMath::Abs(TMath::Abs((*(*itr).first)->eta()) - maxEta) < residual) {
      residual = TMath::Abs(TMath::Abs((*(*itr).first)->eta()) - maxEta);
      maxEtaItr = (*itr).second;
    }
  }

  return (*maxEtaItr)->pt();
}


double PhysicsEvent::GetRecoMaxPt()  {

  double maxPt = -1.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if ((**itr)->pt() > maxPt)
      maxPt = (**itr)->pt();
  }

  return maxPt;
}

double PhysicsEvent::GetRecoMinPt()  {

  double minPt = m_largeNumber;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if ((**itr)->pt() < minPt)
      minPt = (**itr)->pt();
  }

  return minPt;
}

unsigned int PhysicsEvent::GetTruthOutsideDetector()  {

  unsigned int truthOutside = 0;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) > m_etaMax)
      truthOutside++;
  }

  return truthOutside;
}

double PhysicsEvent::GetTruthMaxEta()  {

  double maxEta = -1.0;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) > maxEta)
      maxEta = TMath::Abs((**itr)->eta());
  }

  return maxEta;
}

double PhysicsEvent::GetTruthMinEta()  {

  double minEta = m_largeNumber;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < minEta)
      minEta = TMath::Abs((**itr)->eta());
  }

  return minEta;
}


double PhysicsEvent::GetMatchedMaxEta()  {

  double maxEta = -1.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) > maxEta)
      maxEta = TMath::Abs((**itr)->eta());
  }

  return maxEta;
}

double PhysicsEvent::GetMatchedMinEta()  {

  double minEta = m_largeNumber;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < minEta)
      minEta = TMath::Abs((**itr)->eta());
  }

  return minEta;
}


double PhysicsEvent::GetRecoMaxEta()  {

  double maxEta = -1.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) > maxEta)
      maxEta = TMath::Abs((**itr)->eta());
  }

  return maxEta;
}

double PhysicsEvent::GetRecoMinEta()  {

  double minEta = m_largeNumber;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if (TMath::Abs((**itr)->eta()) < minEta)
      minEta = TMath::Abs((**itr)->eta());
  }

  return minEta;
}


double PhysicsEvent::GetTruthCharge()  {
  Assert("In GetTruthCharge(): Truth not set", m_isTruthSet);
  double charge = 0.0;

  for (std::vector<xAOD::TruthParticleContainer::const_iterator>::iterator itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) 
    charge += (**itr)->charge();

  return charge;
}

double PhysicsEvent::GetMatchedCharge()  {
  Assert("In GetMatchedCharge(): Matched not set", m_isMatchedSet);
  double charge = 0.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) 
    charge += (**itr)->charge();

  return charge;
}

double PhysicsEvent::GetRecoCharge()  {
  Assert("In GetRecoCharge(): Reco not set", m_isRecoSet);
  double charge = 0.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) 
    charge += (**itr)->charge();

  return charge;
}

double PhysicsEvent::GetRecoTruthCharge()  {
  Assert("In GetRecoTruthCharge(): Reco not set", m_isRecoSet);
  double charge = 0.0;

  for (std::vector<xAOD::TrackParticleContainer::const_iterator>::iterator itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    const xAOD::TruthParticle* truth = xAOD::TrackHelper::truthParticle(**itr);
    charge += truth->charge();
  }
  return charge;
}


double PhysicsEvent::GetTruthSmallestDR()  {
  Assert("In GetTruthSmallestDR(): Truth not set", m_isTruthSet);
  std::vector<std::vector<xAOD::TruthParticleContainer::const_iterator> > doublets = 
    combinations<xAOD::TruthParticleContainer::const_iterator>(m_truthItr, 2);

  double smallestDR = m_largeNumber;
    for (std::vector<std::vector<xAOD::TruthParticleContainer::const_iterator> >::iterator doublets_itr = doublets.begin();
	 doublets_itr != doublets.end(); doublets_itr++) {
      double dR = deltaR((*((*doublets_itr)[0]))->phi(),(*((*doublets_itr)[1]))->phi() ,(*((*doublets_itr)[0]))->eta() ,(*((*doublets_itr)[1]))->eta());
      if (dR < smallestDR)
	smallestDR = dR;
    }
    return smallestDR;
}

const xAOD::TruthParticle* PhysicsEvent::GetTruthMaxZ() {
  Assert("In GetTruthMaxZ(): truth not set", m_isTruthSet);

  xAOD::TruthParticleContainer::const_iterator maxZ_itr;

  double maxZ = -m_largeNumber;

  for (auto itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if ((**itr)->hasProdVtx()) {
      if ((**itr)->prodVtx()->z() > maxZ) {
	maxZ = (**itr)->prodVtx()->z();
	maxZ_itr = *itr;
      }
    }
    
    else {
      if (0.0 > maxZ) {
	maxZ = 0.0;
	maxZ_itr = *itr;
      }  
    }
  }

  return *maxZ_itr;

}

const xAOD::TruthParticle* PhysicsEvent::GetTruthMinZ() {
  Assert("In GetTruthMinZ(): truth not set", m_isTruthSet);

  xAOD::TruthParticleContainer::const_iterator minZ_itr;

  double minZ = m_largeNumber;

  for (auto itr = m_truthItr.begin();
       itr != m_truthItr.end();
       itr++) {
    if ((**itr)->hasProdVtx()) {
      if ((**itr)->prodVtx()->z() < minZ) {
	minZ = (**itr)->prodVtx()->z();
	minZ_itr = *itr;
      }
    }
    
    else {
      if (0.0 < minZ) {
	minZ = 0.0;
	minZ_itr = *itr;
      }  
    }
  }

  return *minZ_itr;

}

const xAOD::TrackParticle* PhysicsEvent::GetMatchedMaxZ() {
  Assert("In GetMatchedMaxZ(): matched not set", m_isMatchedSet);

  xAOD::TrackParticleContainer::const_iterator maxZ_itr;

  double maxZ = -m_largeNumber;

  for (auto itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if ((**itr)->z0() > maxZ) {
      maxZ = (**itr)->z0();
      maxZ_itr = *itr;
    }
    
  }

  return *maxZ_itr;

}


const xAOD::TrackParticle* PhysicsEvent::GetMatchedMinZ() {
  Assert("In GetMatchedMinZ(): matched not set", m_isMatchedSet);

  xAOD::TrackParticleContainer::const_iterator minZ_itr;

  double minZ = m_largeNumber;

  for (auto itr = m_matchedItr.begin();
       itr != m_matchedItr.end();
       itr++) {
    if ((**itr)->z0() < minZ) {
      minZ = (**itr)->z0();
      minZ_itr = *itr;
    }
    
  }

  return *minZ_itr;

}

const xAOD::TrackParticle* PhysicsEvent::GetRecoMaxZ() {
  Assert("In GetRecoMaxZ(): reco not set", m_isRecoSet);

  xAOD::TrackParticleContainer::const_iterator maxZ_itr;

  double maxZ = -m_largeNumber;

  for (auto itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if ((**itr)->z0() > maxZ) {
      maxZ = (**itr)->z0();
      maxZ_itr = *itr;
    }
    
  }

  return *maxZ_itr;

}


const xAOD::TrackParticle* PhysicsEvent::GetRecoMinZ() {
  Assert("In GetRecoMinZ(): reco not set", m_isRecoSet);

  xAOD::TrackParticleContainer::const_iterator minZ_itr;

  double minZ = m_largeNumber;

  for (auto itr = m_recoItr.begin();
       itr != m_recoItr.end();
       itr++) {
    if ((**itr)->z0() < minZ) {
      minZ = (**itr)->z0();
      minZ_itr = *itr;
    }
    
  }

  return *minZ_itr;

}

double PhysicsEvent::GetMatchedSmallestDR()  {
  Assert("In GetMatchedSmallestDR(): Matched not set", m_isMatchedSet);
  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > doublets = 
    combinations<xAOD::TrackParticleContainer::const_iterator>(m_matchedItr, 2);

  double smallestDR = m_largeNumber;
    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator doublets_itr = doublets.begin();
	 doublets_itr != doublets.end(); doublets_itr++) {
      double dR = deltaR((*((*doublets_itr)[0]))->phi(),(*((*doublets_itr)[1]))->phi() ,(*((*doublets_itr)[0]))->eta() ,(*((*doublets_itr)[1]))->eta());
      if (dR < smallestDR)
	smallestDR = dR;
    }
    return smallestDR;
}

double PhysicsEvent::GetRecoSmallestDR()  {
  Assert("In GetRecoSmallestDR(): Reco not set", m_isRecoSet);
  std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> > doublets = 
    combinations<xAOD::TrackParticleContainer::const_iterator>(m_recoItr, 2);

  double smallestDR = m_largeNumber;
    for (std::vector<std::vector<xAOD::TrackParticleContainer::const_iterator> >::iterator doublets_itr = doublets.begin();
	 doublets_itr != doublets.end(); doublets_itr++) {
      double dR = deltaR((*((*doublets_itr)[0]))->phi(),(*((*doublets_itr)[1]))->phi() ,(*((*doublets_itr)[0]))->eta() ,(*((*doublets_itr)[1]))->eta());
      if (dR < smallestDR)
	smallestDR = dR;
    }
    return smallestDR;
}

double PhysicsEvent::GetRecoMatchedMaxDR() {

  double maxDR = -1.0;

  for (auto itr = m_recoMatchedPair.begin();
       itr != m_recoMatchedPair.end();
       itr++) {
    double dR = deltaR((*((*itr).first))->phi(), (*((*itr).second))->phi(), (*((*itr).first))->eta(), (*((*itr).second))->eta());
    if (dR > maxDR) {
      maxDR = dR;
    }
  }

  return maxDR;

}

