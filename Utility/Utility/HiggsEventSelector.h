#ifndef HIGGSEVENTSELECTOR_H
#define HIGGSEVENTSELECTOR_H

#include <vector>
#include <utility>
#include <xAODTracking/TrackParticleContainer.h>

class HiggsEventSelector {

 private:
  std::vector<xAOD::TrackParticleContainer::const_iterator>& candidate;
  std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> > m_candidateIsolation;
  bool m_isRankedPtCutSet;
  bool m_isOnShellMassCutSet, m_isOffShellMassCutSet;
  bool m_isCenteredMuonCutSet;
  bool m_isCenteredOffShellMuonCutSet;
  bool m_isDRCutSet;
  bool m_isIsolationCutSet;
  std::vector<double> m_rankedPtCut;
  double m_onShellMassCutLow, m_onShellMassCutHigh;
  double m_offShellMassCutLow, m_offShellMassCutHigh;
  double m_centeredMuonCutEtaLow, m_centeredMuonCutEtaHigh;
  double m_centeredOffShellMuonCutEtaLow, m_centeredOffShellMuonCutEtaHigh;
  int m_centeredMuonCutN;
  double m_dRCut;
  double m_isolationCut;

  std::vector<xAOD::TrackParticleContainer::const_iterator> ChooseOffShellCombination(const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>, std::vector<xAOD::TrackParticleContainer::const_iterator> >&,
										      const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>, std::vector<xAOD::TrackParticleContainer::const_iterator> >&) const;

  std::vector<xAOD::TrackParticleContainer::const_iterator> ChooseOnShellCombination(const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>, std::vector<xAOD::TrackParticleContainer::const_iterator> >&,
										      const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>, std::vector<xAOD::TrackParticleContainer::const_iterator> >&) const;


 public:

  static double m_ZMass;
  static double m_muMass;

  HiggsEventSelector(std::vector<xAOD::TrackParticleContainer::const_iterator>&);
  int Check() const;
  int CheckRankedPtCut() const;
  void CheckOnShellMassCut(const std::vector<xAOD::TrackParticleContainer::const_iterator>&, int&) const;
  void CheckOffShellMassCut(const std::vector<xAOD::TrackParticleContainer::const_iterator>&, int&) const;

  int CheckOnShellMassCut(const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>,
			   std::vector<xAOD::TrackParticleContainer::const_iterator> >&,
			   const std::pair<std::vector<xAOD::TrackParticleContainer::const_iterator>,
			   std::vector<xAOD::TrackParticleContainer::const_iterator> >&, int&) const;
  bool CheckCenteredMuonCut() const;
  void CheckCenteredOffShellMuonCut(const std::vector<xAOD::TrackParticleContainer::const_iterator>&, int&) const;
  bool CheckDRCut() const;
  bool CheckIsolationCut() const;
  void SetRankedPtCut(const std::vector<double>&);
  void SetOnShellMassCut(const double&, const double&);
  void SetOffShellMassCut(const double&, const double&);
  void SetCenteredMuonCut(const double&, const double&, const int&);
  void SetCenteredOffShellMuonCut(const double&, const double&);
  void SetDRCut(const double&);
  void SetIsolationCut(const std::vector<std::pair<xAOD::TrackParticleContainer::const_iterator, double> >&, const double&);
};

#endif //HIGGSEVENTSELECTOR_H
