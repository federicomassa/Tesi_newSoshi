#ifndef HistManager_RunHists_h
#define HistManager_RunHists_h

#include "HistManager/ManageHists.h"
#include "HistManager/TrackHists.h"
#include "HistManager/TruthHists.h"

class RunHists : public ManageHists {

 public:
  RunHists(const TrackHists* const, const TruthHists* const);
  ~RunHists();

  void BookHists();
  #ifndef __MAKECINT__
  void FillHists(float weight) const;
  #endif

 private:
  TString m_name;
  TString m_label;
  const TrackHists* const m_trackHist;
  const TruthHists* const m_truthHist;

  TGraphErrors* m_sysPt_abseta; //!
  TGraphErrors* m_sysQPt_abseta; //!
  TGraphErrors* m_sysEta_abseta; //!
  TGraphErrors* m_sysPhi_abseta; //!
  TGraphErrors* m_sysD0_abseta; //!
  TGraphErrors* m_sysZ0_abseta; //!

  TGraphErrors* m_sigPt_abseta; //!
  TGraphErrors* m_sigQPt_abseta; //!
  TGraphErrors* m_sigEta_abseta; //!
  TGraphErrors* m_sigPhi_abseta; //!
  TGraphErrors* m_sigD0_abseta; //!
  TGraphErrors* m_sigZ0_abseta; //!

  TGraphErrors* m_sysPt_poseta; //!
  TGraphErrors* m_sysQPt_poseta; //!
  TGraphErrors* m_sysEta_poseta; //!
  TGraphErrors* m_sysPhi_poseta; //!
  TGraphErrors* m_sysD0_poseta; //!
  TGraphErrors* m_sysZ0_poseta; //!

  TGraphErrors* m_sigPt_poseta; //!
  TGraphErrors* m_sigQPt_poseta; //!
  TGraphErrors* m_sigEta_poseta; //!
  TGraphErrors* m_sigPhi_poseta; //!
  TGraphErrors* m_sigD0_poseta; //!
  TGraphErrors* m_sigZ0_poseta; //!

  TGraphErrors* m_sysPt_negeta; //!
  TGraphErrors* m_sysQPt_negeta; //!
  TGraphErrors* m_sysEta_negeta; //!
  TGraphErrors* m_sysPhi_negeta; //!
  TGraphErrors* m_sysD0_negeta; //!
  TGraphErrors* m_sysZ0_negeta; //!

  TGraphErrors* m_sigPt_negeta; //!
  TGraphErrors* m_sigQPt_negeta; //!
  TGraphErrors* m_sigEta_negeta; //!
  TGraphErrors* m_sigPhi_negeta; //!
  TGraphErrors* m_sigD0_negeta; //!
  TGraphErrors* m_sigZ0_negeta; //!

  TGraphErrors* m_eff_abseta; //!
  TGraphErrors* m_eff_phi; //!

  TGraphErrors* m_PixHits; //!
  TGraphErrors* m_SCTHits; //!
  TGraphErrors* m_SiHits; //!

  TGraphErrors* m_chiSqPerDof; //!

};

#endif // end HistManager_RunHists_h
