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
  TGraphErrors* m_sysPhi_abseta; //!
  TGraphErrors* m_sysD0_abseta; //!
  TGraphErrors* m_sysZ0_abseta; //!

  TGraphErrors* m_sigPt_abseta; //!
  TGraphErrors* m_sigQPt_abseta; //!
  TGraphErrors* m_sigPhi_abseta; //!
  TGraphErrors* m_sigD0_abseta; //!
  TGraphErrors* m_sigZ0_abseta; //!

  TGraphErrors* m_eff_abseta; //!
  TGraphErrors* m_eff_phi; //!

};

#endif // end HistManager_RunHists_h
