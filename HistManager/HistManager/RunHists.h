#ifndef HistManager_RunHists_h
#define HistManager_RunHists_h

#include "HistManager/ManageHists.h"
#include "HistManager/TrackHists.h"

class RunHists : public ManageHists {

 public:
 RunHists(const TrackHists* const);
  ~RunHists();

  void BookHists();
  #ifndef __MAKECINT__
  void FillHists(float weight) const;
  #endif

 private:
  TString m_name;
  TString m_label;
  const TrackHists* const m_trackHist;
  TH1F* m_testHist; //!
};

#endif
