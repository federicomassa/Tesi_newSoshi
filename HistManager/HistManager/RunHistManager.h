#ifndef HistManager_RunHistManager_H
#define HistManager_RunHistManager_H

#include "HistManager/RunHists.h"
#include "HistManager/TrackHistManager.h"
#include "HistManager/TruthHists.h"
#include "EventLoop/Worker.h"
#include <TGraphErrors.h>
#include <iostream>
#include <vector>

class RunHistManager {

 public:

  RunHistManager(const TrackHistManager* const trackHistManager, const TruthHists* const truthHist);
  ~RunHistManager();

  void Init(EL::Worker* wk);
#ifndef __MAKECINT__
  void FillHists(float) const;
#endif // not __MAKECINT__
  void BookHists();
  
  TGraphErrors* declareGraphErrors(TString cat, TString name, TString xaxis, TString yaxis);
  
  private:
  EL::Worker* m_wk;
  TString m_name;
  const std::vector<TrackHists*>& m_trackHistList;
  const TruthHists* const m_truthHist;
  std::vector<RunHists*> m_runHistList;

  TGraphErrors* m_fakeProb_abseta; //!

};

#endif //HistManager_TrackHistManager_H
