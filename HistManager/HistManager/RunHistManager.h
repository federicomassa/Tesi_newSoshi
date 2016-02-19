#ifndef HistManager_RunHistManager_H
#define HistManager_RunHistManager_H

#include "HistManager/RunHists.h"
#include "HistManager/TrackHistManager.h"
#include <iostream>
#include <vector>

class RunHistManager {

 public:

  RunHistManager(TrackHistManager* trackHistManager);
  ~RunHistManager();

  void Init(EL::Worker* wk);
#ifndef __MAKECINT__
  void FillHists(float) const;
#endif // not __MAKECINT__


  private:

  TString m_name;
  const std::vector<TrackHists*>& m_trackHistList;
  std::vector<RunHists*> m_runHistList;

};

#endif //HistManager_TrackHistManager_H
