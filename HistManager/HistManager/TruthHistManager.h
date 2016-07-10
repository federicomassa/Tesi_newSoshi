#ifndef HistManager_TruthHistManager_H
#define HistManager_TruthHistManager_H

#include "HistManager/TruthHists.h"

#include <iostream>
#include <vector>

class TruthHistManager {

 public:

  TruthHistManager(TString name);
  ~TruthHistManager();

  void Init(EL::Worker* wk);
  void FillHists(const xAOD::TruthParticle* tp, float weight) const;


  private:

    TString m_name;

    TruthHists* m_all;          //!

    std::vector<TruthHists*> m_histList;

 public:
    const std::vector<TruthHists*>& GetList() const {
      return m_histList;
    }
    TString GetName() const {return m_name;}
};

#endif //HistManager_TruthHistManager_H
