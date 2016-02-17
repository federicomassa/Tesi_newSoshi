#ifndef HistManager_ProcessHists_H
#define HistManager_ProcessHists_H

#include <HistManager/ProcessFeaturesInterface.h>
#include <HistManager/ManageHists.h>

class ProcessHists : public ManageHists {

  public:

    ProcessHists(TString name);
    ~ProcessHists();

    void BookHists();
    void FillHists(ProcessFeatures &proc, float weight) const; 

  private:

    TString m_name;
    TString m_label;

  public:

    TH1F* m_numLep; //!
    TH1F* m_ptLep1; //!
    TH1F* m_ptLep2; //!
    TH1F* m_etaLep1; //!
    TH1F* m_etaLep2; //!
    TH1F* m_phiLep1; //!
    TH1F* m_phiLep2; //!

    TH1F* m_massLep12; //!
    TH1F* m_ptLep12; //!
    TH1F* m_dRLep12; //!
    TH1F* m_dEtaLep12; //!

    TH1F* m_numJet; //!
    TH1F* m_ptJet1; //!
    TH1F* m_etaJet1; //!
    TH1F* m_phiJet1; //!

};

#endif //HistManager_ProcessHists_H
