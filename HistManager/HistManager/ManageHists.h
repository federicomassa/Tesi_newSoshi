#ifndef HistManager_ManageHists_H
#define HistManager_ManageHists_H

#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <EventLoop/Worker.h>

class ManageHists{

  public:

    ManageHists();
    virtual ~ManageHists();

    virtual void Init(EL::Worker* wk);
    virtual void BookHists() { return; }

    TH1F* declare1D(TString cat, TString name, TString xaxis, Int_t xbins, Float_t xlow, Float_t xhigh);
    TH1F* declare1D(TString cat, TString name, TString xaxis, Int_t xbins, Double_t* xbinVec);

    TH2F* declare2D(TString cat, TString name, TString xaxis, TString yaxis, Int_t xbins, Float_t xlow, Float_t xhigh, Int_t ybins, Float_t ylow, Float_t yhigh);
    TH2F* declare2D(TString cat, TString name, TString xaxis, TString yaxis, Int_t xbins, Double_t* xbinVec, Int_t ybins, Float_t ylow, Float_t yhigh);

    TH3F* declare3D(TString cat, TString name, TString xaxis, TString yaxis, TString zaxis, Int_t xbins, Float_t xlow, Float_t xhigh, Int_t ybins, Float_t ylow, Float_t yhigh, Int_t zbins, Float_t zlow, Float_t zhigh);

  private:

    EL::Worker* m_wk;  //!

};

#endif //HistManager_ManageHists_H
