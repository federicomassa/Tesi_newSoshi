/******************************************
 *
 * Base class used to book a set of histograms.  
 * Many sets can be created and compared for 
 * efficiecny studies or population studies.
 *
 * G. Facini
 * Sep 9 08:58:12 CEST 2014
 *
 ******************************************/

#include "HistManager/ManageHists.h"
#include "HistManager/TrackHelper.h"

ManageHists::ManageHists() {}

ManageHists::~ManageHists() {}

void ManageHists::Init(EL::Worker* wk) {

  m_wk = wk;
  BookHists();
  //m_wk = 0;

}

TH1F* ManageHists::declare1D(TString cat, TString name, TString xaxis, Int_t xbins, Float_t xlow, Float_t xhigh) {
  TH1F* tmp = new TH1F(TString(cat+"__"+name),name,xbins,xlow,xhigh);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->Sumw2();
  m_wk->addOutput(tmp);
  return tmp;
}; // declare1D

TH1F* ManageHists::declare1D(TString cat, TString name, TString xaxis, Int_t xbins, Double_t* xbinVec) {
  TH1F* tmp = new TH1F(TString(cat+"__"+name),name,xbins,xbinVec);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->Sumw2();
  m_wk->addOutput(tmp);
  return tmp;
}; // declare1D



TH2F* ManageHists::declare2D(TString cat, TString name, TString xaxis, TString yaxis, Int_t xbins, Float_t xlow, Float_t xhigh, Int_t ybins, Float_t ylow, Float_t yhigh) {
  TH2F* tmp = new TH2F(TString(cat+"__"+name),name,xbins,xlow,xhigh,ybins,ylow,yhigh);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->GetYaxis()->SetTitle(yaxis);
  tmp->Sumw2();
  m_wk->addOutput(tmp);
  return tmp;
}; // declare2D

TH2F* ManageHists::declare2D(TString cat, TString name, TString xaxis, TString yaxis, Int_t xbins, Double_t* xbinVec, Int_t ybins, Float_t ylow, Float_t yhigh) {
  TH2F* tmp = new TH2F(TString(cat+"__"+name),name,xbins,xbinVec,ybins,ylow,yhigh);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->GetYaxis()->SetTitle(yaxis);
  tmp->Sumw2();
  m_wk->addOutput(tmp);
  return tmp;
}; // declare2D

TH3F* ManageHists::declare3D(TString cat, TString name, TString xaxis, TString yaxis, TString zaxis, Int_t xbins, Float_t xlow, Float_t xhigh, Int_t ybins, Float_t ylow, Float_t yhigh, Int_t zbins, Float_t zlow, Float_t zhigh) {
  TH3F* tmp = new TH3F(TString(cat+"__"+name),name,xbins,xlow,xhigh,ybins,ylow,yhigh,zbins,zlow,zhigh);
  tmp->GetXaxis()->SetTitle(xaxis);
  tmp->GetYaxis()->SetTitle(yaxis);
  tmp->GetZaxis()->SetTitle(zaxis);
  tmp->Sumw2();
  m_wk->addOutput(tmp);
  return tmp;
}; // declare3D
