#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TAxis.h>
#include <TH1F.h>


void hitPlots(const char* filename) {

  TFile* f = new TFile(filename);

  // HITS PLOTS

  TGraphErrors* PixHits = dynamic_cast<TGraphErrors*>(f->Get("RunHist_reco_all_hard_all__PixHits"));
  TGraphErrors* SCTHits = dynamic_cast<TGraphErrors*>(f->Get("RunHist_reco_all_hard_all__SCTHits"));
  TGraphErrors* SiHits = dynamic_cast<TGraphErrors*>(f->Get("RunHist_reco_all_hard_all__SiHits"));

  TMultiGraph* mg = new TMultiGraph;
  mg->SetTitle("Number of hits per track");
  
  mg->Add(PixHits);
  mg->Add(SCTHits);
  mg->Add(SiHits);
  
  TCanvas* c = new TCanvas();

  PixHits->SetLineColor(kRed);
  PixHits->SetMarkerColor(kRed);
  
  SCTHits->SetLineColor(kBlue);
  SCTHits->SetMarkerColor(kBlue);
  
  SiHits->SetLineColor(kGreen);
  SiHits->SetMarkerColor(kGreen);
  
  mg->Draw("APE");
  mg->GetXaxis()->SetTitle("|#eta|");
  mg->GetYaxis()->SetTitle("hits");
  mg->GetYaxis()->SetRangeUser(0.0,25.0);

  TLegend* lg = new TLegend(0.15, 0.89,0.35,0.75);
  lg->AddEntry(PixHits, "PixHits", "lp");
  lg->AddEntry(SCTHits, "SCTHits", "lp");
  lg->AddEntry(SiHits, "SiHits", "lp");

  lg->SetFillColor(0);
  lg->Draw();

  //////////////////


  // EFF PLOTS

  TFile* f1 = new TFile("TestInclBrl4/hist-50K_pu50.root");
  TFile* f2 = new TFile("TestInclBrl4/hist-50K_pu200.root");
  //  TFile* f3 = new TFile("TestInclBrl4/hist-50K_pu200.root");


  TH1F* Eff1 = new TH1F("Eff1", "Efficiency vs |#eta|;|#eta|", 20,0,4);
  TH1F* Eff2 = new TH1F("Eff2", "Efficiency", 20,0,4);
  //  TH1F* Eff3 = new TH1F("Eff3", "Efficiency", 20,0,4);


  TH1F* num1 = dynamic_cast<TH1F*>(f1->Get("EventHist_all__reco_abseta"));
  TH1F* den1 = dynamic_cast<TH1F*>(f1->Get("EventHist_all__visibleTruth_abseta"));
  
  TH1F* num2 = dynamic_cast<TH1F*>(f2->Get("EventHist_all__reco_abseta"));
  TH1F* den2 = dynamic_cast<TH1F*>(f2->Get("EventHist_all__visibleTruth_abseta"));

  //  TH1F* num3 = dynamic_cast<TH1F*>(f3->Get("EventHist_all__reco_abseta"));
  // TH1F* den3 = dynamic_cast<TH1F*>(f3->Get("EventHist_all__visibleTruth_abseta"));

  Eff1->Divide(num1,den1,1,1,"B");
  Eff2->Divide(num2,den2,1,1,"B");
  //  Eff3->Divide(num3,den3,1,1,"B");

  TCanvas* c_eff = new TCanvas();

  Eff1->SetLineColor(kRed);
  Eff1->SetMarkerColor(kRed);
  Eff1->GetYaxis()->SetRangeUser(0.0,1.0);
  
  Eff2->SetLineColor(kBlue);
  Eff2->SetMarkerColor(kBlue);
  
  //  Eff3->SetLineColor(kGreen);
  //  Eff3->SetMarkerColor(kGreen);
  
  Eff1->Draw();
  Eff2->Draw("SAME");
  //  Eff3->Draw("SAME");


  TLegend* lg_eff = new TLegend(0.12, 0.32,0.32,0.12);
  lg_eff->AddEntry(Eff1, "InclBrl4: pu50", "lp");
  lg_eff->AddEntry(Eff2, "InclBrl4: pu200", "lp");
  //  lg_eff->AddEntry(Eff3, "InclBrl4: pu200", "lp");

  lg_eff->SetFillColor(0);
  lg_eff->Draw();

  //////////////////
  

}
