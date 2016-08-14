#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>

#include "plotStyle.h"


void nTracks(const char* fname) {

  TFile* f = new TFile(fname);
  
  TH1F* nRecotrk1GeV = (TH1F*) f->Get("EventHist_all__nRecotrk1GeV");

  TH1F* nPrimaryChargedTruth1GeV = (TH1F*) f->Get("EventHist_all__nPrimaryChargedTruth1GeV");
  TH1F* nPrimaryChargedTruth = (TH1F*) f->Get("EventHist_all__nPrimaryChargedTruth");
  TH1F* nChargedTruth1GeV = (TH1F*) f->Get("EventHist_all__nChargedTruth1GeV");
  TH1F* nChargedTruth = (TH1F*) f->Get("EventHist_all__nChargedTruth");

  TH1F* nSiHits_matched = (TH1F*) f->Get("TrackHist_matched_all__nSiHits");
  TH1F* nSiHits_pileup = (TH1F*) f->Get("TrackHist_pileup_all__nSiHits");

  gStyle->SetOptStat(10);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);

  std::vector<pair<TH1F*, const char*>> hists;
  hists.push_back(std::make_pair(nRecotrk1GeV, "nRecotrk1GeV.png"));
  hists.push_back(std::make_pair(nPrimaryChargedTruth1GeV, "nPrimaryChargedTruth1GeV.png"));
  hists.push_back(std::make_pair(nPrimaryChargedTruth, "nPrimaryChargedTruth.png"));
  hists.push_back(std::make_pair(nChargedTruth1GeV, "nChargedTruth1GeV.png"));
  hists.push_back(std::make_pair(nChargedTruth, "nChargedTruth.png"));

  plotStyle(nRecotrk1GeV, "", 10, "Number of reconstructed tracks with p_{T} > 1 GeV in event", "", kBlue, 0, 20); 
  nRecotrk1GeV->SetEntries(nRecotrk1GeV->GetEffectiveEntries());

  plotStyle(nPrimaryChargedTruth1GeV, "", 10, "Number of primary charged particles per event with p_{T} > 1 GeV", "", kBlue,0,20); 
  nPrimaryChargedTruth1GeV->SetEntries(nPrimaryChargedTruth1GeV->GetEffectiveEntries());
  plotStyle(nPrimaryChargedTruth, "", 10, "Number of primary charged particles per event", "", kBlue,0,60); 
  nPrimaryChargedTruth->SetEntries(nPrimaryChargedTruth->GetEffectiveEntries());
  plotStyle(nChargedTruth1GeV, "", 10, "Number of charged particles per event with p_{T} > 1 GeV", "", kBlue,0,30); 
  plotStyle(nChargedTruth, "", 10, "Number of charged particles per event", "", kBlue, 0, 160); 



  for (auto itr = hists.begin();
       itr != hists.end();
       itr++) {

    TCanvas* c = new TCanvas();
      
    

    ((*itr).first)->Draw();
    
    c->Print((*itr).second, "png");
    
    delete c;

  }

  ///////////////

  TCanvas* c = new TCanvas();

  nRecotrk1GeV->SetFillColorAlpha(kRed,0.5);
  nRecotrk1GeV->SetLineColor(kRed);
  nRecotrk1GeV->SetFillStyle(3002);
  nRecotrk1GeV->GetXaxis()->SetTitle("Number of particles/tracks per event with p_{T} > 1 GeV");
  nRecotrk1GeV->Draw("hist");

  nPrimaryChargedTruth1GeV->SetLineColor(kBlue);
  nPrimaryChargedTruth1GeV->SetFillColorAlpha(kBlue,0.5);
  nPrimaryChargedTruth1GeV->SetFillStyle(3004);
  nPrimaryChargedTruth1GeV->Draw("hist same");

  TLegend* leg = new TLegend(0.55,0.55,0.85,0.35);
  leg->AddEntry(nRecotrk1GeV, "Reconstructed tracks", "lp");
  leg->AddEntry(nPrimaryChargedTruth1GeV, "Primary charged particles", "lp");

  leg->Draw();

  c->Print("nPrimaryChargedTruth_recotrk1GeV.png","png");

  delete c;
  delete leg;
  /////////////////

  gStyle->SetOptStat(0);

  c = new TCanvas();

  nSiHits_matched->SetFillColorAlpha(kRed,0.5);
  nSiHits_matched->SetLineColor(kRed);
  nSiHits_matched->SetFillStyle(3002);
  nSiHits_matched->SetTitle("");
  nSiHits_matched->GetXaxis()->SetTitle("Number of total silicon hits");
  nSiHits_matched->Draw("hist");

  nSiHits_pileup->SetLineColor(kBlue);
  nSiHits_pileup->SetFillColorAlpha(kBlue,0.5);
  nSiHits_pileup->SetFillStyle(3004);
  nSiHits_pileup->Draw("hist same");

  leg = new TLegend(0.15,0.85,0.35,0.70);
  leg->AddEntry(nSiHits_matched, "Matched tracks", "lp");
  leg->AddEntry(nSiHits_pileup, "Pile-up tracks", "lp");

  leg->Draw();

  c->Print("nSiHits_matchedVsPileup.png","png");




}
