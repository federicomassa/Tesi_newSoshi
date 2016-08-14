#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>
#include <iostream>

#include "plotStyle.h"


void muonTracking(const char* fname1, const char* fname2, const char* fname3) {
  
  std::vector<TFile*> f;
  f.push_back(new TFile(fname1));
  f.push_back(new TFile(fname2));
  f.push_back(new TFile(fname3));

  std::vector<TGraphErrors*> sigPt,
    sigQPt,
    sigEta,
    sigPhi,
    sigD0,
    sigZ0;

  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kGreen);
  colors.push_back(kBlue);

  for (int i = 0; i < 3; i++) {
    sigPt.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigPt_abseta"));
    sigQPt.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigQPt_abseta"));
    sigEta.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigEta_abseta"));
    sigPhi.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigPhi_abseta"));
    sigD0.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigD0_abseta"));
    sigZ0.push_back((TGraphErrors*) f[i]->Get("RunHist_reco_all_hard_all__sigZ0_abseta"));

    plotStyle(sigPt.back(), "", "|#eta|", "#sigma(p_{T}) [MeV]", colors[i], 0, 4); 
    plotStyle(sigQPt.back(), "", "|#eta|", "#sigma(q/p_{T}) x p^{T}_{truth}", colors[i], 0, 4); 
    plotStyle(sigEta.back(), "", "|#eta|", "#sigma(#eta)/#eta_{truth}", colors[i], 0, 4); 
    plotStyle(sigPhi.back(), "", "|#eta|", "#sigma(#phi) [rad]", colors[i], 0, 4); 
    plotStyle(sigD0.back(), "", "|#eta|", "#sigma(d_{0}) [mm]", colors[i], 0, 4); 
    plotStyle(sigZ0.back(), "", "|#eta|", "#sigma(z_{0}) [mm]", colors[i], 0, 4); 

  }

  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);

  std::vector<pair<std::vector<TGraphErrors*>, const char*> > hists;
  hists.push_back(std::make_pair(sigPt, "sigPt.png"));
  hists.push_back(std::make_pair(sigQPt, "sigQPt.png"));
  hists.push_back(std::make_pair(sigEta, "sigEta.png"));
  hists.push_back(std::make_pair(sigPhi, "sigPhi.png"));
  hists.push_back(std::make_pair(sigD0, "sigD0.png"));
  hists.push_back(std::make_pair(sigZ0, "sigZ0.png"));


  for (auto itr = hists.begin();
       itr != hists.end();
       itr++) {

    TCanvas* c = new TCanvas();
      
    TMultiGraph* mg = new TMultiGraph;
    TLegend* leg = new TLegend(0.25,0.85, 0.45, 0.65);
    
    leg->AddEntry((*itr).first[0], "ExtBrl4", "lp");
    leg->AddEntry((*itr).first[1], "IExtBrl4", "lp");
    leg->AddEntry((*itr).first[2], "InclBrl4", "lp");
    
    mg->Add(((*itr).first)[0]);
    mg->Add(((*itr).first)[1]);
    mg->Add(((*itr).first)[2]);

    mg->SetTitle("");
    mg->Draw("APE");
    mg->GetXaxis()->SetTitle((*itr).first[0]->GetXaxis()->GetTitle());
    mg->GetXaxis()->SetRangeUser(0,4);
    mg->GetYaxis()->SetTitle((*itr).first[0]->GetYaxis()->GetTitle());
    mg->GetYaxis()->SetTitleOffset(1.35);

    mg->Draw("APE");
    
    leg->Draw();

    c->Print((*itr).second, "png");
    
    delete c;

  }



}
