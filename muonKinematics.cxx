#include <TH1F.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>

#include "plotStyle.h"


void muonKinematics(const char* fname) {

  TFile* f = new TFile(fname);
  
  TH1F* muonPt = (TH1F*) f->Get("TruthHist_hard_all__pt");
  TH1F* muonEta = (TH1F*) f->Get("TruthHist_hard_all__eta");
  TH1F* muonPhi = (TH1F*) f->Get("TruthHist_hard_all__phi");
  TH1F* truthOutsideDetector = (TH1F*) f->Get("EventHist_all__truthOutsideDetector");

  gStyle->SetOptStat(10);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);

  std::vector<pair<TH1F*, const char*>> hists;
  hists.push_back(std::make_pair(muonPt, "muonPt.png"));
  hists.push_back(std::make_pair(muonEta, "muonEta.png"));
  hists.push_back(std::make_pair(muonPhi, "muonPhi.png"));
  hists.push_back(std::make_pair(truthOutsideDetector, "truthOutsideDetector.png"));

  plotStyle(muonPt, "", 10, "p_{T,#mu} [GeV]", "A.U.", kBlue, 0, 150); 
  plotStyle(muonEta, "", 10, "#eta_{#mu}", "A.U.", kBlue); 
  plotStyle(muonPhi, "", 10, "#phi_{#mu}", "A.U.", kBlue); 
  plotStyle(truthOutsideDetector, "", 10, "Number of muons per event lying outside of the acceptance", "A.U.", kBlue, 0, 5); 



  for (auto itr = hists.begin();
       itr != hists.end();
       itr++) {

    TCanvas* c = new TCanvas();
      
    if ((*itr).first == muonPhi) {
      (*itr).first->GetYaxis()->SetRangeUser(0.0,6500);
    }

    

    ((*itr).first)->Draw();
    

    if ((*itr).first == truthOutsideDetector)
      (c->cd())->SetLogy(1);

    c->Print((*itr).second, "png");
    
    delete c;

  }


}
