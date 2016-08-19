#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>
#include <iostream>

#include "plotStyle.h"

void draw(std::vector<TH1F*>& v, const std::string& xaxis,  const std::string& imageName) {

  TCanvas* c = new TCanvas();

  if (v.size() != 3) {
    std::cout << "DEVE ESSERE 3" << std::endl;
    exit(1);
  }

  // offShell->SetFillColorAlpha(kRed,0.5);
  // offShell->SetFillStyle(3002);
  v[0]->SetTitle("");
  v[0]->Scale(1/v[0]->GetEntries());
  v[0]->GetYaxis()->SetTitle("Entries (normalized)");
  v[0]->GetXaxis()->SetTitle(xaxis.c_str());
  //  v[0]->GetYaxis()->SetRangeUser(0.0,yhigh);
  v[0]->Draw();

  v[1]->Scale(1/v[1]->GetEntries());
  v[1]->Draw("same");

  v[2]->Scale(1/v[2]->GetEntries());
  v[2]->Draw("same");


  TLegend* leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(v[0], "ExtBrl4", "l");
  leg->AddEntry(v[1], "IExtBrl4", "l");
  leg->AddEntry(v[2], "InclBrl4", "l");

  leg->Draw();

  c->Print((imageName+".png").c_str(),"png");

  delete c;
  delete leg;



}

void massResolution(const char* fname1, const char* fname2, const char* fname3) {

  std::vector<TFile*> f;
  f.push_back(new TFile(fname1));
  f.push_back(new TFile(fname2));
  f.push_back(new TFile(fname3));

  std::vector<TH1F*> onShell, offShell, tot,
    onShell27, offShell27, tot27,
    onShell32, offShell32, tot32,
    onShell4, offShell4, tot4;

  for (int i = 0; i < 3; i++) {

    onShell.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass"));
    offShell.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass"));
    tot.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass"));


    onShell27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass27"));
    offShell27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass27"));
    tot27.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass27"));

    onShell32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass32"));
    offShell32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass32"));
    tot32.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass32"));

    onShell4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOnShellMass4"));
    offShell4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoOffShellMass4"));
    tot4.push_back((TH1F*) f[i]->Get("EventHist_all__sigRecoMass4"));

    
    int color;
    if (i == 0)
      color = kRed;
    else if (i == 1)
      color = kGreen;
    else if (i == 2)
      color = kBlue;
  

    plotStyle(onShell.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -20, 20); 
    plotStyle(offShell.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -10, 10); 
    plotStyle(tot.back(), "", 10, "m_{4#mu} [GeV]", "", color, -20, 20); 


  
    plotStyle(onShell27.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -20, 20); 
    plotStyle(onShell32.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -20, 20); 
    plotStyle(onShell4.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -10, 10); 
    
    plotStyle(offShell27.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -10, 10); 
    plotStyle(offShell32.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -10, 10); 
    plotStyle(offShell4.back(), "", 10, "m_{#mu#mu} [GeV]", "", color, -10, 10); 
    
    plotStyle(tot27.back(), "", 10, "m_{4#mu} [GeV]", "", color, -20, 20); 
    plotStyle(tot32.back(), "", 10, "m_{4} [GeV]", "", color, -20, 20); 
    plotStyle(tot4.back(), "", 10, "m_{4} [GeV]", "", color, -30, 30); 
    

  }


  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);


  draw(onShell, "#mu#mu mass - truth on-shell Z boson mass [GeV]", "sigRecoOnShell");
  draw(offShell, "#mu#mu mass - truth off-shell Z boson mass [GeV]", "sigRecoOffShell");
  draw(tot, "4#mu mass - truth Higgs boson mass [GeV]", "sigReco");


  draw(onShell27, "#mu#mu mass - truth on-shell Z boson mass [GeV]", "sigRecoOnShell27");
  draw(onShell32, "#mu#mu mass - truth on-shell Z boson mass [GeV]", "sigRecoOnShell32");
  draw(onShell4, "#mu#mu mass - truth on-shell Z boson mass [GeV]", "sigRecoOnShell4");


  draw(offShell27, "#mu#mu mass - truth off-shell Z boson mass [GeV]", "sigRecoOffShell27");
  draw(offShell32, "#mu#mu mass - truth off-shell Z boson mass [GeV]", "sigRecoOffShell32");
  draw(offShell4, "#mu#mu mass - truth off-shell Z boson mass [GeV]", "sigRecoOffShell4");


  draw(tot27, "4#mu mass - truth Higgs boson mass [GeV]", "sigReco27");
  draw(tot32, "4#mu mass - truth Higgs boson mass [GeV]", "sigReco32");
  draw(tot4, "4#mu mass - truth Higgs boson mass [GeV]", "sigReco4");




}
