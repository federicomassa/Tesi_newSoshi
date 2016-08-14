#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>

#include "plotStyle.h"


void massResolution(const char* higgsFileName) {

  TFile* higgsFile= new TFile(higgsFileName);

  
  TH1F* onShell = (TH1F*) higgsFile->Get("EventHist_all__recoOnShellMass");
  TH1F* offShell = (TH1F*) higgsFile->Get("EventHist_all__recoOffShellMass");
  TH1F* tot = (TH1F*) higgsFile->Get("EventHist_all__recoMass");




  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);


  plotStyle(onShell, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(offShell, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(tot, "", 10, "m_{4#mu} [GeV]", "", kBlue, 0, 180); 

  ///////////////

  TCanvas* c = new TCanvas();

  onShell->SetLineColor(kBlue);
  // offShell->SetFillColorAlpha(kRed,0.5);
  // offShell->SetFillStyle(3002);
  onShell->SetTitle("");
  onShell->GetYaxis()->SetTitle("Entries");
  onShell->GetXaxis()->SetTitle("Reconstructed Z boson mass [GeV]");
  onShell->Draw();

  // onShell->SetFillColorAlpha(kBlue,0.5);
  offShell->SetLineColor(kRed);
  //  onShell->SetFillStyle(3004);
  offShell->Draw("same");

  TLegend* leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(onShell, "on-shell Z", "lp");
  leg->AddEntry(offShell, "off-shell Z", "lp");

  leg->Draw();

  c->Print("recoOnShellVsOffShell.png","png");

  delete c;
  delete leg;
  /////////////////

  c = new TCanvas();
  
  
  tot->SetLineColor(kGreen);
  // tot->SetFillColorAlpha(kGreen-8,0.5);
  // tot->SetFillStyle(3002);
  tot->SetTitle("");
  tot->GetYaxis()->SetTitle("Entries");
  tot->GetXaxis()->SetTitle("Reconstructed four lepton mass [GeV]");
  tot->GetXaxis()->SetRangeUser(40,170);
  tot->Draw();
  
  c->Print("recoMass.png","png");


}
