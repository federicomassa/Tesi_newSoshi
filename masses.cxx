#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>

#include "plotStyle.h"


void masses(const char* fname) {

  TFile* f = new TFile(fname);
  
  TH1F* onShell = (TH1F*) f->Get("EventHist_all__truthOnShellMass");
  TH1F* offShell = (TH1F*) f->Get("EventHist_all__truthOffShellMass");
  TH1F* tot = (TH1F*) f->Get("EventHist_all__truthMassWithPhotons");

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

  offShell->SetLineColor(kRed);
  offShell->SetFillColorAlpha(kRed,0.5);
  offShell->SetFillStyle(3002);
  offShell->SetTitle("");
  offShell->GetYaxis()->SetTitle("A.U.");
  offShell->GetXaxis()->SetTitle("Generated Z boson mass [GeV]");
  offShell->Draw("hist");

  onShell->SetFillColorAlpha(kBlue,0.5);
  onShell->SetLineColor(kBlue);
  onShell->SetFillStyle(3004);
  onShell->Draw("hist same");

  TLegend* leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(onShell, "on-shell Z", "lp");
  leg->AddEntry(offShell, "off-shell Z", "lp");

  leg->Draw();

  c->Print("truthOnShellVsOffShell.png","png");

  delete c;
  delete leg;
  /////////////////

  c = new TCanvas();
  
  
  tot->SetLineColor(kGreen);
  tot->SetFillColorAlpha(kGreen-8,0.5);
  tot->SetFillStyle(3002);
  tot->SetTitle("");
  tot->GetYaxis()->SetTitle("A.U.");
  tot->GetXaxis()->SetTitle("Generated four lepton mass [GeV]");
  tot->GetXaxis()->SetRangeUser(40,170);
  tot->Draw("hist");
  
  c->Print("gen4muMass.png","png");


}
