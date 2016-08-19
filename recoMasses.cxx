#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>
#include <vector>
#include <string>

#include "plotStyle.h"


void recoMasses(const char* higgsFileName) {

  TFile* higgsFile= new TFile(higgsFileName);

  
  TH1F* onShell = (TH1F*) higgsFile->Get("EventHist_all__recoOnShellMass");
  TH1F* offShell = (TH1F*) higgsFile->Get("EventHist_all__recoOffShellMass");
  TH1F* tot = (TH1F*) higgsFile->Get("EventHist_all__recoMass");
  TH1F* pt = (TH1F*) higgsFile->Get("EventHist_all__recoPt");
  TH1F* eta = (TH1F*) higgsFile->Get("EventHist_all__recoEta");
  TH1F* phi = (TH1F*) higgsFile->Get("EventHist_all__recoPhi");
  TH1F* truthPt = (TH1F*) higgsFile->Get("EventHist_all__truthPtWithPhotons");
  TH1F* truthEta = (TH1F*) higgsFile->Get("EventHist_all__truthEtaWithPhotons");
  TH1F* truthPhi = (TH1F*) higgsFile->Get("EventHist_all__truthPhiWithPhotons");
 

  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);


  plotStyle(onShell, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(offShell, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(tot, "", 10, "m_{4#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(pt, "", 10, "m_{4#mu} [GeV]", "", kBlue, 0, 180); 
  plotStyle(eta, "", 10, "m_{4#mu} [GeV]", "", kBlue, 0, 7);
  plotStyle(phi, "", 10, "m_{4#mu} [GeV]", "", kBlue, -3.1416, 3.1416); 

  ///////////////

  TCanvas* c = new TCanvas();

  onShell->SetLineColor(kBlue);
  // offShell->SetFillColorAlpha(kRed,0.5);
  // offShell->SetFillStyle(3002);
  onShell->SetTitle("");
  onShell->GetYaxis()->SetTitle("Entries");
  onShell->GetXaxis()->SetTitle("Reconstructed #mu#mu mass [GeV]");
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

  delete c;
  ///////============

  ///////////////

  c = new TCanvas();

  pt->SetLineColor(kBlue);
  // truthPt->SetFillColorAlpha(kRed,0.5);
  // truthPt->SetFillStyle(3002);
  pt->SetTitle("");
  pt->GetYaxis()->SetTitle("Entries (normalised)");
  pt->GetXaxis()->SetTitle("Reconstructed 4#mu p_{T} [GeV]");
  pt->Scale(1/pt->GetEntries());
  pt->Draw();

  truthPt->SetFillColorAlpha(kRed,0.5);
  truthPt->SetLineColor(kRed);
  truthPt->Scale(1/truthPt->GetEntries());
  truthPt->SetFillStyle(3004);
  truthPt->Draw("hist same");

  leg = new TLegend(0.65,0.85,0.85,0.70);
  leg->AddEntry(pt, "Reco", "lp");
  leg->AddEntry(truthPt, "Truth", "lp");

  leg->Draw();

  c->Print("recoPt.png","png");

  delete c;
  delete leg;

  ///////////////

  c = new TCanvas();

  eta->SetLineColor(kBlue);
  // truthEta->SetFillColorAlpha(kRed,0.5);
  // truthEta->SetFillStyle(3002);
  eta->SetTitle("");
  eta->GetYaxis()->SetTitle("Entries (normalised)");
  eta->GetXaxis()->SetTitle("Reconstructed 4#mu |#eta|");
  eta->Scale(1/eta->GetEntries());
  eta->GetYaxis()->SetRangeUser(0.0,1.1*eta->GetMaximum());
  eta->Draw();

  truthEta->SetFillColorAlpha(kRed,0.5);
  truthEta->SetLineColor(kRed);
  truthEta->Scale(1/truthEta->GetEntries());
  truthEta->SetFillStyle(3004);
  truthEta->Draw("hist same");

  leg = new TLegend(0.65,0.85,0.85,0.70);
  leg->AddEntry(eta, "Reco", "lp");
  leg->AddEntry(truthEta, "Truth", "lp");

  leg->Draw();

  c->Print("recoEta.png","png");

  delete c;
  delete leg;


  ///////////////

  c = new TCanvas();

  phi->SetLineColor(kBlue);
  // truthPhi->SetFillColorAlpha(kRed,0.5);
  // truthPhi->SetFillStyle(3002);
  phi->SetTitle("");
  phi->GetYaxis()->SetTitle("Entries (normalised)");
  phi->GetXaxis()->SetTitle("Reconstructed 4#mu #phi [rad]");
  phi->Scale(1/phi->GetEntries());
  phi->Draw();

  truthPhi->SetFillColorAlpha(kRed,0.5);
  truthPhi->SetLineColor(kRed);
  truthPhi->Scale(1/truthPhi->GetEntries());
  truthPhi->SetFillStyle(3004);
  truthPhi->Draw("hist same");

  leg = new TLegend(0.15,0.35,0.35,0.20);
  leg->AddEntry(phi, "Reco", "lp");
  leg->AddEntry(truthPhi, "Truth", "lp");

  leg->Draw();

  c->Print("recoPhi.png","png");

  delete c;
  delete leg;




}
