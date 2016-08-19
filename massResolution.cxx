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

  TH1F* onShell27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOnShellMass27");
  TH1F* offShell27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOffShellMass27");
  TH1F* tot27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoMass27");
  TH1F* pt27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPt27");
  TH1F* eta27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoEta27");
  TH1F* phi27 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPhi27");

  TH1F* onShell32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOnShellMass32");
  TH1F* offShell32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOffShellMass32");
  TH1F* tot32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoMass32");
  TH1F* pt32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPt32");
  TH1F* eta32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoEta32");
  TH1F* phi32 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPhi32");

  TH1F* onShell4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOnShellMass4");
  TH1F* offShell4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoOffShellMass4");
  TH1F* tot4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoMass4");
  TH1F* pt4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPt4");
  TH1F* eta4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoEta4");
  TH1F* phi4 = (TH1F*) higgsFile->Get("EventHist_all__sigRecoPhi4");


  gStyle->SetOptStat(0);
  gStyle->SetStatX(0.8);
  gStyle->SetStatY(0.8);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);


  plotStyle(onShell27, "", 10, "m_{#mu#mu} [GeV]", "", kBlack, -20, 20); 
  plotStyle(onShell32, "", 10, "m_{#mu#mu} [GeV]", "", kRed, -20, 20); 
  plotStyle(onShell4, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, -20, 20); 

  ///////////////

  TCanvas* c = new TCanvas();

  // offShell->SetFillColorAlpha(kRed,0.5);
  // offShell->SetFillStyle(3002);
  onShell27->SetTitle("");
  onShell27->Scale(1/onShell27->GetEntries());
  onShell27->GetYaxis()->SetTitle("Entries (normalized)");
  onShell27->GetXaxis()->SetTitle("On-shell #mu#mu mass - truth on-shell Z mass [GeV]");
  onShell27->GetYaxis()->SetRangeUser(0.0,0.16);
  onShell27->Draw();

  onShell32->Scale(1/onShell32->GetEntries());
  onShell32->Draw("same");

  onShell4->Scale(1/onShell4->GetEntries());
  onShell4->Draw("same");


  TLegend* leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(onShell27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(onShell32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(onShell4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoOnShellMass.png","png");

  delete c;
  delete leg;
  /////////////////

  plotStyle(offShell27, "", 10, "m_{#mu#mu} [GeV]", "", kBlack, -10, 10); 
  plotStyle(offShell32, "", 10, "m_{#mu#mu} [GeV]", "", kRed, -10, 10); 
  plotStyle(offShell4, "", 10, "m_{#mu#mu} [GeV]", "", kBlue, -10, 10); 


  c = new TCanvas();

  // offShell->SetFillColorAlpha(kRed,0.5);
  // offShell->SetFillStyle(3002);
  offShell27->SetTitle("");
  offShell27->Scale(1/offShell27->GetEntries());
  offShell27->GetYaxis()->SetTitle("Entries (normalized)");
  offShell27->GetXaxis()->SetTitle("Off-shell #mu#mu mass - truth off-shell Z mass [GeV]");
  //  offShell27->GetYaxis()->SetRangeUser(0.0,0.18);
  offShell27->Draw();

  offShell32->Scale(1/offShell32->GetEntries());
  offShell32->Draw("same");

  offShell4->Scale(1/offShell4->GetEntries());
  offShell4->Draw("same");


  leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(offShell27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(offShell32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(offShell4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoOffShellMass.png","png");

  delete c;
  delete leg;
  
  ////////////===========

  plotStyle(tot27, "", 10, "m_{4#mu} [GeV]", "", kBlack, -20, 20); 
  plotStyle(tot32, "", 10, "m_{4#mu} [GeV]", "", kRed, -20, 20); 
  plotStyle(tot4, "", 10, "m_{4#mu} [GeV]", "", kBlue, -20, 20); 


  c = new TCanvas();

  // tot->SetFillColorAlpha(kRed,0.5);
  // tot->SetFillStyle(3002);
  tot27->SetTitle("");
  tot27->Scale(1/tot27->GetEntries());
  tot27->GetYaxis()->SetTitle("Entries (normalized)");
  tot27->GetXaxis()->SetTitle("4#mu mass - truth Higgs mass [GeV]");
  tot27->GetYaxis()->SetRangeUser(0.0,0.16);
  tot27->Draw();

  tot32->Scale(1/tot32->GetEntries());
  tot32->Draw("same");

  tot4->Scale(1/tot4->GetEntries());
  tot4->Draw("same");


  leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(tot27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(tot32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(tot4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoMass.png","png");

  delete c;
  delete leg;
  
  ///////==================

  plotStyle(pt27, "", 10, "m_{4#mu} [GeV]", "", kBlack, -20, 20); 
  plotStyle(pt32, "", 10, "m_{4#mu} [GeV]", "", kRed, -20, 20); 
  plotStyle(pt4, "", 10, "m_{4#mu} [GeV]", "", kBlue, -20, 20); 


  c = new TCanvas();

  // pt->SetFillColorAlpha(kRed,0.5);
  // pt->SetFillStyle(3002);
  pt27->SetTitle("");
  pt27->Scale(1/pt27->GetEntries());
  pt27->GetYaxis()->SetTitle("Entries (normalized)");
  pt27->GetXaxis()->SetTitle("4#mu p_{T} - truth Higgs p_{T} [GeV]");
  pt27->GetYaxis()->SetRangeUser(0.0,1.1*pt27->GetMaximum());
  pt27->Draw();

  pt32->Scale(1/pt32->GetEntries());
  pt32->Draw("same");

  pt4->Scale(1/pt4->GetEntries());
  pt4->Draw("same");


  leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(pt27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(pt32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(pt4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoPt.png","png");

  delete c;
  delete leg;
  
  //////=========

  plotStyle(eta27, "", 10, "m_{4#mu} [GeV]", "", kBlack, -0.4,0.4); 
  plotStyle(eta32, "", 10, "m_{4#mu} [GeV]", "", kRed, -0.4,0.4); 
  plotStyle(eta4, "", 10, "m_{4#mu} [GeV]", "", kBlue, -0.4,0.4); 


  c = new TCanvas();

  // eta->SetFillColorAlpha(kRed,0.5);
  // eta->SetFillStyle(3002);
  eta27->SetTitle("");
  eta27->Scale(1/eta27->GetEntries());
  eta27->GetYaxis()->SetTitle("Entries (normalized)");
  eta27->GetXaxis()->SetTitle("4#mu #eta - truth Higgs #eta");
  eta27->GetYaxis()->SetRangeUser(0.0,1.1*eta27->GetMaximum());
  eta27->Draw();

  eta32->Scale(1/eta32->GetEntries());
  eta32->Draw("same");

  eta4->Scale(1/eta4->GetEntries());
  eta4->Draw("same");


  leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(eta27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(eta32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(eta4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoEta.png","png");

  delete c;
  delete leg;

  /////===========

  plotStyle(phi27, "", 10, "m_{4#mu} [GeV]", "", kBlack, -0.2,0.2); 
  plotStyle(phi32, "", 10, "m_{4#mu} [GeV]", "", kRed, -0.2,0.2); 
  plotStyle(phi4, "", 10, "m_{4#mu} [GeV]", "", kBlue, -0.2,0.2); 


  c = new TCanvas();

  // phi->SetFillColorAlpha(kRed,0.5);
  // phi->SetFillStyle(3002);
  phi27->SetTitle("");
  phi27->Scale(1/phi27->GetEntries());
  phi27->GetYaxis()->SetTitle("Entries (normalized)");
  phi27->GetXaxis()->SetTitle("4#mu #phi - truth Higgs #phi [rad]");
  phi27->GetYaxis()->SetRangeUser(0.0,1.1*phi27->GetMaximum());
  phi27->Draw();

  phi32->Scale(1/phi32->GetEntries());
  phi32->Draw("same");

  phi4->Scale(1/phi4->GetEntries());
  phi4->Draw("same");


  leg = new TLegend(0.15,0.85,0.35,0.65);
  leg->AddEntry(phi27, "|#eta^{#mu}_{max}| < 2.7", "l");
  leg->AddEntry(phi32, "2.7 < |#eta^{#mu}_{max}| < 3.2", "l");
  leg->AddEntry(phi4, "3.2 < |#eta^{#mu}_{max}| < 4.0", "l");

  leg->Draw();

  c->Print("sigRecoPhi.png","png");

  delete c;
  delete leg;




}
