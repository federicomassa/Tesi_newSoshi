#include <TCanvas.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>
#include <string>
#include <sstream>
#include <vector>
#include "utils/checkPtr.h"
#include "utils/assert.h"


using namespace std;

void fixedPt(const string& particle, const double& pt, const TString& fixOption, const TString& value) {
  //fixOption = "PILEUP"
  //          = "LAYOUT"

  //value     = "50", "140", ... for fixOption = "PILEUP"
  //          = "ExtBrl4", "IExtBrl4", ... for fixOption = "LAYOUT"

  string particle_lower, particle_upper;
  bool fixPileup, fixLayout;

  if (particle.compare("pi") == 0) {
    particle_lower = particle;
    particle_upper = "Pi";
  }

  else if (particle.compare("mu") == 0) {
    particle_lower = particle;
    particle_upper = "Mu";
  }

  else {
    cout << "ERROR: particle not recognized" << endl;
    return;
  }


  if (fixOption.EqualTo("PILEUP")) {
    fixPileup = true;
    fixLayout = false;
  }
  else if (fixOption.EqualTo("LAYOUT")) {
    fixPileup = false;
    fixLayout = true;
  }
  else {
    cout << "ERROR: fixOption not recognized" << endl;
    return;
  }
    
  string pt_str;

  stringstream ss;
  ss << pt/1000; //to GeV
  ss >> pt_str;
  ss.clear();
  
  string imagePath;
  if (fixPileup) 
    imagePath = "Plots/Mixed/" + particle_lower + pt_str + "pu" + value + "_"; 
  else if (fixLayout) 
    imagePath = "Plots/" + value + "/" + particle_lower + pt_str + "_";

  string legendLabel;
  if (fixPileup) 
    legendLabel = particle_lower + pt_str + "pu" + value + ": ";
  else if (fixLayout) 
    legendLabel = particle_lower + pt_str + " " + value + ": pu";

  string basePath, prePath, postPath;

  if (fixPileup) {
    basePath = "AnalysisResults";
    prePath = basePath + "/";
    postPath = "/dataset" + particle_upper + "_" + 
      pt_str + "GeV/hist-" + particle_lower + 
      pt_str + "pu" + value + ".root"; 
  }
  else if (fixLayout) {
    basePath = "AnalysisResults";
    prePath = basePath + "/" + value;
    postPath = "/dataset" + particle_upper + "_" + 
      pt_str + "GeV/hist-" + particle_lower + 
      pt_str + "pu";  
  }

  vector<string> Layouts;
  Layouts.push_back("ExtBrl4");
  Layouts.push_back("IExtBrl4");
  Layouts.push_back("InclBrl4");

  vector<string> Pileup;
  Pileup.push_back("50");
  Pileup.push_back("140");

  vector<string>* IterVariable = 0;

  if (fixPileup)
    IterVariable = &Layouts;
  else if (fixLayout)
    IterVariable = &Pileup;

  vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kGreen);

  vector<int> markerStyle;
  markerStyle.push_back(1);
  markerStyle.push_back(22);
  markerStyle.push_back(33);

  if (Layouts.size() >= Pileup.size()) {
    Assert("Layouts, colors and markerStyle vectors mismatched", 
	   Layouts.size() == colors.size() &&
	   Layouts.size() == markerStyle.size());
  }
  else {
    Assert("Pileup, colors and markerStyle vectors mismatched", 
	   Pileup.size() == colors.size() &&
	   Pileup.size() == markerStyle.size());
  }

  vector<string> Path;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    if (fixPileup) 
      Path.push_back(prePath + IterVariable->at(i) + postPath);
    else if (fixLayout)
      Path.push_back(prePath + postPath + IterVariable->at(i) + ".root");
  }

  vector<TFile*> inFiles;
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    inFiles.push_back(new TFile(Path[i].c_str()));

  /////// DRAW SIGMA /////
  
  // =============== sigma pt vs eta ================ //

  vector<TGraphErrors*> sigPt_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigPt_abseta")));
    CheckPtr(sigPt_abseta[i]);
  }

  TCanvas* canv_sigPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sigPt_abseta = new TMultiGraph;
  mg_sigPt_abseta   -> SetName("sigPt_abseta_combined");
  mg_sigPt_abseta   -> SetTitle("#sigma(p_{T}) vs |#eta|;Truth |#eta|;#sigma(p_{T}) [MeV]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigPt_abseta[i]     -> SetLineColor(colors[i]);
    sigPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sigPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sigPt_abseta->Add(sigPt_abseta[i]);
  
  mg_sigPt_abseta     -> Draw("APE");
  mg_sigPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //mg_sigPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,4500.0);
  
  TLegend* leg_sigPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sigPt_abseta->AddEntry(sigPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sigPt_abseta  -> SetFillColor(0);
  leg_sigPt_abseta  -> Draw();
  
  canv_sigPt_abseta -> Print((imagePath + "sigPt_abseta.png").c_str(), "png");
  
  // =============== sigma qpt vs eta ================ //

  vector<TGraphErrors*> sigQPt_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigQPt_abseta")));
    CheckPtr(sigQPt_abseta[i]);
  }

  TCanvas* canv_sigQPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sigQPt_abseta = new TMultiGraph;
  mg_sigQPt_abseta   -> SetName("sigQPt_abseta_combined");
  mg_sigQPt_abseta   -> SetTitle("p_{T}^{truth} x #sigma(q/p_{T}) vs |#eta|;Truth |#eta|;p_{T}^{truth} x #sigma(q/p_{T})");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigQPt_abseta[i]     -> SetLineColor(colors[i]);
    sigQPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sigQPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sigQPt_abseta->Add(sigQPt_abseta[i]);
  
  mg_sigQPt_abseta     -> Draw("APE");
  mg_sigQPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigQPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.13);
  
  TLegend* leg_sigQPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sigQPt_abseta->AddEntry(sigQPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sigQPt_abseta  -> SetFillColor(0);
  leg_sigQPt_abseta  -> Draw();
  
  canv_sigQPt_abseta -> Print((imagePath + "sigQPt_abseta.png").c_str(), "png");
  

  // =============== sigma phi vs eta ================ //

  vector<TGraphErrors*> sigPhi_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigPhi_abseta")));
    CheckPtr(sigPhi_abseta[i]);
  }

  TCanvas* canv_sigPhi_abseta = new TCanvas();
  
  TMultiGraph* mg_sigPhi_abseta = new TMultiGraph;
  mg_sigPhi_abseta   -> SetName("sigPhi_abseta_combined");
  mg_sigPhi_abseta   -> SetTitle("#sigma(#phi) vs |#eta|;Truth |#eta|;#sigma(#phi)");

  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigPhi_abseta[i]     -> SetLineColor(colors[i]);
    sigPhi_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sigPhi_abseta[i]     -> SetMarkerColor(colors[i]);
  }
    
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sigPhi_abseta->Add(sigPhi_abseta[i]);
  
  mg_sigPhi_abseta     -> Draw("APE");
  mg_sigPhi_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigPhi_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.0014);
  
  TLegend* leg_sigPhi_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sigPhi_abseta->AddEntry(sigPhi_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sigPhi_abseta  -> SetFillColor(0);
  leg_sigPhi_abseta  -> Draw();
  
  canv_sigPhi_abseta -> Print((imagePath + "sigPhi_abseta.png").c_str(), "png");
  
  // =============== sigma d0 vs eta ================ //

  vector<TGraphErrors*> sigD0_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigD0_abseta")));
    CheckPtr(sigD0_abseta[i]);
  }

  TCanvas* canv_sigD0_abseta = new TCanvas();
  
  TMultiGraph* mg_sigD0_abseta = new TMultiGraph;
  mg_sigD0_abseta   -> SetName("sigD0_abseta_combined");
  mg_sigD0_abseta   -> SetTitle("#sigma(D0) vs |#eta|;Truth |#eta|;#sigma(D0) [mm]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigD0_abseta[i]     -> SetLineColor(colors[i]);
    sigD0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sigD0_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sigD0_abseta->Add(sigD0_abseta[i]);
  
  mg_sigD0_abseta     -> Draw("APE");
  mg_sigD0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigD0_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.11);
  
  TLegend* leg_sigD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sigD0_abseta->AddEntry(sigD0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sigD0_abseta  -> SetFillColor(0);
  leg_sigD0_abseta  -> Draw();
  
  canv_sigD0_abseta -> Print((imagePath + "sigD0_abseta.png").c_str(), "png");
  
  // =============== sigma Z0 vs eta ================ //

  vector<TGraphErrors*> sigZ0_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigZ0_abseta")));
    CheckPtr(sigZ0_abseta[i]);
  }

  TCanvas* canv_sigZ0_abseta = new TCanvas();
  
  TMultiGraph* mg_sigZ0_abseta = new TMultiGraph;
  mg_sigZ0_abseta   -> SetName("sigZ0_abseta_combined");
  mg_sigZ0_abseta   -> SetTitle("#sigma(Z0) vs |#eta|;Truth |#eta|;#sigma(Z0) [mm]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sigZ0_abseta[i]     -> SetLineColor(colors[i]);
    sigZ0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sigZ0_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sigZ0_abseta->Add(sigZ0_abseta[i]);
  
  mg_sigZ0_abseta     -> Draw("APE");
  mg_sigZ0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigZ0_abseta     -> GetYaxis() -> SetRangeUser(0.2,0.7);
  
  TLegend* leg_sigZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sigZ0_abseta->AddEntry(sigZ0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sigZ0_abseta  -> SetFillColor(0);
  leg_sigZ0_abseta  -> Draw();
  
  canv_sigZ0_abseta -> Print((imagePath + "sigZ0_abseta.png").c_str(), "png");
  

  /////// DRAW SYSTEMATIC BIAS /////
  
  // =============== sys pt vs eta ================ //

  vector<TGraphErrors*> sysPt_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysPt_abseta")));
    CheckPtr(sigPt_abseta[i]);
  }

  TCanvas* canv_sysPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sysPt_abseta = new TMultiGraph;
  mg_sysPt_abseta   -> SetName("sysPt_abseta_combined");
  mg_sysPt_abseta   -> SetTitle("bias(p_{T}) vs |#eta|;Truth |#eta|;bias(p_{T}) [MeV]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysPt_abseta[i]     -> SetLineColor(colors[i]);
    sysPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sysPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sysPt_abseta->Add(sysPt_abseta[i]);
  
  mg_sysPt_abseta     -> Draw("APE");
  mg_sysPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysPt_abseta     -> GetYaxis() -> SetRangeUser(-2000.0,400.0);
  
  TLegend* leg_sysPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sysPt_abseta->AddEntry(sysPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sysPt_abseta  -> SetFillColor(0);
  leg_sysPt_abseta  -> Draw();
  
  canv_sysPt_abseta -> Print((imagePath + "sysPt_abseta.png").c_str(), "png");
  
  // =============== sys qpt vs eta ================ //

  vector<TGraphErrors*> sysQPt_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysQPt_abseta")));
    CheckPtr(sigQPt_abseta[i]);
  }

  TCanvas* canv_sysQPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sysQPt_abseta = new TMultiGraph;
  mg_sysQPt_abseta   -> SetName("sysQPt_abseta_combined");
  mg_sysQPt_abseta   -> SetTitle("p_{T}^{truth} x bias(q/p_{T}) vs |#eta|;Truth |#eta|;p_{T}^{truth} x bias(q/p_{T})");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysQPt_abseta[i]     -> SetLineColor(colors[i]);
    sysQPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sysQPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sysQPt_abseta->Add(sysQPt_abseta[i]);
  
  mg_sysQPt_abseta     -> Draw("APE");
  mg_sysQPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysQPt_abseta     -> GetYaxis() -> SetRangeUser(-0.008,0.012);
  
  TLegend* leg_sysQPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sysQPt_abseta->AddEntry(sysQPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sysQPt_abseta  -> SetFillColor(0);
  leg_sysQPt_abseta  -> Draw();
  
  canv_sysQPt_abseta -> Print((imagePath + "sysQPt_abseta.png").c_str(), "png");
  

  // =============== sys phi vs eta ================ //

  vector<TGraphErrors*> sysPhi_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysPhi_abseta")));
    CheckPtr(sigPhi_abseta[i]);
  }

  TCanvas* canv_sysPhi_abseta = new TCanvas();
  
  TMultiGraph* mg_sysPhi_abseta = new TMultiGraph;
  mg_sysPhi_abseta   -> SetName("sysPhi_abseta_combined");
  mg_sysPhi_abseta   -> SetTitle("bias(#phi) vs |#eta|;Truth |#eta|;bias(#phi)");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysPhi_abseta[i]     -> SetLineColor(colors[i]);
    sysPhi_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sysPhi_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sysPhi_abseta->Add(sysPhi_abseta[i]);
  
  mg_sysPhi_abseta     -> Draw("APE");
  mg_sysPhi_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysPhi_abseta     -> GetYaxis() -> SetRangeUser(-0.13e-3,0.06e-3);
  
  TLegend* leg_sysPhi_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sysPhi_abseta->AddEntry(sysPhi_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sysPhi_abseta  -> SetFillColor(0);
  leg_sysPhi_abseta  -> Draw();
  
  canv_sysPhi_abseta -> Print((imagePath + "sysPhi_abseta.png").c_str(), "png");
  
  // =============== sys d0 vs eta ================ //

  vector<TGraphErrors*> sysD0_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysD0_abseta")));
    CheckPtr(sigD0_abseta[i]);
  }

  TCanvas* canv_sysD0_abseta = new TCanvas();
  
  TMultiGraph* mg_sysD0_abseta = new TMultiGraph;
  mg_sysD0_abseta   -> SetName("sysD0_abseta_combined");
  mg_sysD0_abseta   -> SetTitle("bias(D0) vs |#eta|;Truth |#eta|;bias(D0) [mm]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysD0_abseta[i]     -> SetLineColor(colors[i]);
    sysD0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sysD0_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sysD0_abseta->Add(sysD0_abseta[i]);
  
  mg_sysD0_abseta     -> Draw("APE");
  mg_sysD0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysD0_abseta     -> GetYaxis() -> SetRangeUser(-0.015,0.005);
  
  TLegend* leg_sysD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sysD0_abseta->AddEntry(sysD0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sysD0_abseta  -> SetFillColor(0);
  leg_sysD0_abseta  -> Draw();
  
  canv_sysD0_abseta -> Print((imagePath + "sysD0_abseta.png").c_str(), "png");
  
  // =============== sys Z0 vs eta ================ //

  vector<TGraphErrors*> sysZ0_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysZ0_abseta")));
    CheckPtr(sigZ0_abseta[i]);
  }

  TCanvas* canv_sysZ0_abseta = new TCanvas();
  
  TMultiGraph* mg_sysZ0_abseta = new TMultiGraph;
  mg_sysZ0_abseta   -> SetName("sysZ0_abseta_combined");
  mg_sysZ0_abseta   -> SetTitle("bias(Z0) vs |#eta|;Truth |#eta|;bias(Z0) [mm]");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    sysZ0_abseta[i]     -> SetLineColor(colors[i]);
    sysZ0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    sysZ0_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_sysZ0_abseta->Add(sysZ0_abseta[i]);
  
  mg_sysZ0_abseta     -> Draw("APE");
  mg_sysZ0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysZ0_abseta     -> GetYaxis() -> SetRangeUser(-0.25,0.3);
  
  TLegend* leg_sysZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_sysZ0_abseta->AddEntry(sysZ0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_sysZ0_abseta  -> SetFillColor(0);
  leg_sysZ0_abseta  -> Draw();
  
  canv_sysZ0_abseta -> Print((imagePath + "sysZ0_abseta.png").c_str(), "png");

  

  /////// DRAW EFFICIENCY /////
  
  // =============== eff vs eta ================ //
  // we choose classical definition of efficiency, see comments to RunHists.cxx

  vector<TGraphErrors*> eff_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    eff_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_true__eff_abseta")));
    CheckPtr(eff_abseta[i]);
  }

  TCanvas* canv_eff_abseta = new TCanvas();
  
  TMultiGraph* mg_eff_abseta = new TMultiGraph;
  mg_eff_abseta   -> SetName("eff_abseta_combined");
  mg_eff_abseta   -> SetTitle("Efficiency vs |#eta|;Truth |#eta|;Efficiency");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    eff_abseta[i]     -> SetLineColor(colors[i]);
    eff_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    eff_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_eff_abseta->Add(eff_abseta[i]);
  
  mg_eff_abseta     -> Draw("APE");
  mg_eff_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  mg_eff_abseta     -> GetYaxis() -> SetRangeUser(0.92,1.0);
  
  TLegend* leg_eff_abseta = new TLegend(0.55,0.35,0.85,0.15);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_eff_abseta->AddEntry(eff_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_eff_abseta  -> SetFillColor(0);
  leg_eff_abseta  -> Draw();
  
  canv_eff_abseta -> Print((imagePath + "eff_abseta.png").c_str(), "png");

  // =============== eff vs phi ================ //

  vector<TGraphErrors*> eff_phi;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    eff_phi.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_true__eff_phi")));
    CheckPtr(eff_phi[i]);
  }

  TCanvas* canv_eff_phi = new TCanvas();
  
  TMultiGraph* mg_eff_phi = new TMultiGraph;
  mg_eff_phi   -> SetName("eff_phi_combined");
  mg_eff_phi   -> SetTitle("Efficiency vs |#phi|;Truth #phi;Efficiency");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    eff_phi[i]     -> SetLineColor(colors[i]);
    eff_phi[i]     -> SetMarkerStyle(markerStyle[i]);
    eff_phi[i]     -> SetMarkerColor(colors[i]);
  }

  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_eff_phi->Add(eff_phi[i]);
  
  mg_eff_phi     -> Draw("APE");
  mg_eff_phi     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_eff_phi     -> GetYaxis() -> SetRangeUser(0.8,0.92);
  
  TLegend* leg_eff_phi = new TLegend(0.15,0.88,0.35,0.7);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_eff_phi->AddEntry(eff_phi[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_eff_phi  -> SetFillColor(0);
  leg_eff_phi  -> Draw();
  
  canv_eff_phi -> Print((imagePath + "eff_phi.png").c_str(), "png");

  ///////////////////////////// DRAW FAKES /////////////////////////////////


  vector<TGraphErrors*> fakeProb_abseta;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    fakeProb_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHistManager__fakeProb_abseta")));
    CheckPtr(fakeProb_abseta[i]);
  }

  TCanvas* canv_fakeProb_abseta = new TCanvas();
  
  TMultiGraph* mg_fakeProb_abseta = new TMultiGraph;
  mg_fakeProb_abseta   -> SetName("fakeProb_abseta_combined");
  mg_fakeProb_abseta   -> SetTitle("Fake probability vs |#eta|;Truth |#eta|;Fake probability");
  
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    fakeProb_abseta[i]     -> SetLineColor(colors[i]);
    fakeProb_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
    fakeProb_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < IterVariable->size(); i++)
    mg_fakeProb_abseta->Add(fakeProb_abseta[i]);
  
  mg_fakeProb_abseta     -> Draw("APE");
  mg_fakeProb_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  mg_fakeProb_abseta     -> GetYaxis() -> SetRangeUser(0.0, 0.035);
  
  TLegend* leg_fakeProb_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_fakeProb_abseta->AddEntry(fakeProb_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_fakeProb_abseta  -> SetFillColor(0);
  leg_fakeProb_abseta  -> Draw();
  
  canv_fakeProb_abseta -> Print((imagePath + "fakeProb_abseta.png").c_str(), "png");

  ///////////////////////////// DRAW EVENT HISTS /////////////////////////////////

  gStyle->SetOptStat(11);

 
    vector<TH1F*> nRecoTracks;
  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    nRecoTracks.push_back(dynamic_cast<TH1F*>(inFiles[i]->Get("EventHist_all__nRecotrk")));
    CheckPtr(nRecoTracks[i]);
  }

  TCanvas* canv_nRecoTracks = new TCanvas();
  
  if (IterVariable->size() >= 1)
    nRecoTracks[0]     -> SetTitle("Number of reconstructed tracks / event");

  for (unsigned int i = 0; i < IterVariable->size(); i++) {
    nRecoTracks[i]     -> SetLineColor(colors[i]);
    nRecoTracks[i]     -> SetMarkerStyle(markerStyle[i]);
    nRecoTracks[i]     -> SetMarkerColor(colors[i]);
    if (i == 0)
      nRecoTracks[i]     -> DrawCopy();
    else
      nRecoTracks[i]     -> DrawCopy("same");
  }
  
  

  
  TLegend* leg_nRecoTracks = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < IterVariable->size(); i++) 
    leg_nRecoTracks->AddEntry(nRecoTracks[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
  leg_nRecoTracks  -> SetFillColor(0);
  leg_nRecoTracks  -> Draw("same");
  
  canv_nRecoTracks -> Print((imagePath + "nRecoTracks.png").c_str(), "png");

  
}
