#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TLegend.h>
#include <TAxis.h>
#include <string>
#include <sstream>
#include <vector>
#include "utils/checkPtr.h"
#include "utils/assert.h"

using namespace std;

bool debug = true;

void PrintMessage(const std::string& message) {
  if (debug) {
    std::cout << message << std::endl;
  }
}


void fixedPtAndPileup(const string& particle_lowercase, const double& pt_MeV, const string& layout) {

  PrintMessage("Init");

  string particle_lower, particle_upper;
  
  if (particle_lowercase.compare("pi") == 0) {
    particle_lower = particle_lowercase;
    particle_upper = "Pi";
  }

  else if (particle_lowercase.compare("mu") == 0) {
    particle_lower = particle_lowercase;
    particle_upper = "Mu";
  }
    
  string pt_str;

  stringstream ss;
  ss << pt_MeV/1000; //to GeV
  ss >> pt_str;
  ss.clear();

  const string basePath = "AnalysisResults";
  const string prePath = basePath + "/";
  const string postPath = "/dataset" + particle_upper + "_" + 
    pt_str + "GeV/hist-" + particle_lower + 
    pt_str + "pu"; 

  std::vector<string> pileup;
  pileup.push_back("50");
  pileup.push_back("140");

  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);

  std::vector<int> markerStyle;
  markerStyle.push_back(1);
  markerStyle.push_back(33);

  Assert("Pileup, colors and markerStyle not matched",
	 pileup.size() == colors.size() &&
	 pileup.size() == markerStyle.size());

  std::vector<string> Path;
  for (unsigned int i = 0; i < pileup.size(); i++)
    Path.push_back(prePath + layout + postPath + pileup[i] + ".root");

  std::vector<TFile*> inFiles;
  for (unsigned int i = 0; i < pileup.size(); i++) 
    inFiles.push_back(new TFile(Path[i].c_str()));

  /////// DRAW SIGMA /////
  
  // =============== sigma pt vs eta ================ //

  PrintMessage("sigma pt");

  std::vector<TGraphErrors*> sigPt_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sigPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigPt_abseta")));
    CheckPtr(sigPt_abseta[i]);
  }

  PrintMessage("after retrieve");

  TCanvas* canv_sigPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sigPt_abseta = new TMultiGraph;
  mg_sigPt_abseta   -> SetName("sigPt_abseta_combined");
  mg_sigPt_abseta   -> SetTitle("#sigma(p_{T}) vs |#eta|;Truth |#eta|;#sigma(p_{T}) [MeV]");
  
  PrintMessage("after mg decl");

  for (unsigned int i = 0; i < pileup.size(); i++) {
  sigPt_abseta[i]     -> SetLineColor(colors[i]);
  sigPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sigPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }


  PrintMessage("after esthetics");

  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sigPt_abseta->Add(sigPt_abseta[i]);
  
  mg_sigPt_abseta     -> Draw("APE");
  mg_sigPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //mg_sigPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,4500.0);
  

  TLegend* leg_sigPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sigPt_abseta->AddEntry(sigPt_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sigPt_abseta  -> SetFillColor(0);
  leg_sigPt_abseta  -> Draw();
  
  PrintMessage("after leg draw");
  
  canv_sigPt_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sigPt_abseta.png").c_str(), "png");

  
  // =============== sigma qpt vs eta ================ //

  PrintMessage("sigma qpt");

  std::vector<TGraphErrors*> sigQPt_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sigQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigQPt_abseta")));
    CheckPtr(sigQPt_abseta[i]);
  }

  TCanvas* canv_sigQPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sigQPt_abseta = new TMultiGraph;
  mg_sigQPt_abseta   -> SetName("sigQPt_abseta_combined");
  mg_sigQPt_abseta   -> SetTitle("p_{T}^{truth} x #sigma(q/p_{T}) vs |#eta|;Truth |#eta|;p_{T}^{truth} x #sigma(q/p_{T})");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sigQPt_abseta[i]     -> SetLineColor(colors[i]);
  sigQPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sigQPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sigQPt_abseta->Add(sigQPt_abseta[i]);
  
  mg_sigQPt_abseta     -> Draw("APE");
  mg_sigQPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigQPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.13);
  
  TLegend* leg_sigQPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sigQPt_abseta->AddEntry(sigQPt_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sigQPt_abseta  -> SetFillColor(0);
  leg_sigQPt_abseta  -> Draw();
  
  canv_sigQPt_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sigQPt_abseta.png").c_str(), "png");
  

  // =============== sigma phi vs eta ================ //

  std::vector<TGraphErrors*> sigPhi_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sigPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigPhi_abseta")));
    CheckPtr(sigPhi_abseta[i]);
  }

  TCanvas* canv_sigPhi_abseta = new TCanvas();
  
  TMultiGraph* mg_sigPhi_abseta = new TMultiGraph;
  mg_sigPhi_abseta   -> SetName("sigPhi_abseta_combined");
  mg_sigPhi_abseta   -> SetTitle("#sigma(#phi) vs |#eta|;Truth |#eta|;#sigma(#phi)");
 
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sigPhi_abseta[i]     -> SetLineColor(colors[i]);
  sigPhi_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sigPhi_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sigPhi_abseta->Add(sigPhi_abseta[i]);
  
  mg_sigPhi_abseta     -> Draw("APE");
  mg_sigPhi_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigPhi_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.0014);
  
  TLegend* leg_sigPhi_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sigPhi_abseta->AddEntry(sigPhi_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sigPhi_abseta  -> SetFillColor(0);
  leg_sigPhi_abseta  -> Draw();
  
  canv_sigPhi_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sigPhi_abseta.png").c_str(), "png");
  
  // =============== sigma d0 vs eta ================ //

  std::vector<TGraphErrors*> sigD0_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sigD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigD0_abseta")));
    CheckPtr(sigD0_abseta[i]);
  }

  TCanvas* canv_sigD0_abseta = new TCanvas();
  
  TMultiGraph* mg_sigD0_abseta = new TMultiGraph;
  mg_sigD0_abseta   -> SetName("sigD0_abseta_combined");
  mg_sigD0_abseta   -> SetTitle("#sigma(D0) vs |#eta|;Truth |#eta|;#sigma(D0) [mm]");
    
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sigD0_abseta[i]     -> SetLineColor(colors[i]);
  sigD0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sigD0_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sigD0_abseta->Add(sigD0_abseta[i]);
  
  mg_sigD0_abseta     -> Draw("APE");
  mg_sigD0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigD0_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.11);
  
  TLegend* leg_sigD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sigD0_abseta->AddEntry(sigD0_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sigD0_abseta  -> SetFillColor(0);
  leg_sigD0_abseta  -> Draw();
  
  canv_sigD0_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sigD0_abseta.png").c_str(), "png");
  
  // =============== sigma Z0 vs eta ================ //

  std::vector<TGraphErrors*> sigZ0_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sigZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sigZ0_abseta")));
    CheckPtr(sigZ0_abseta[i]);
  }

  TCanvas* canv_sigZ0_abseta = new TCanvas();
  
  TMultiGraph* mg_sigZ0_abseta = new TMultiGraph;
  mg_sigZ0_abseta   -> SetName("sigZ0_abseta_combined");
  mg_sigZ0_abseta   -> SetTitle("#sigma(Z0) vs |#eta|;Truth |#eta|;#sigma(Z0) [mm]");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sigZ0_abseta[i]     -> SetLineColor(colors[i]);
  sigZ0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sigZ0_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sigZ0_abseta->Add(sigZ0_abseta[i]);
  
  mg_sigZ0_abseta     -> Draw("APE");
  mg_sigZ0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sigZ0_abseta     -> GetYaxis() -> SetRangeUser(0.2,0.7);
  
  TLegend* leg_sigZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sigZ0_abseta->AddEntry(sigZ0_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sigZ0_abseta  -> SetFillColor(0);
  leg_sigZ0_abseta  -> Draw();
  
  canv_sigZ0_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sigZ0_abseta.png").c_str(), "png");
  

  /////// DRAW SYSTEMATIC BIAS /////
  
  // =============== sys pt vs eta ================ //

  std::vector<TGraphErrors*> sysPt_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sysPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysPt_abseta")));
    CheckPtr(sigPt_abseta[i]);
  }

  TCanvas* canv_sysPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sysPt_abseta = new TMultiGraph;
  mg_sysPt_abseta   -> SetName("sysPt_abseta_combined");
  mg_sysPt_abseta   -> SetTitle("bias(p_{T}) vs |#eta|;Truth |#eta|;bias(p_{T}) [MeV]");

  for (unsigned int i = 0; i < pileup.size(); i++) {
  sysPt_abseta[i]     -> SetLineColor(colors[i]);
  sysPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sysPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sysPt_abseta->Add(sysPt_abseta[i]);
  
  mg_sysPt_abseta     -> Draw("APE");
  mg_sysPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysPt_abseta     -> GetYaxis() -> SetRangeUser(-2000.0,400.0);
  
  TLegend* leg_sysPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sysPt_abseta->AddEntry(sysPt_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sysPt_abseta  -> SetFillColor(0);
  leg_sysPt_abseta  -> Draw();
  
  canv_sysPt_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sysPt_abseta.png").c_str(), "png");
  
  // =============== sys qpt vs eta ================ //

  std::vector<TGraphErrors*> sysQPt_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sysQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysQPt_abseta")));
    CheckPtr(sigQPt_abseta[i]);
  }

  TCanvas* canv_sysQPt_abseta = new TCanvas();
  
  TMultiGraph* mg_sysQPt_abseta = new TMultiGraph;
  mg_sysQPt_abseta   -> SetName("sysQPt_abseta_combined");
  mg_sysQPt_abseta   -> SetTitle("p_{T}^{truth} x bias(q/p_{T}) vs |#eta|;Truth |#eta|;p_{T}^{truth} x bias(q/p_{T})");

  for (unsigned int i = 0; i < pileup.size(); i++) {
  sysQPt_abseta[i]     -> SetLineColor(colors[i]);
  sysQPt_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sysQPt_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sysQPt_abseta->Add(sysQPt_abseta[i]);
  
  mg_sysQPt_abseta     -> Draw("APE");
  mg_sysQPt_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysQPt_abseta     -> GetYaxis() -> SetRangeUser(-0.008,0.012);
  
  TLegend* leg_sysQPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sysQPt_abseta->AddEntry(sysQPt_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sysQPt_abseta  -> SetFillColor(0);
  leg_sysQPt_abseta  -> Draw();
  
  canv_sysQPt_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sysQPt_abseta.png").c_str(), "png");
  

  // =============== sys phi vs eta ================ //

  std::vector<TGraphErrors*> sysPhi_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sysPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysPhi_abseta")));
    CheckPtr(sigPhi_abseta[i]);
  }

  TCanvas* canv_sysPhi_abseta = new TCanvas();
  
  TMultiGraph* mg_sysPhi_abseta = new TMultiGraph;
  mg_sysPhi_abseta   -> SetName("sysPhi_abseta_combined");
  mg_sysPhi_abseta   -> SetTitle("bias(#phi) vs |#eta|;Truth |#eta|;bias(#phi)");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sysPhi_abseta[i]     -> SetLineColor(colors[i]);
  sysPhi_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sysPhi_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sysPhi_abseta->Add(sysPhi_abseta[i]);
  
  mg_sysPhi_abseta     -> Draw("APE");
  mg_sysPhi_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysPhi_abseta     -> GetYaxis() -> SetRangeUser(-0.13e-3,0.06e-3);
  
  TLegend* leg_sysPhi_abseta = new TLegend(0.15,0.35,0.35,0.15);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sysPhi_abseta->AddEntry(sysPhi_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sysPhi_abseta  -> SetFillColor(0);
  leg_sysPhi_abseta  -> Draw();
  
  canv_sysPhi_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sysPhi_abseta.png").c_str(), "png");
  
  // =============== sys d0 vs eta ================ //

  std::vector<TGraphErrors*> sysD0_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sysD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysD0_abseta")));
    CheckPtr(sigD0_abseta[i]);
  }

  TCanvas* canv_sysD0_abseta = new TCanvas();
  
  TMultiGraph* mg_sysD0_abseta = new TMultiGraph;
  mg_sysD0_abseta   -> SetName("sysD0_abseta_combined");
  mg_sysD0_abseta   -> SetTitle("bias(D0) vs |#eta|;Truth |#eta|;bias(D0) [mm]");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sysD0_abseta[i]     -> SetLineColor(colors[i]);
  sysD0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sysD0_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sysD0_abseta->Add(sysD0_abseta[i]);
  
  mg_sysD0_abseta     -> Draw("APE");
  mg_sysD0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysD0_abseta     -> GetYaxis() -> SetRangeUser(-0.015,0.005);
  
  TLegend* leg_sysD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sysD0_abseta->AddEntry(sysD0_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sysD0_abseta  -> SetFillColor(0);
  leg_sysD0_abseta  -> Draw();
  
  canv_sysD0_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sysD0_abseta.png").c_str(), "png");
  
  // =============== sys Z0 vs eta ================ //

  std::vector<TGraphErrors*> sysZ0_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    sysZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__sysZ0_abseta")));
    CheckPtr(sigZ0_abseta[i]);
  }

  TCanvas* canv_sysZ0_abseta = new TCanvas();
  
  TMultiGraph* mg_sysZ0_abseta = new TMultiGraph;
  mg_sysZ0_abseta   -> SetName("sysZ0_abseta_combined");
  mg_sysZ0_abseta   -> SetTitle("bias(Z0) vs |#eta|;Truth |#eta|;bias(Z0) [mm]");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  sysZ0_abseta[i]     -> SetLineColor(colors[i]);
  sysZ0_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  sysZ0_abseta[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_sysZ0_abseta->Add(sysZ0_abseta[i]);
  
  mg_sysZ0_abseta     -> Draw("APE");
  mg_sysZ0_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_sysZ0_abseta     -> GetYaxis() -> SetRangeUser(-0.25,0.3);
  
  TLegend* leg_sysZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_sysZ0_abseta->AddEntry(sysZ0_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_sysZ0_abseta  -> SetFillColor(0);
  leg_sysZ0_abseta  -> Draw();
  
  canv_sysZ0_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "sysZ0_abseta.png").c_str(), "png");

  

  /////// DRAW EFFICIENCY /////
  
  // =============== eff vs eta ================ //

  std::vector<TGraphErrors*> eff_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    eff_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__eff_abseta")));
    CheckPtr(eff_abseta[i]);
  }

  TCanvas* canv_eff_abseta = new TCanvas();
  
  TMultiGraph* mg_eff_abseta = new TMultiGraph;
  mg_eff_abseta   -> SetName("eff_abseta_combined");
  mg_eff_abseta   -> SetTitle("Efficiency vs |#eta|;Truth |#eta|;Efficiency");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  eff_abseta[i]     -> SetLineColor(colors[i]);
  eff_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  eff_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_eff_abseta->Add(eff_abseta[i]);
  
  mg_eff_abseta     -> Draw("APE");
  mg_eff_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  mg_eff_abseta     -> GetYaxis() -> SetRangeUser(0.92,1.0);
  
  TLegend* leg_eff_abseta = new TLegend(0.55,0.35,0.85,0.15);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_eff_abseta->AddEntry(eff_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_eff_abseta  -> SetFillColor(0);
  leg_eff_abseta  -> Draw();
  
  canv_eff_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "eff_abseta.png").c_str(), "png");

  // =============== eff vs phi ================ //

  std::vector<TGraphErrors*> eff_phi;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    eff_phi.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__eff_phi")));
    CheckPtr(eff_phi[i]);
  }

  TCanvas* canv_eff_phi = new TCanvas();
  
  TMultiGraph* mg_eff_phi = new TMultiGraph;
  mg_eff_phi   -> SetName("eff_phi_combined");
  mg_eff_phi   -> SetTitle("Efficiency vs |#phi|;Truth #phi;Efficiency");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  eff_phi[i]     -> SetLineColor(colors[i]);
  eff_phi[i]     -> SetMarkerStyle(markerStyle[i]);
  eff_phi[i]     -> SetMarkerColor(colors[i]);
  }

  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_eff_phi->Add(eff_phi[i]);
  
  mg_eff_phi     -> Draw("APE");
  mg_eff_phi     -> GetYaxis() -> SetTitleOffset(1.2);
  //  mg_eff_phi     -> GetYaxis() -> SetRangeUser(0.8,0.92);
  
  TLegend* leg_eff_phi = new TLegend(0.15,0.88,0.35,0.7);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_eff_phi->AddEntry(eff_phi[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_eff_phi  -> SetFillColor(0);
  leg_eff_phi  -> Draw();
  
  canv_eff_phi -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "eff_phi.png").c_str(), "png");

  ///////////////////////////// DRAW FAKES /////////////////////////////////


  std::vector<TGraphErrors*> fakeProb_abseta;
  for (unsigned int i = 0; i < pileup.size(); i++) {
    fakeProb_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHistManager__fakeProb_abseta")));
    CheckPtr(fakeProb_abseta[i]);
  }

  TCanvas* canv_fakeProb_abseta = new TCanvas();
  
  TMultiGraph* mg_fakeProb_abseta = new TMultiGraph;
  mg_fakeProb_abseta   -> SetName("fakeProb_abseta_combined");
  mg_fakeProb_abseta   -> SetTitle("Fake probability vs |#eta|;Truth |#eta|;Fake probability");
  
  for (unsigned int i = 0; i < pileup.size(); i++) {
  fakeProb_abseta[i]     -> SetLineColor(colors[i]);
  fakeProb_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
  fakeProb_abseta[i]     -> SetMarkerColor(colors[i]);
  }
  
  for (unsigned int i = 0; i < pileup.size(); i++)
    mg_fakeProb_abseta->Add(fakeProb_abseta[i]);
  
  mg_fakeProb_abseta     -> Draw("APE");
  mg_fakeProb_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  mg_fakeProb_abseta     -> GetYaxis() -> SetRangeUser(0.0, 0.035);
  
  TLegend* leg_fakeProb_abseta = new TLegend(0.15,0.85,0.35,0.55);
  for (unsigned int i = 0; i < pileup.size(); i++) 
    leg_fakeProb_abseta->AddEntry(fakeProb_abseta[i], (particle_lower + pt_str + "pu" + pileup[i] + ": " + layout).c_str(), "lp");
  
  leg_fakeProb_abseta  -> SetFillColor(0);
  leg_fakeProb_abseta  -> Draw();
  
  canv_fakeProb_abseta -> Print(("Plots/" + layout + "/" + particle_lower + pt_str + "fakeProb_abseta.png").c_str(), "png");

  
}
