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
using namespace std;

void fixedPtAndPileup(const string& particle, const double& pt, const int& pileup) {

  string particle_lower, particle_upper;
  
  if (particle.compare("pi") == 0) {
    particle_lower = particle;
    particle_upper = "Pi";
  }

  else if (particle.compare("mu") == 0) {
    particle_lower = particle;
    particle_upper = "Mu";
  }


    
  string pt_str, pileup_str;

  stringstream ss;
  ss << pt/1000; //to GeV
  ss >> pt_str;
  ss.clear();

  ss << pileup;
  ss >> pileup_str;
  ss.clear();

  const string basePath = "AnalysisResults";
  const string prePath = basePath + "/";
  const string postPath = "/dataset" + particle_upper + "_" + 
    pt_str + "GeV/hist-" + particle_lower + 
    pt_str + "pu" + pileup_str + ".root"; 

  std::vector<string> Layouts;
  Layouts.push_back("ExtBrl4");
  Layouts.push_back("IExtBrl4");
  Layouts.push_back("InclBrl4");
  
  std::vector<string> Path;
  for (unsigned int i = 0; i < Layouts.size(); i++)
    Path.push_back(prePath + Layouts[i] + postPath);

  std::vector<TFile*> inFiles;
  for (unsigned int i = 0; i < Layouts.size(); i++) 
    inFiles.push_back(new TFile(Path[i].c_str()));

  
  /////// DRAW EFFICIENCY /////
  
  std::vector<TGraphErrors*> eff_abseta;
  for (unsigned int i = 0; i < Layouts.size(); i++) {
    eff_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all__eff_abseta")));
    CheckPtr(eff_abseta[i]);
  }

  TCanvas* canv_eff_abseta = new TCanvas();
  
  TMultiGraph* mg_eff_abseta = new TMultiGraph;
  mg_eff_abseta   -> SetName("eff_abseta_combined");
  mg_eff_abseta   -> SetTitle("Efficiency vs |#eta|;Truth |#eta|;Efficiency");
  

  eff_abseta[0]     -> SetLineColor(kRed);
  eff_abseta[1]     -> SetLineColor(kBlue);
  eff_abseta[2]     -> SetLineColor(kGreen);

  eff_abseta[0]     -> SetMarkerStyle(1);
  eff_abseta[1]     -> SetMarkerStyle(22);
  eff_abseta[2]     -> SetMarkerStyle(33);

  eff_abseta[0]     -> SetMarkerColor(kRed);
  eff_abseta[1]     -> SetMarkerColor(kBlue);
  eff_abseta[2]     -> SetMarkerColor(kGreen);
  
  for (unsigned int i = 0; i < Layouts.size(); i++)
    mg_eff_abseta->Add(eff_abseta[i]);
  
  mg_eff_abseta     -> Draw("APE");
  mg_eff_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
  mg_eff_abseta     -> GetYaxis() -> SetRangeUser(0.92,1.0);
  
  TLegend* leg_eff_abseta = new TLegend(0.55,0.35,0.85,0.15);
  for (unsigned int i = 0; i < Layouts.size(); i++) 
    leg_eff_abseta->AddEntry(eff_abseta[i], (particle_lower + pt_str + "pu" + pileup_str + ": " + Layouts[i]).c_str(), "lp");
  
  leg_eff_abseta  -> SetFillColor(0);
  leg_eff_abseta  -> Draw();
  
  canv_eff_abseta -> Print(("Plots/" + particle_lower + pt_str + "pu" + pileup_str + "eff_abseta.png").c_str(), "png");


}
