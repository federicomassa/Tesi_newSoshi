#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
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

void fixedPt(const string& particle, const double& pt, const string& fixOption, const string& value,
	     const bool& doSigmaPlots = true, 
	     const bool& doBiasPlots = true, 
	     const bool& doEffPlots = true, 
	     const bool& doFakePlots = true, 
	     const bool& doEventPlots = true,
	     const bool& doHitsPlots = true) {

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


  if (fixOption.compare("PILEUP") == 0) {
    fixPileup = true;
    fixLayout = false;
  }
  else if (fixOption.compare("LAYOUT") == 0) {
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
  //  Pileup.push_back("10");
  Pileup.push_back("50");
  Pileup.push_back("140");
  Pileup.push_back("200");
  //Pileup.push_back("300");

  vector<string>* IterVariable = 0;

  if (fixPileup)
    IterVariable = &Layouts;
  else if (fixLayout)
    IterVariable = &Pileup;

  vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kGreen);
  colors.push_back(kBlack);
  colors.push_back(kPink - 8);

  vector<int> markerStyle;
  markerStyle.push_back(1);
  markerStyle.push_back(22);
  markerStyle.push_back(33);
  markerStyle.push_back(34);
  markerStyle.push_back(35);

  if (Layouts.size() >= Pileup.size()) {
    Assert("Layouts, colors and markerStyle vectors mismatched", 
	   Layouts.size() <= colors.size() &&
	   Layouts.size() <= markerStyle.size());
  }
  else {
    Assert("Pileup, colors and markerStyle vectors mismatched", 
	   Pileup.size() <= colors.size() &&
	   Pileup.size() <= markerStyle.size());
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


  if (doSigmaPlots) {
    /////// DRAW SIGMA /////
  
    // =============== sigma pt vs eta ================ //

    vector<TGraphErrors*> sigPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigPt_abseta")));
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
    //    mg_sigPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,4500.0);
    mg_sigPt_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
  
    TLegend* leg_sigPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigPt_abseta->AddEntry(sigPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigPt_abseta  -> SetFillColor(0);
    leg_sigPt_abseta  -> Draw();
  
    canv_sigPt_abseta -> Print((imagePath + "sigPt_abseta.png").c_str(), "png");
    //    delete canv_sigPt_abseta;
  
    // =============== sigma qpt vs eta ================ //

    vector<TGraphErrors*> sigQPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigQPt_abseta")));
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
    mg_sigQPt_abseta     -> GetYaxis() -> SetRangeUser(0.0,1.0);
    mg_sigQPt_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
  
    TLegend* leg_sigQPt_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigQPt_abseta->AddEntry(sigQPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigQPt_abseta  -> SetFillColor(0);
    leg_sigQPt_abseta  -> Draw();
  
    canv_sigQPt_abseta -> Print((imagePath + "sigQPt_abseta.png").c_str(), "png");
    //    delete canv_sigQPt_abseta;
  
    // =============== sigma eta vs eta ================ //

    vector<TGraphErrors*> sigEta_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigEta_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigEta_abseta")));
      CheckPtr(sigEta_abseta[i]);
    }

    TCanvas* canv_sigEta_abseta = new TCanvas();
  
    TMultiGraph* mg_sigEta_abseta = new TMultiGraph;
    mg_sigEta_abseta   -> SetName("sigEta_abseta_combined");
    mg_sigEta_abseta   -> SetTitle("#sigma(#eta)/|#eta_{truth}| vs |#eta|;Truth |#eta|;#sigma(#eta)/|#eta_{truth}|");
  
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigEta_abseta[i]     -> SetLineColor(colors[i]);
      sigEta_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
      sigEta_abseta[i]     -> SetMarkerColor(colors[i]);
    }

    for (unsigned int i = 0; i < IterVariable->size(); i++)
      mg_sigEta_abseta->Add(sigEta_abseta[i]);
  
    mg_sigEta_abseta     -> Draw("APE");
    mg_sigEta_abseta     -> GetYaxis() -> SetTitleOffset(1.4);
    mg_sigEta_abseta     -> GetYaxis() -> SetRangeUser(0.0002,0.002);
    mg_sigEta_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
  
    TLegend* leg_sigEta_abseta = new TLegend(0.68,0.88,0.88,0.68);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigEta_abseta->AddEntry(sigEta_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigEta_abseta  -> SetFillColor(0);
    leg_sigEta_abseta  -> Draw();
  
    canv_sigEta_abseta -> Print((imagePath + "sigEta_abseta.png").c_str(), "png");
    //    delete canv_sigEta_abseta;
  


    // =============== sigma phi vs eta ================ //

    vector<TGraphErrors*> sigPhi_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigPhi_abseta")));
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
    mg_sigPhi_abseta     -> GetYaxis() -> SetRangeUser(0.0,0.0026);
    mg_sigPhi_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);

    TLegend* leg_sigPhi_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigPhi_abseta->AddEntry(sigPhi_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigPhi_abseta  -> SetFillColor(0);
    leg_sigPhi_abseta  -> Draw();
  
    canv_sigPhi_abseta -> Print((imagePath + "sigPhi_abseta.png").c_str(), "png");
    //    delete canv_sigPhi_abseta;
  
    // =============== sigma d0 vs eta ================ //

    vector<TGraphErrors*> sigD0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigD0_abseta")));
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
    mg_sigD0_abseta     -> GetYaxis() -> SetRangeUser(0.01,0.09);
    mg_sigD0_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
  
    TLegend* leg_sigD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigD0_abseta->AddEntry(sigD0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigD0_abseta  -> SetFillColor(0);
    leg_sigD0_abseta  -> Draw();
  
    canv_sigD0_abseta -> Print((imagePath + "sigD0_abseta.png").c_str(), "png");
    //    delete canv_sigD0_abseta;
  
    // =============== sigma Z0 vs eta ================ //

    vector<TGraphErrors*> sigZ0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigZ0_abseta")));
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
    mg_sigZ0_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_sigZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sigZ0_abseta->AddEntry(sigZ0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sigZ0_abseta  -> SetFillColor(0);
    leg_sigZ0_abseta  -> Draw();
  
    canv_sigZ0_abseta -> Print((imagePath + "sigZ0_abseta.png").c_str(), "png");
    //    delete canv_sigZ0_abseta;
  }  

  if (doBiasPlots) {
    /////// DRAW SYSTEMATIC BIAS /////
  
    // =============== sys pt vs eta ================ //

    vector<TGraphErrors*> sysPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysPt_abseta")));
      CheckPtr(sysPt_abseta[i]);
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
    mg_sysPt_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_sysPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysPt_abseta->AddEntry(sysPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysPt_abseta  -> SetFillColor(0);
    leg_sysPt_abseta  -> Draw();
  
    canv_sysPt_abseta -> Print((imagePath + "sysPt_abseta.png").c_str(), "png");
    //    delete canv_sysPt_abseta;
  
    // =============== sys qpt vs eta ================ //

    vector<TGraphErrors*> sysQPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysQPt_abseta")));
      CheckPtr(sysQPt_abseta[i]);
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
    mg_sysQPt_abseta     -> GetYaxis() -> SetRangeUser(-0.1,0.1);
    mg_sysQPt_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_sysQPt_abseta = new TLegend(0.15,0.35,0.35,0.15);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysQPt_abseta->AddEntry(sysQPt_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysQPt_abseta  -> SetFillColor(0);
    leg_sysQPt_abseta  -> Draw();
  
    canv_sysQPt_abseta -> Print((imagePath + "sysQPt_abseta.png").c_str(), "png");
    //    delete canv_sysQPt_abseta;
  
    // =============== sys eta vs eta ================ //

    vector<TGraphErrors*> sysEta_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysEta_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysEta_abseta")));
      CheckPtr(sysEta_abseta[i]);
    }

    TCanvas* canv_sysEta_abseta = new TCanvas();
  
    TMultiGraph* mg_sysEta_abseta = new TMultiGraph;
    mg_sysEta_abseta   -> SetName("sysEta_abseta_combined");
    mg_sysEta_abseta   -> SetTitle("bias(#eta)/|#eta_{truth}| vs |#eta|;Truth |#eta|;bias(#eta)/|#eta_{truth}");
  
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysEta_abseta[i]     -> SetLineColor(colors[i]);
      sysEta_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
      sysEta_abseta[i]     -> SetMarkerColor(colors[i]);
    }

  
    for (unsigned int i = 0; i < IterVariable->size(); i++)
      mg_sysEta_abseta->Add(sysEta_abseta[i]);
  
    mg_sysEta_abseta     -> Draw("APE");
    mg_sysEta_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
    mg_sysEta_abseta     -> GetYaxis() -> SetRangeUser(-0.05E-3,0.15E-3);
    mg_sysEta_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_sysEta_abseta = new TLegend(0.68,0.88,0.88,0.68);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysEta_abseta->AddEntry(sysEta_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysEta_abseta  -> SetFillColor(0);
    leg_sysEta_abseta  -> Draw();
  
    canv_sysEta_abseta -> Print((imagePath + "sysEta_abseta.png").c_str(), "png");
    //    delete canv_sysEta_abseta;
  


    // =============== sys phi vs eta ================ //

    vector<TGraphErrors*> sysPhi_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysPhi_abseta")));
      CheckPtr(sysPhi_abseta[i]);
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
    mg_sysPhi_abseta     -> GetYaxis() -> SetRangeUser(-0.15e-3,0.25e-3);
    mg_sysPhi_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);

    TLegend* leg_sysPhi_abseta = new TLegend(0.15,0.88,0.31,0.68);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysPhi_abseta->AddEntry(sysPhi_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysPhi_abseta  -> SetFillColor(0);
    leg_sysPhi_abseta  -> Draw();
  
    canv_sysPhi_abseta -> Print((imagePath + "sysPhi_abseta.png").c_str(), "png");
    //    delete canv_sysPhi_abseta;
  
    // =============== sys d0 vs eta ================ //

    vector<TGraphErrors*> sysD0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysD0_abseta")));
      CheckPtr(sysD0_abseta[i]);
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
    mg_sysD0_abseta     -> GetYaxis() -> SetRangeUser(-0.015,0.006);
    mg_sysD0_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_sysD0_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysD0_abseta->AddEntry(sysD0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysD0_abseta  -> SetFillColor(0);
    leg_sysD0_abseta  -> Draw();
  
    canv_sysD0_abseta -> Print((imagePath + "sysD0_abseta.png").c_str(), "png");
    //    delete canv_sysD0_abseta;
  
    // =============== sys Z0 vs eta ================ //

    vector<TGraphErrors*> sysZ0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysZ0_abseta")));
      CheckPtr(sysZ0_abseta[i]);
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
    mg_sysZ0_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
      
    TLegend* leg_sysZ0_abseta = new TLegend(0.7,0.28,0.88,0.12);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_sysZ0_abseta->AddEntry(sysZ0_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_sysZ0_abseta  -> SetFillColor(0);
    leg_sysZ0_abseta  -> Draw();
  
    canv_sysZ0_abseta -> Print((imagePath + "sysZ0_abseta.png").c_str(), "png");
    //    delete canv_sysZ0_abseta;

  }

  if (doEffPlots) {
    /////// DRAW EFFICIENCY /////
  
    // =============== eff vs eta ================ //
    // we choose classical definition of efficiency, see comments to RunHists.cxx

    vector<TGraphErrors*> eff_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      eff_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_true_hard_all__eff_abseta")));
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
    mg_eff_abseta     -> GetYaxis() -> SetRangeUser(0.85,1.0);
    mg_eff_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_eff_abseta = new TLegend(0.15,0.30,0.30,0.15);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_eff_abseta->AddEntry(eff_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_eff_abseta  -> SetFillColor(0);
    leg_eff_abseta  -> Draw();
  
    canv_eff_abseta -> Print((imagePath + "eff_abseta.png").c_str(), "png");
    //    delete canv_eff_abseta;

    // =============== eff vs phi ================ //

    vector<TGraphErrors*> eff_phi;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      eff_phi.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_true_hard_all__eff_phi")));
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
    mg_eff_phi     -> GetXaxis() -> SetRangeUser(-3.2,3.2);  

    TLegend* leg_eff_phi = new TLegend(0.80,0.98,0.98,0.78);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_eff_phi->AddEntry(eff_phi[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_eff_phi  -> SetFillColor(0);
    leg_eff_phi  -> Draw();
  
    canv_eff_phi -> Print((imagePath + "eff_phi.png").c_str(), "png");
    //    delete canv_eff_phi;

  }

  if (doFakePlots) {
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
    mg_fakeProb_abseta     -> GetYaxis() -> SetRangeUser(0.0, 0.008);
    mg_fakeProb_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_fakeProb_abseta = new TLegend(0.15,0.85,0.35,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_fakeProb_abseta->AddEntry(fakeProb_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_fakeProb_abseta  -> SetFillColor(0);
    leg_fakeProb_abseta  -> Draw();
  
    canv_fakeProb_abseta -> Print((imagePath + "fakeProb_abseta.png").c_str(), "png");
    //    delete canv_fakeProb_abseta;
  }

  if (doEventPlots) {
    ///////////////////////////// DRAW EVENT HISTS /////////////////////////////////

    gStyle->SetOptStat(11);

 
    vector<TH1F*> nRecoTracks;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nRecoTracks.push_back(dynamic_cast<TH1F*>(inFiles[i]->Get("EventHist_all__nRecotrkwide")));
      CheckPtr(nRecoTracks[i]);
    }

    TCanvas* canv_nRecoTracks = new TCanvas();
  
    if (IterVariable->size() >= 1)
      nRecoTracks[0]     -> SetTitle("Number of reconstructed tracks / event");

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nRecoTracks[i]     -> SetLineColor(colors[i]);
      nRecoTracks[i]     -> SetMarkerStyle(markerStyle[i]);
      nRecoTracks[i]     -> SetMarkerColor(colors[i]);
      nRecoTracks[i]     -> GetXaxis() -> SetRangeUser(-0.5,10.5);
      if (i == 0) {
	nRecoTracks[i]     -> DrawCopy();
      }
      else
	nRecoTracks[i]     -> DrawCopy("same");
    }
  
    
    TLegend* leg_nRecoTracks = new TLegend(0.55,0.85,0.75,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_nRecoTracks->AddEntry(nRecoTracks[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_nRecoTracks  -> SetFillColor(0);
    leg_nRecoTracks  -> Draw("same");
    
    canv_nRecoTracks -> Print((imagePath + "nRecoTracks.png").c_str(), "png");
    //    delete canv_nRecoTracks;
  
    ////////////////// N charged truth /////////////////
    // TODO: understand why NchargedTruth does NOT include muons (|pdgId| != 13)

    vector<TH1F*> nPrimaryChargedTruth;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nPrimaryChargedTruth.push_back(dynamic_cast<TH1F*>(inFiles[i]->Get("EventHist_all__nPrimaryChargedTruth")));
      CheckPtr(nPrimaryChargedTruth[i]);
    }

    TCanvas* canv_nPrimaryChargedTruth = new TCanvas();
  
    if (IterVariable->size() >= 1)
      nPrimaryChargedTruth[0]     -> SetTitle("Number of primary charged truth in cone / event");

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nPrimaryChargedTruth[i]     -> SetLineColor(colors[i]);
      nPrimaryChargedTruth[i]     -> SetMarkerStyle(markerStyle[i]);
      nPrimaryChargedTruth[i]     -> SetMarkerColor(colors[i]);
      nPrimaryChargedTruth[i]     -> GetXaxis() -> SetRangeUser(-0.5,29.5);
      if (i == 0) {
	nPrimaryChargedTruth[i]     -> DrawCopy();
      }
      else
	nPrimaryChargedTruth[i]     -> DrawCopy("same");
    }
  
    
    TLegend* leg_nPrimaryChargedTruth = new TLegend(0.55,0.85,0.75,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_nPrimaryChargedTruth->AddEntry(nPrimaryChargedTruth[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_nPrimaryChargedTruth  -> SetFillColor(0);
    leg_nPrimaryChargedTruth  -> Draw("same");
    
    canv_nPrimaryChargedTruth -> Print((imagePath + "nPrimaryChargedTruth.png").c_str(), "png");
    //    delete canv_nPrimaryChargedTruth;
  
    ////////////////// N charged truth /////////////////

    vector<TH1F*> nChargedTruth;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nChargedTruth.push_back(dynamic_cast<TH1F*>(inFiles[i]->Get("EventHist_all__nChargedTruth")));
      CheckPtr(nChargedTruth[i]);
    }

    TCanvas* canv_nChargedTruth = new TCanvas();
  
    if (IterVariable->size() >= 1)
      nChargedTruth[0]     -> SetTitle("Number of charged truth in cone / event");

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nChargedTruth[i]     -> SetLineColor(colors[i]);
      nChargedTruth[i]     -> SetMarkerStyle(markerStyle[i]);
      nChargedTruth[i]     -> SetMarkerColor(colors[i]);
      nChargedTruth[i]     -> GetXaxis() -> SetRangeUser(-0.5,59.5);
      if (i == 0) {
	nChargedTruth[i]     -> DrawCopy();
      }
      else
	nChargedTruth[i]     -> DrawCopy("same");
    }
  
    
    TLegend* leg_nChargedTruth = new TLegend(0.55,0.85,0.75,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_nChargedTruth->AddEntry(nChargedTruth[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_nChargedTruth  -> SetFillColor(0);
    leg_nChargedTruth  -> Draw("same");
    
    canv_nChargedTruth -> Print((imagePath + "nChargedTruth.png").c_str(), "png");
    //    delete canv_nChargedTruth;
  


  }

  if (doHitsPlots) {

    gStyle->SetOptStat(1111);

    // PixHits vs eta
    vector<TH2F*> nPixHits;
    vector<TCanvas*> canv_nPixHits;

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nPixHits.push_back(dynamic_cast<TH2F*>(inFiles[i]->Get("TrackHist_reco_all__eta_nPixHits")));
      CheckPtr(nPixHits[i]);
      
      if (fixPileup)
	nPixHits[i] -> SetTitle(("Number of Pixel hits vs |#eta|: " + IterVariable->at(i) + " " + particle_lower + pt_str + "pu" + value + 
				";|#eta|;").c_str());
      else if (fixLayout)
	nPixHits[i] -> SetTitle(("Number of Pixel hits vs |#eta|: " + value + " " + particle_lower + pt_str + "pu" + IterVariable->at(i) + 
				 ";|#eta|;").c_str());
      	
      canv_nPixHits.push_back(new TCanvas());
      canv_nPixHits[i]->cd();
      nPixHits[i]     -> GetXaxis() -> SetRangeUser(0.0,4.0);
      nPixHits[i]     -> DrawCopy("COL0");

      string hitsPlotPath;
      
      if (fixPileup) 
	hitsPlotPath = "Plots/" + IterVariable->at(i) + "/" + particle_lower + pt_str + "pu" + value + "_";
      else if (fixLayout)
	hitsPlotPath = "Plots/" + value + "/" + particle_lower + pt_str + "pu" + IterVariable->at(i) + "_";

      canv_nPixHits[i] -> Print((hitsPlotPath + "nPixHits.png").c_str(), "png");
    }
    
    // SCTHits vs eta
    vector<TH2F*> nSCTHits;
    vector<TCanvas*> canv_nSCTHits;

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nSCTHits.push_back(dynamic_cast<TH2F*>(inFiles[i]->Get("TrackHist_reco_all__eta_nSCTHits")));
      CheckPtr(nSCTHits[i]);
      
      if (fixPileup)
	nSCTHits[i] -> SetTitle(("Number of SCT hits vs |#eta|: " + IterVariable->at(i) + " " + particle_lower + pt_str + "pu" + value + 
				";|#eta|;").c_str());
      else if (fixLayout)
	nSCTHits[i] -> SetTitle(("Number of SCT hits vs |#eta|: " + value + " " + particle_lower + pt_str + "pu" + IterVariable->at(i) + 
				 ";|#eta|;").c_str());
      	
      canv_nSCTHits.push_back(new TCanvas());
      canv_nSCTHits[i]->cd();
      nSCTHits[i]     -> GetXaxis() -> SetRangeUser(0.0,4.0);
      nSCTHits[i]     -> DrawCopy("COL0");

      string hitsPlotPath;

      if (fixPileup) 
	hitsPlotPath = "Plots/" + IterVariable->at(i) + "/" + particle_lower + pt_str + "pu" + value + "_";
      else if (fixLayout)
	hitsPlotPath = "Plots/" + value + "/" + particle_lower + pt_str + "pu" + IterVariable->at(i) + "_";

      canv_nSCTHits[i] -> Print((hitsPlotPath + "nSCTHits.png").c_str(), "png");
    }

    // Si Hits vs eta
    vector<TH2F*> nSiHits;
    vector<TCanvas*> canv_nSiHits;

    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      nSiHits.push_back(dynamic_cast<TH2F*>(inFiles[i]->Get("TrackHist_reco_all__eta_nSiHits")));
      CheckPtr(nSiHits[i]);
      
      if (fixPileup)
	nSiHits[i] -> SetTitle(("Number of Si hits vs |#eta|: " + IterVariable->at(i) + " " + particle_lower + pt_str + "pu" + value + 
				";|#eta|;").c_str());
      else if (fixLayout)
	nSiHits[i] -> SetTitle(("Number of Si hits vs |#eta|: " + value + " " + particle_lower + pt_str + "pu" + IterVariable->at(i) + 
				 ";|#eta|;").c_str());
      	
      canv_nSiHits.push_back(new TCanvas());
      canv_nSiHits[i]->cd();
      nSiHits[i]     -> GetXaxis() -> SetRangeUser(0.0,4.0);
      nSiHits[i]     -> DrawCopy("COL0");

      string hitsPlotPath;

      if (fixPileup) 
	hitsPlotPath = "Plots/" + IterVariable->at(i) + "/" + particle_lower + pt_str + "pu" + value + "_";
      else if (fixLayout)
	hitsPlotPath = "Plots/" + value + "/" + particle_lower + pt_str + "pu" + IterVariable->at(i) + "_";
      
      canv_nSiHits[i] -> Print((hitsPlotPath + "nSiHits.png").c_str(), "png");
    }

  } //doHits
}
