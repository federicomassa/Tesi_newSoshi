#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TMath.h>
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
#include "../utils/checkPtr.h"
#include "../utils/assert.h"


using namespace std;

vector<double> makeVector(double a, double b, double c, double d) {

  vector<double> v;
  v.push_back(a);
  v.push_back(b);
  v.push_back(c);
  v.push_back(d);

  return v;
}

void fixedPt(const string& particle, const double& pt, const string& fixOption, const string& value,
	     const bool& doSigmaPlots = true, 
	     const bool& doBiasPlots = true, 
	     const bool& doEffPlots = true, 
	     const bool& doFakePlots = true, 
	     const bool& doEventPlots = true,
	     const bool& doHitsPlots = true) {

  //fixOption = "PILEUP"
  //          = "LAYOUT"
  //          = "PILEUP_LAYOUT" --> pt = dummy

  //value     = "50", "140", ... for fixOption = "PILEUP"
  //          = "ExtBrl4", "IExtBrl4", ... for fixOption = "LAYOUT"

  string particle_lower, particle_upper;
  bool fixPileup, fixLayout, fixPileup_Layout;

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
    fixPileup_Layout = false;
  }
  else if (fixOption.compare("LAYOUT") == 0) {
    fixPileup = false;
    fixLayout = true;
    fixPileup_Layout = false;
  }
  else if (fixOption.compare("PILEUP_LAYOUT") == 0) {
    fixPileup = false;
    fixLayout = false;
    fixPileup_Layout = true;
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
    imagePath = "Mixed/" + particle_lower + pt_str + "pu" + value + "_"; 
  else if (fixLayout) 
    imagePath = value + "/" + particle_lower + pt_str + "_";
  else if (fixPileup_Layout) //value = layout
    imagePath = value + "/" + particle_lower + "pu200_";

  string legendLabel, legendPostLabel;
  if (fixPileup) 
    legendLabel = "p_{T} = " + pt_str + " GeV, <#mu> = " + value + ": ";
  else if (fixLayout) 
    legendLabel = "p_{T} = " + pt_str + " GeV, " + value + ": <#mu> = ";
  else if (fixPileup_Layout) {
    legendLabel = "p_{T} = "; 
    legendPostLabel = "GeV, <#mu> = 200: " + value;
  }

  string basePath, prePath, postPath, postpostPath, postpostpostPath;

  if (fixPileup) {
    basePath = "../TrackingResults";
    prePath = basePath + "/";
    postPath = "/dataset" + particle_upper + "_" + 
      pt_str + "GeV/hist-" + particle_lower + 
      pt_str + "pu" + value + ".root"; 
  }
  else if (fixLayout) {
    basePath = "../TrackingResults";
    prePath = basePath + "/" + value;
    postPath = "/dataset" + particle_upper + "_" + 
      pt_str + "GeV/hist-" + particle_lower + 
      pt_str + "pu";  
  }
  else if (fixPileup_Layout) {
    basePath = "../TrackingResults";
    prePath = basePath + "/" + value;
    postPath = "/dataset" + particle_upper + "_";

    postpostPath = "GeV/hist-" + particle_lower;
    postpostpostPath = "pu200.root";
  }

  vector<string> Layouts;
  Layouts.push_back("ExtBrl4");
  Layouts.push_back("ExtBrl4");
  Layouts.push_back("ExtBrl4");

  vector<string> Pileup;
  Pileup.push_back("0");
  Pileup.push_back("50");
  Pileup.push_back("140");
  Pileup.push_back("200");
  //Pileup.push_back("300");

  vector<string> Pt;
  Pt.push_back("15");
  Pt.push_back("15");
  Pt.push_back("15");
  Pt.push_back("15");

  vector<string>* IterVariable = 0;

  if (fixPileup)
    IterVariable = &Layouts;
  else if (fixLayout)
    IterVariable = &Pileup;
  else if (fixPileup_Layout)
    IterVariable = &Pt;

  vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kGreen);
  colors.push_back(kBlue);
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
    else if (fixPileup_Layout)
      Path.push_back(prePath + postPath + IterVariable->at(i) + postpostPath + IterVariable->at(i) + postpostpostPath);
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

    vector<TGraphErrors*> sigQPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigQPt_abseta")));
      CheckPtr(sigQPt_abseta[i]);
    }

    vector<TGraphErrors*> sigEta_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigEta_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigEta_abseta")));
      CheckPtr(sigEta_abseta[i]);
    }

    vector<TGraphErrors*> sigPhi_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigPhi_abseta")));
      CheckPtr(sigPhi_abseta[i]);
    }

    vector<TGraphErrors*> sigD0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigD0_abseta")));
      CheckPtr(sigD0_abseta[i]);
    }

    vector<TGraphErrors*> sigZ0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sigZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sigZ0_abseta")));
      CheckPtr(sigZ0_abseta[i]);
    }

    //////////////

    vector<TGraphErrors*> sysPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysPt_abseta")));
      CheckPtr(sysPt_abseta[i]);
    }

    vector<TGraphErrors*> sysQPt_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysQPt_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysQPt_abseta")));
      CheckPtr(sysQPt_abseta[i]);
    }

    vector<TGraphErrors*> sysEta_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysEta_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysEta_abseta")));
      CheckPtr(sysEta_abseta[i]);
    }

    vector<TGraphErrors*> sysPhi_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysPhi_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysPhi_abseta")));
      CheckPtr(sysPhi_abseta[i]);
    }

    vector<TGraphErrors*> sysD0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysD0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysD0_abseta")));
      CheckPtr(sysD0_abseta[i]);
    }

    vector<TGraphErrors*> sysZ0_abseta;
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      sysZ0_abseta.push_back(dynamic_cast<TGraphErrors*>(inFiles[i]->Get("RunHist_reco_all_hard_all__sysZ0_abseta")));
      CheckPtr(sysZ0_abseta[i]);
    }

    vector<vector<TGraphErrors*> > plots;
    plots.push_back(sigPt_abseta);
    plots.push_back(sigQPt_abseta);
    plots.push_back(sigEta_abseta);
    plots.push_back(sigPhi_abseta);
    plots.push_back(sigD0_abseta);
    plots.push_back(sigZ0_abseta);
    plots.push_back(sysPt_abseta);
    plots.push_back(sysQPt_abseta);
    plots.push_back(sysEta_abseta);
    plots.push_back(sysPhi_abseta);
    plots.push_back(sysD0_abseta);
    plots.push_back(sysZ0_abseta);

    vector<string> title;
    title.push_back(";Truth |#eta|;#sigma(p_{T}) [MeV]");
    title.push_back(";Truth |#eta|;p_{T}^{truth} x #sigma(q/p_{T})");
    title.push_back(";Truth |#eta|;#sigma(#eta)/|#eta|");
    title.push_back(";Truth |#eta|;#sigma(#phi)");
    title.push_back(";Truth |#eta|;#sigma(d_{0}) [mm]");
    title.push_back(";Truth |#eta|;#sigma(z_{0}) [mm]");

    title.push_back(";Truth |#eta|;bias(p_{T}) [MeV]");
    title.push_back(";Truth |#eta|;p_{T}^{truth} x bias(q/p_{T})");
    title.push_back(";Truth |#eta|;bias(#eta)/|#eta|");
    title.push_back(";Truth |#eta|;bias(#phi)");
    title.push_back(";Truth |#eta|;bias(d_{0}) [mm]");
    title.push_back(";Truth |#eta|;bias(z_{0}) [mm]");

    vector<string> name;
    name.push_back("sigPt_abseta.png");
    name.push_back("sigQPt_abseta.png");
    name.push_back("sigEta_abseta.png");
    name.push_back("sigPhi_abseta.png");
    name.push_back("sigD0_abseta.png");
    name.push_back("sigZ0_abseta.png");

    name.push_back("sysPt_abseta.png");
    name.push_back("sysQPt_abseta.png");
    name.push_back("sysEta_abseta.png");
    name.push_back("sysPhi_abseta.png");
    name.push_back("sysD0_abseta.png");
    name.push_back("sysZ0_abseta.png");


    vector<vector<double> > legendPosition;
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.55,0.85,0.85,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));

    legendPosition.push_back(makeVector(0.15,0.45,0.45,0.15));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.55,0.45,0.85,0.15));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));
    legendPosition.push_back(makeVector(0.15,0.85,0.45,0.55));

    for (vector<vector<TGraphErrors*> >::iterator itr = plots.begin();
	 itr != plots.end();
	 itr++) {

      TCanvas* canv_abseta = new TCanvas();
      
      TMultiGraph* mg_abseta = new TMultiGraph;
      
      double minimum = 1000000000;
      double maximum = -1000000000;

      for (unsigned int i = 0; i < IterVariable->size(); i++) {
	(*itr)[i]     -> SetLineColor(colors[i]);
	(*itr)[i]     -> SetMarkerStyle(markerStyle[i]);
	(*itr)[i]     -> SetMarkerColor(colors[i]);

	if (TMath::MaxElement((*itr)[i] -> GetN(), (*itr)[i]->GetY()) > maximum)
	  maximum = TMath::MaxElement((*itr)[i] -> GetN(), (*itr)[i]->GetY());

	if (TMath::MinElement((*itr)[i] -> GetN(), (*itr)[i]->GetY()) < minimum)
	  minimum = TMath::MinElement((*itr)[i] -> GetN(), (*itr)[i]->GetY());	

	mg_abseta->Add((*itr)[i]);
      }
     
      
      mg_abseta   -> SetTitle(title[itr-plots.begin()].c_str());
      mg_abseta     -> Draw("APE");
      mg_abseta     -> GetYaxis() -> SetTitleOffset(1.3);
      //    mg_abseta     -> GetYaxis() -> SetRangeUser(0.0,4500.0);
      mg_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);
      //mg_abseta     -> GetYaxis() -> SetRangeUser(minimum - 0.1*(maximum - minimum), maximum + 0.1*(maximum - minimum));
      mg_abseta     -> GetYaxis() -> SetTitleSize(0.04);
      mg_abseta     -> GetXaxis() -> SetTitleSize(0.04);
      
      TLegend* leg_abseta = new TLegend(legendPosition[itr-plots.begin()][0], legendPosition[itr-plots.begin()][1],
					legendPosition[itr-plots.begin()][2], legendPosition[itr-plots.begin()][3]);

      for (unsigned int i = 0; i < IterVariable->size(); i++) {
	if (!fixPileup_Layout)
	  leg_abseta->AddEntry((*itr)[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
	else
	  leg_abseta->AddEntry((*itr)[i], (legendLabel + IterVariable->at(i) + legendPostLabel).c_str(), "lp");
      }
      
      leg_abseta  -> SetBorderSize(0);
      leg_abseta  -> SetFillColor(0);
      leg_abseta  -> Draw();
      
      canv_abseta -> Print((imagePath + name[itr-plots.begin()]).c_str(), "png");
      delete canv_abseta;
      
      
      
    }
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
    mg_eff_abseta   -> SetTitle(";Truth |#eta|;Efficiency");
  
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
    mg_eff_phi   -> SetTitle(";Truth #phi;Efficiency");
  
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
    mg_fakeProb_abseta   -> SetTitle(";Truth |#eta|;Fake probability");
  
    for (unsigned int i = 0; i < IterVariable->size(); i++) {
      fakeProb_abseta[i]     -> SetLineWidth(2);
      fakeProb_abseta[i]     -> SetLineColor(colors[i]);
      fakeProb_abseta[i]     -> SetMarkerStyle(markerStyle[i]);
      fakeProb_abseta[i]     -> SetMarkerColor(colors[i]);
    }
  
    for (unsigned int i = 0; i < IterVariable->size(); i++)
      mg_fakeProb_abseta->Add(fakeProb_abseta[i]);
    
  
    mg_fakeProb_abseta     -> Draw("APE");
    mg_fakeProb_abseta     -> GetYaxis() -> SetTitleOffset(1.2);
    //mg_fakeProb_abseta     -> GetYaxis() -> SetRangeUser(0.0, 0.008);
    mg_fakeProb_abseta     -> GetXaxis() -> SetRangeUser(0.0,4.0);  

    TLegend* leg_fakeProb_abseta = new TLegend(0.15,0.85,0.45,0.55);
    for (unsigned int i = 0; i < IterVariable->size(); i++) 
      leg_fakeProb_abseta->AddEntry(fakeProb_abseta[i], (legendLabel + IterVariable->at(i)).c_str(), "lp");
  
    leg_fakeProb_abseta  -> SetFillColor(0);
    leg_fakeProb_abseta  -> SetBorderSize(0);
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
