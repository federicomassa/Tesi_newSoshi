#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>


void fakeProb(const char* latexName = "latexFake.txt") {


  ofstream latex(latexName);

  std::string basePath = "../TrackingResults/";

  std::vector<std::vector<std::vector<TFile*> > > file;
  
  std::vector<std::string> layouts;
  layouts.push_back("ExtBrl4");
  layouts.push_back("ExtBrl4");
  layouts.push_back("ExtBrl4");

  std::vector<std::string> pt;
  pt.push_back("15");
  pt.push_back("15");
  pt.push_back("15");
  pt.push_back("15");

  std::vector<std::string> pileup;
  pileup.push_back("0");
  pileup.push_back("50");
  pileup.push_back("140");
  pileup.push_back("200");


  for (std::vector<std::string>::iterator layoutItr = layouts.begin();
       layoutItr != layouts.end();
       layoutItr++) {

    std::vector<std::vector<TFile*> > layoutFile;
    //std::vector<std::vector<TFile*> > layoutFile;

    for (std::vector<std::string>::iterator ptItr = pt.begin();
	 ptItr != pt.end();
	 ptItr++) {

      std::vector<TFile*> ptFile;
      //std::vector<TFile*> ptFile;

      for (std::vector<std::string>::iterator pileupItr = pileup.begin();
	   pileupItr != pileup.end();
	   pileupItr++) {

	std::string path = basePath + (*layoutItr) + "/datasetPi_" +
	  (*ptItr) + "GeV/hist-pi" + (*ptItr) + "pu" + (*pileupItr) + ".root";


	TFile* pileupFile = new TFile(path.c_str());
	//std::string pileupFile = path;

	ptFile.push_back(pileupFile);
	    
      }

      layoutFile.push_back(ptFile);

    }

    file.push_back(layoutFile);
    
  }
  
  //////////////////////////////////////
  ///DIFFERENT ORDER TO MAKE THE TABLE RIGHT


  for (std::vector<std::string>::iterator ptItr = pt.begin();
       ptItr != pt.end();
       ptItr++) {

    latex << "\\multirow{" << pileup.size() << "}{*}{" << (*ptItr) << " GeV} ";

    for (std::vector<std::string>::iterator pileupItr = pileup.begin();
	 pileupItr != pileup.end();
	 pileupItr++) {

      latex << " & " << (*pileupItr) << " & ";

      for (std::vector<std::string>::iterator layoutItr = layouts.begin();
	   layoutItr != layouts.end();
	   layoutItr++) {
	
	//eta is dummy, any reco plot is okay
	TH1F* reco = (TH1F*) file[layoutItr - layouts.begin()][ptItr - pt.begin()][pileupItr - pileup.begin()]->Get("TrackHist_reco_all__eta");
	TH1F* fake = (TH1F*) file[layoutItr - layouts.begin()][ptItr - pt.begin()][pileupItr - pileup.begin()]->Get("TrackHist_reco_fake__eta");

	double eff = fake->GetEntries()/reco->GetEntries();
	double effError = TMath::Sqrt(eff*(1.0-eff)/reco->GetEntries());

	eff = eff*100;
	effError = effError*100;

	int precision = -TMath::FloorNint(TMath::Log10(effError));

	
	latex << std::fixed << std::setprecision(precision) << eff << " $\\pm$ " << effError << " \\%";

	if (layoutItr != --layouts.end())
	  latex << " & ";
														
      }

      if (pileupItr == --pileup.end())
	latex << "\\\\ \\hline" << std::endl;
      else
	latex << "\\\\ \\cline{2-5}" << std::endl;

    }
  }



  

}
