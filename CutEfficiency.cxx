#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <iostream>

void CutEfficiency(const char* fname) {

  TFile* f = new TFile(fname);
  TH1F* recoCutCode = (TH1F*) f->Get("EventHist_all__recoCutCode");
  TH1F* recoCutCode27 = (TH1F*) f->Get("EventHist_all__recoCutCode27");
  TH1F* recoCutCode32 = (TH1F*) f->Get("EventHist_all__recoCutCode32");
  TH1F* recoCutCode4 = (TH1F*) f->Get("EventHist_all__recoCutCode4");

  std::vector<std::string> cutName;
  cutName.push_back("Hits              ");
  cutName.push_back("PtMin             ");
  cutName.push_back("dZ                ");
  cutName.push_back("isolation         ");
  cutName.push_back("Pt1               ");
  cutName.push_back("Pt2               ");
  cutName.push_back("Pt3               ");
  cutName.push_back("Pt4               ");
  cutName.push_back("DR                ");
  cutName.push_back("Eta offshell      ");
  cutName.push_back("Onshell mass      ");
  cutName.push_back("Offshell mass     ");
  {
    ///////// NO ETA SELECTION ////////

    double entriesOld = recoCutCode->GetEntries() - recoCutCode->GetBinContent(1);
    std::cout <<    "Initial:          " << entriesOld << std::endl;
    
    double entries = entriesOld - recoCutCode->GetBinContent(3);
    
    for (auto itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries << '\t' << entries/entriesOld << std::endl;
      entriesOld = entries;
      entries -= recoCutCode->GetBinContent(4 + itr - cutName.begin());
    }

    std::cout << std::endl;
    //////////////////////////////////
    
  }


  {
    ///////// ETA 2.7 ////////

    double entriesOld = recoCutCode27->GetEntries() - recoCutCode27->GetBinContent(1);
    std::cout <<    "Initial:          " << entriesOld << std::endl;
    
    double entries = entriesOld - recoCutCode27->GetBinContent(3);
    
    for (auto itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries << '\t' << entries/entriesOld << std::endl;
      entriesOld = entries;
      entries -= recoCutCode27->GetBinContent(4 + itr - cutName.begin());
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }

  {
    ///////// ETA 3.2 ////////

    double entriesOld = recoCutCode32->GetEntries() - recoCutCode32->GetBinContent(1);
    std::cout <<    "Initial:          " << entriesOld << std::endl;
    
    double entries = entriesOld - recoCutCode32->GetBinContent(3);
    
    for (auto itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries << '\t' << entries/entriesOld << std::endl;
      entriesOld = entries;
      entries -= recoCutCode32->GetBinContent(4 + itr - cutName.begin());
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }


  {
    ///////// ETA 4.0 ////////

    double entriesOld = recoCutCode4->GetEntries() - recoCutCode4->GetBinContent(1);
    std::cout <<    "Initial:          " << entriesOld << std::endl;
    
    double entries = entriesOld - recoCutCode4->GetBinContent(3);
    
    for (auto itr = cutName.begin();
	 itr != cutName.end();
	 itr++) {
      std::cout << (*itr) << entries << '\t' << entries/entriesOld << std::endl;
      entriesOld = entries;
      entries -= recoCutCode4->GetBinContent(4 + itr - cutName.begin());
    }
    
    std::cout << std::endl;
    //////////////////////////////////
    
  }

  

}
