#include <TH1F.h>
#include <THStack.h>

void provaStack2() {

  TH1F* h = new TH1F("h", "Ciao!", 100, -5, 5);
  h->FillRandom("gaus", 1000);

  TH1F* h2 = new TH1F("h2", "Ciao!", 100, -5, 5);
  h2->FillRandom("gaus", 5000);

  THStack* hs = new THStack("hs", "Prova stack");

  h->SetFillColor(kRed);
  h->SetMarkerColor(kRed);

  hs->Add(h);

  h2->SetFillColor(kBlue);
  h2->SetMarkerColor(kBlue);

  hs->Add(h2);

  hs->Draw();

}
