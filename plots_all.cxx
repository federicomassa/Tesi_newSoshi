#include "fixedPt.cxx"

void plots_all(const bool& doSigmaPlots = true,
	       const bool& doBiasPlots = true,
	       const bool& doEffPlots = true,
	       const bool& doFakePlots = true,
	       const bool& doEventPlots = true,
	       const bool& doHitsPlots = true) {
  fixedPt("pi", 15000, "PILEUP", "10", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, doHitsPlots);
  fixedPt("pi", 15000, "PILEUP", "50", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, doHitsPlots);
  fixedPt("pi", 15000, "PILEUP", "140", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, doHitsPlots);
  fixedPt("pi", 15000, "PILEUP", "200", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, doHitsPlots);
  fixedPt("pi", 15000, "PILEUP", "300", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, doHitsPlots);

  fixedPt("pi", 15000, "LAYOUT", "ExtBrl4", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, false);
  fixedPt("pi", 15000, "LAYOUT", "IExtBrl4", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, false);
  fixedPt("pi", 15000, "LAYOUT", "InclBrl4", doSigmaPlots, doBiasPlots, doEffPlots, doFakePlots, doEventPlots, false);
}
