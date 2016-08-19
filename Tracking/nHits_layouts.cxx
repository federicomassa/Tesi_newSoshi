#include "nHits.cxx"

#include <TMath.h>
#include <fstream>
#include <iomanip>

void nHits_layouts() {

  ofstream f("latexHits.txt");
  double pixMean, SCTMean, siMean;
  double pixMeanError, SCTMeanError, siMeanError;
  double pixPrecision, SCTPrecision, siPrecision;
  
  nHits("ExtBrl4", &pixMean, &pixMeanError, &SCTMean, &SCTMeanError, &siMean, &siMeanError);
  pixPrecision = -TMath::FloorNint(TMath::Log10(pixMeanError));
  SCTPrecision = -TMath::FloorNint(TMath::Log10(SCTMeanError));
  siPrecision = -TMath::FloorNint(TMath::Log10(siMeanError));
  f << "\\textbf{ExtBrl4} & " << std::fixed << std::setprecision(pixPrecision) << pixMean << " $\\pm$ " << pixMeanError << " & " << std::setprecision(SCTPrecision) << SCTMean << " $\\pm$ " << SCTMeanError << " & " << std::setprecision(siPrecision) << siMean << " $\\pm$ " << siMeanError << " \\\\ \\hline" << std::endl;

  nHits("ExtBrl4", &pixMean, &pixMeanError, &SCTMean, &SCTMeanError, &siMean, &siMeanError);
  pixPrecision = -TMath::FloorNint(TMath::Log10(pixMeanError));
  SCTPrecision = -TMath::FloorNint(TMath::Log10(SCTMeanError));
  siPrecision = -TMath::FloorNint(TMath::Log10(siMeanError));

  f << "\\textbf{ExtBrl4} & " << std::fixed << std::setprecision(pixPrecision) << pixMean << " $\\pm$ " << pixMeanError << " & " << std::setprecision(SCTPrecision) << SCTMean << " $\\pm$ " << SCTMeanError << " & " << std::setprecision(siPrecision) << siMean << " $\\pm$ " << siMeanError << " \\\\ \\hline" << std::endl;

  nHits("ExtBrl4", &pixMean, &pixMeanError, &SCTMean, &SCTMeanError, &siMean, &siMeanError);
  pixPrecision = -TMath::FloorNint(TMath::Log10(pixMeanError));
  SCTPrecision = -TMath::FloorNint(TMath::Log10(SCTMeanError));
  siPrecision = -TMath::FloorNint(TMath::Log10(siMeanError));

  f << "\\textbf{ExtBrl4} & " << std::fixed << std::setprecision(pixPrecision) << pixMean << " $\\pm$ " << pixMeanError << " & " << std::setprecision(SCTPrecision) << SCTMean << " $\\pm$ " << SCTMeanError << " & " << std::setprecision(siPrecision) << siMean << " $\\pm$ " << siMeanError << " \\\\ \\hline" << std::endl;

  f.close();

}
