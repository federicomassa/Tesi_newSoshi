#include <HistManager/TrackHists.h>
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/Vertex.h"
#include <HistManager/TrackHelper.h>
#include <HistManager/assert.h>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <cstdlib>

const bool debug = false;

void PrintMessage(const char* message) {
  if (debug)
    std::cout << message << std::endl;
}

TrackHists::TrackHists(TString name) {
  m_baseName = name;
  m_name = "TrackHist_" + name;
  m_label = "Track "; // don't forget the space
  m_isPseudo = false;
  if (m_name.Contains("pseudo")) { m_label = "Pseudotrack "; m_isPseudo = true; }
}

TrackHists::~TrackHists() {}

void TrackHists::BookHists() {

  std::cout << "Creating histograms for " << m_name << std::endl;

  m_muValue = declare1D(m_name, "muValue", m_label + "pilueup  #mu",30,0.0,30.0);
  m_lumibk  = declare1D(m_name, "lumibk",  m_label + "lumi block",100,0,1000);
  m_pixClus = declare1D(m_name, "pixClus",  m_label + "N pixel clusters",200,0,4000);

  // from truth particle itself
  m_prodR         = declare1D(m_name, "prodR",         "production R [mm]", 110,    0.0, 1100.0);
  m_prodZ         = declare1D(m_name, "prodZ",         "production Z [mm]", 310,    0.0, 3100.0);
  m_truthPt       = declare1D(m_name, "truthPt",       "Truth p_{T} [GeV]", 200,    0.0, 200.0); 
  m_truthPtnarrow = declare1D(m_name, "truthPtnarrow", "Truth p_{T} [GeV]", 200,    0.0,   50.0); 
  m_truthEta      = declare1D(m_name, "truthEta",      "Truth #eta",         80,   -4.0,    4.0); 
  m_truthAbseta   = declare1D(m_name, "truthAbseta",       "|#eta|", 20, 0.0, 4.0); 
  m_truthPhi      = declare1D(m_name, "truthPhi",      "Truth #phi",         32,   -3.2,    3.2); 
  m_truthD0       = declare1D(m_name, "truthD0",       "Truth d_{0} [mm]",	100,   -0.1,    0.1); 
  m_truthD0_wide  = declare1D(m_name, "truthD0_wide",  "Truth d_{0} [mm]",  100,   -2.0,    2.0); 
  m_truthZ0       = declare1D(m_name, "truthZ0",       "Truth z_{0} [mm]",  100, -200.0,  200.0); 

  m_biasPt  = declare1D(m_name, "biasPt",  "#sigma(p_{T}) [MeV]",        100, -8000.0,8000.0); 
  m_biasQPt = declare1D(m_name, "biasQPt", "p_{T}^{truth} x #sigma(q/p_{T})", 100,    -0.2,   0.2); 
  m_biasEta = declare1D(m_name, "biasEta", "#sigma(#eta)", 			         100,   -2e-3,  2e-3); 
  m_biasPhi = declare1D(m_name, "biasPhi", "#sigma(#phi)", 			         100,	  -1.0e-3,  1.0e-3); 
  m_biasD0  = declare1D(m_name, "biasD0",  "#sigma(d_{0}) [mm]",         100,    -0.2,   0.2); 
  m_biasZ0  = declare1D(m_name, "biasZ0",  "#sigma(z_{0}) [mm]",         100,    -0.8,   0.8); 

  for (int i = 0; i < m_etaVectorSize; i++) {

    std::stringstream ss;

    std::string count = "";
    std::string etaLimitDown = "";
    std::string etaLimitUp = "";

    ss << i+1;
    ss >> count;
    ss.clear();

    ss << 0.2*double(i+1);
    ss >> etaLimitUp;
    ss.clear();

    ss << 0.2*double(i);
    ss >> etaLimitDown;
    ss.clear();

    const std::string PtTitle = "#Delta(p_{T})( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";
    const std::string QPtTitle = "#Delta(q/p_{T}) x p_{T}( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";
    const std::string EtaTitle = "#Delta(#eta)/|#eta_{truth}|( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";
    const std::string PhiTitle = "#Delta(#phi)( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";
    const std::string D0Title = "#Delta(D0)( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";
    const std::string Z0Title = "#Delta(Z0)( " + etaLimitDown + " < |#eta| < " + etaLimitUp + ")";

    const std::string PtUnits = "[MeV]";
    const std::string D0Units = "[mm]";
    const std::string Z0Units = "[mm]";

    const std::string PtName = "biasPt_bin";
    const std::string QPtName = "biasQPt_bin";
    const std::string EtaName = "biasEta_bin";
    const std::string PhiName = "biasPhi_bin";
    const std::string D0Name = "biasD0_bin";
    const std::string Z0Name = "biasZ0_bin";

    m_biasPt_abseta.push_back(declare1D(m_name, (PtName + count).c_str(),  (PtTitle + PtUnits).c_str(), 200,-16000.0,16000.0));
    m_biasQPt_abseta.push_back(declare1D(m_name, (QPtName + count).c_str(),  (QPtTitle).c_str(), 200,-0.4,0.4));
    m_biasEta_abseta.push_back(declare1D(m_name, (EtaName + count).c_str(),  (EtaTitle).c_str(), 200,-0.08,0.08));
    m_biasPhi_abseta.push_back(declare1D(m_name, (PhiName + count).c_str(),  (PhiTitle).c_str(), 200,-5e-3,5e-3));
    m_biasD0_abseta.push_back(declare1D(m_name, (D0Name + count).c_str(),  (D0Title + D0Units).c_str(), 200,-1.0,1.0));
    m_biasZ0_abseta.push_back(declare1D(m_name, (Z0Name + count).c_str(),  (Z0Title + Z0Units).c_str(), 100,-0.8,0.8));

    m_biasPt_poseta.push_back(declare1D(m_name, (PtName + count + "_pos").c_str(),  (PtTitle + PtUnits).c_str(), 200,-16000.0,16000.0));
    m_biasQPt_poseta.push_back(declare1D(m_name, (QPtName + count + "_pos").c_str(),  (QPtTitle).c_str(), 200,-0.4,0.4));
    m_biasEta_poseta.push_back(declare1D(m_name, (EtaName + count + "_pos").c_str(),  (EtaTitle).c_str(), 200,-0.4,0.4));
    m_biasPhi_poseta.push_back(declare1D(m_name, (PhiName + count + "_pos").c_str(),  (PhiTitle).c_str(), 200,-5e-3,5e-3));
    m_biasD0_poseta.push_back(declare1D(m_name, (D0Name + count + "_pos").c_str(),  (D0Title + D0Units).c_str(), 200,-1.0,1.0));
    m_biasZ0_poseta.push_back(declare1D(m_name, (Z0Name + count + "_pos").c_str(),  (Z0Title + Z0Units).c_str(), 200,-200.0,200.0));
 
    m_biasPt_negeta.push_back(declare1D(m_name, (PtName + count + "_neg").c_str(),  (PtTitle + PtUnits).c_str(), 200,-16000.0,16000.0));
    m_biasQPt_negeta.push_back(declare1D(m_name, (QPtName + count + "_neg").c_str(),  (QPtTitle).c_str(), 200,-0.4,0.4));
    m_biasEta_negeta.push_back(declare1D(m_name, (EtaName + count + "_neg").c_str(),  (EtaTitle).c_str(), 200,-0.4,0.4));
    m_biasPhi_negeta.push_back(declare1D(m_name, (PhiName + count + "_neg").c_str(),  (PhiTitle).c_str(), 200,-5e-3,5e-3));
    m_biasD0_negeta.push_back(declare1D(m_name, (D0Name + count + "_neg").c_str(),  (D0Title + D0Units).c_str(), 200,-1.0,1.0));
    m_biasZ0_negeta.push_back(declare1D(m_name, (Z0Name + count + "_neg").c_str(),  (Z0Title + Z0Units).c_str(), 200,-200.0,200.0));
 
  }

  /* track parameterization */
  m_qoverp      = declare1D(m_name, "qoverp",      "#frac{q}{p} [GeV]", 200,  -2e-3,   2e-3); 
  m_pt          = declare1D(m_name, "pt",          "p_{T} [GeV]",       200,    0.0, 2000.0); 
  m_ptnarrow    = declare1D(m_name, "ptnarrow",    "p_{T} [GeV]",       200,    0.0,   50.0); 
  m_eta         = declare1D(m_name, "eta",         "#eta",               40,   -4.0,    4.0); 
  m_abseta      = declare1D(m_name, "abseta",      "#eta",               20,    0.0,    4.0); 
  m_phi         = declare1D(m_name, "phi",         "#phi",               32,   -3.2,    3.2); 
  m_d0          = declare1D(m_name, "d0",          "d_{0} [mm]",        100,   -0.2,    0.2); 
  m_d0_manual   = declare1D(m_name, "d0_manual",   "d_{0} [mm]",        100,   -0.2,    0.2); 
  m_d0_truthd0  = declare2D(m_name, "d0_truthd0",  "track d_{0} [mm]", "bias d_{0} [mm]", 100, -0.2, 0.2, 100, -0.2, 0.2);
  m_z0          = declare1D(m_name, "z0",          "z_{0} [mm]",        100, -200.0,  200.0); 
  m_z0_abseta   = declare2D(m_name, "z0_abseta",   "|#eta|", "z_{0} [mm]", 20,0.0,4.0, 100, -200, 200);
  m_z0_truthz0  = declare2D(m_name, "z0_truthz0",  "track z_{0} [mm]", "bias z_{0} [mm]", 100, -200.0, 200.0, 100, -200.0, 200.0); 
  m_z0Corr1     = declare1D(m_name, "z0Corr1",     "z_{0} at PV [mm]",   80,   -8.0,    8.0); 
  m_z0Corr2     = declare1D(m_name, "z0Corr2",     "z_{0} at PV [mm]",   80,   -8.0,    8.0); 

  /* track parameterization errors */
  m_qoverp_err = declare1D(m_name, "qoverp_err", "#sigma #frac{q}{p} [GeV]", 200, 0.0,2e-6); 
  m_eta_err    = declare1D(m_name, "eta_err",    "#sigma #eta",              100, 0.0, 0.1); 
  m_d0_err     = declare1D(m_name, "d0_err",     "#sigma d_{0}",             100, 0.0, 0.04); 
  m_d0signif   = declare1D(m_name, "d0signif",   "d_{0} Significance",       100,-5.0, 5.0); 
  m_z0_err     = declare1D(m_name, "z0_err",     "#sigma z_{0}",             100, 0.0, 0.4); 
  m_phi_err    = declare1D(m_name, "phi_err",    "#sigma #phi",              100,-3.5, 3.5); 

  m_d0sign      = declare1D(m_name, "d0sign",      "signed IP [mm]",         100,   -10.0,   15.0); 
  m_d0signsig   = declare1D(m_name, "d0signsig",   "signed IP significance", 150,   -10.0,   20.0); 
  m_d0signsigP  = declare1D(m_name, "d0signsigP",  "signed IP significance", 150,   -10.0,   20.0); 
  m_d0signsigN  = declare1D(m_name, "d0signsigN",  "signed IP significance", 150,   -10.0,   20.0); 

  /* track properties */
  m_nPixHits = declare1D(m_name, "nPixHits", "N Pixel Hits",   16, -0.5, 15.5); 
  m_nSCTHits = declare1D(m_name, "nSCTHits", "N SCT Hits",	   25, -0.5, 24.5); 
  m_nSiHits  = declare1D(m_name, "nSiHits",  "N Silicon Hits", 30, -0.5, 29.5);

  m_eta_nPixHits = declare2D(m_name, "eta_nPixHits", "|#eta|", "N Pixel Hits",   20, 0, 4.0, 25, -0.5, 24.5); 
  m_eta_nSCTHits = declare2D(m_name, "eta_nSCTHits", "|#eta|", "N SCT Hits",	   20, 0, 4.0, 20, -0.5, 19.5); 
  m_eta_nSiHits  = declare2D(m_name, "eta_nSiHits",  "|#eta|", "N Silicon Hits", 20, 0, 4.0, 30, -0.5, 29.5);

  m_nGangedPix         = declare1D(m_name, "nGangedPix",         m_label + "N Pixel Ganged Hits",     8, -0.5,  7.5); 
  m_nGangedPixFF       = declare1D(m_name, "nGangedPixFF",       m_label + "N Pixel Ganged FF Hits",  6, -0.5,  5.5); 
  m_nPixLay            = declare1D(m_name, "nPixLay",            m_label + "N Pixel Layers",         10, -0.5,  9.5); 
  m_nPixSharedHits     = declare1D(m_name, "nPixSharedHits",     m_label + "N Pixel Shared Hits",     8, -0.5,  7.5); 
  m_nPixSplitHits      = declare1D(m_name, "nPixSplitHits",      m_label + "N Pixel Split Hits",      8, -0.5,  7.5); 
  m_nPixSpoiltHits     = declare1D(m_name, "nPixSpoiltHits",     m_label + "N Pixel Spoilt Hits",     8, -0.5,  7.5); 
  m_nPixOutliers       = declare1D(m_name, "nPixOutliers",       m_label + "N Pixel Outliers",        8, -0.5,  7.5); 
  m_nPixHoles          = declare1D(m_name, "nPixHoles",          m_label + "N Pixel Holes",           8, -0.5,  7.5); 
  m_nPixelDeadSensors  = declare1D(m_name, "nPixelDeadSensors",  m_label + "N Pixel Dead Sensors",    8, -0.5,  7.5); 
  m_nSCTSharedHits     = declare1D(m_name, "nSCTSharedHits",     m_label + "N SCT Shared Hits",      11, -0.5, 10.5); 
  m_nSCTSpoiltHits     = declare1D(m_name, "nSCTSpoiltHits",     m_label + "N SCT Spoilt Hits",      11, -0.5, 10.5); 
  m_nSCTOutliers       = declare1D(m_name, "nSCTOutliers",	     m_label + "N SCT Outliers",          8, -0.5,  7.5); 
  m_nSCTHoles          = declare1D(m_name, "nSCTHoles",	         m_label + "N SCT Holes",             8, -0.5,  7.5); 
  m_nSCTDoubleHoles    = declare1D(m_name, "nSCTDoubleHoles",    m_label + "N SCT Double Holes",      8, -0.5,  7.5); 
  m_nSCTDeadSensors    = declare1D(m_name, "nSCTDeadSensors",    m_label + "N SCT Dead Sensors",      5, -0.5,  4.5);
  m_nSiDeadSensors     = declare1D(m_name, "nSiDeadSensors",     m_label + "N Silicon Dead Sensors",  9, -0.5,  8.5);

  m_nTRTHits     = declare1D(m_name, "nTRTHits",     m_label + "N TRT Hits",                41, -0.5, 40.5 ); 
  m_nTRTOutliers = declare1D(m_name, "nTRTOutliers", m_label + "N TRT Outliers",            11, -0.5, 10.5 ); 
  m_nTRTHoles    = declare1D(m_name, "nTRTHoles",    m_label + "N TRT Holes",               11, -0.5, 10.5 ); 
  m_nTRTHTHits   = declare1D(m_name, "nTRTHTHits",   m_label + "N High Threshold TRT Hits", 31, -0.5, 30.5 ); 

  m_chiSqPerDof  = declare1D(m_name, "chiSqPerDof",  m_label + "#chi^{2}/DOF", 120,  0.0, 12.0);
  m_Dof          = declare1D(m_name, "Dof",          m_label + "DOF", 100, 0,100);
  m_nOutliers    = declare1D(m_name, "nOutliers",    m_label + "N Outliers",    20, -0.5, 19.5);
  m_stdDevChi2OS = declare1D(m_name, "stdDevChi2OS", m_label + "100x Std. Dev. #chi^{2} from surface", 100, 0.0, 500.0);

  m_eProbComb = declare1D(m_name, "eProbComb", m_label + "e Prob Comb", 20, 0.0, 1.0); 
  m_eProbHT   = declare1D(m_name, "eProbHT",   m_label + "e Prob HT",   20, 0.0, 1.0); 
  m_eProbToT  = declare1D(m_name, "eProbToT",  m_label + "e Prob ToT",  20, 0.0, 1.0); 
  m_eProbBrem = declare1D(m_name, "eProbBrem", m_label + "e Prob Brem", 20, 0.0, 1.0); 

  m_truthMatchProb = declare1D(m_name, "truthMatchProb", m_label + "Truth Match \"Prob\"", 20, 0.0, 1.0);

  m_matchingDR = declare1D(m_name, "matchingDR", "#Delta R", 100, 0.0, 0.1); 

  m_truthAbseta_dr001 = declare1D(m_name, "truthAbseta_dr001", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr002 = declare1D(m_name, "truthAbseta_dr002", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr003 = declare1D(m_name, "truthAbseta_dr003", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr004 = declare1D(m_name, "truthAbseta_dr004", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr005 = declare1D(m_name, "truthAbseta_dr005", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr006 = declare1D(m_name, "truthAbseta_dr006", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr007 = declare1D(m_name, "truthAbseta_dr007", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr008 = declare1D(m_name, "truthAbseta_dr008", "|#eta|", 25, 0.0, 5.0); 
  m_truthAbseta_dr009 = declare1D(m_name, "truthAbseta_dr009", "|#eta|", 25, 0.0, 5.0); 


  // Cluster study
  m_IBLCharge       = declare1D(m_name, "IBLCharge",       "charge at IBL",       100,    0.,     250000.0 ); 
  m_IBLHitSize      = declare1D(m_name, "IBLHitSize",      "IBL Hit: Size",       20,   0.,     20. ); 
  m_IBLHitSizeZ     = declare1D(m_name, "IBLHitSizeZ",     "IBL Hit: Size(z)",    8,    0.,     8. ); 
  m_IBLHitSizePhi   = declare1D(m_name, "IBLHitSizePhi",   "IBL Hit: Size(#phi)", 10,   0.,     10. ); 
  m_IBLlocalPhi     = declare1D(m_name, "IBLlocalPhi",     "IBL local #phi",   50,  0.0, 1.0);
  m_IBLlocalTheta   = declare1D(m_name, "IBLlocalTheta",   "IBL local #theta", 32, -1.6, 1.6);
  m_IBLlocalPhi_sizePhi = declare2D(m_name, "IBLlocalPhi_sizePhi", "IBL local #phi",   "IBL Size #phi", 50,  0.0, 1.0, 10, -0.5, 9.5);
  m_IBLlocalTheta_sizeZ = declare2D(m_name, "IBLlocalTheta_sizeZ", "IBL local #theta", "IBL Size z",    32, -1.6, 1.6, 10, -0.5, 9.5);

  m_IBLbiasedResidualX   = declare1D(m_name, "IBLbiasedResidualX", "IBL biased residual x[mm]", 100, -0.2, 0.2);
  m_IBLbiasedResidualY   = declare1D(m_name, "IBLbiasedResidualY", "IBL biased residual y[mm]", 100, -0.2, 0.2);
  m_IBLunbiasedResidualX = declare1D(m_name, "IBLunbiasedResidualX", "IBL unbiased residual x[mm]", 100, -0.2, 0.2);
  m_IBLunbiasedResidualY = declare1D(m_name, "IBLunbiasedResidualY", "IBL unbiased residual y[mm]", 100, -0.2, 0.2);
  m_IBLResidualX         = declare1D(m_name, "IBLResidualX", "IBL residual x[mm]", 100, -0.2, 0.2);
  m_IBLResidualY         = declare1D(m_name, "IBLResidualY", "IBL residual y[mm]", 100, -0.2, 0.2);

  m_BLayerCharge       = declare1D(m_name, "BLayerCharge",     "charge at B-Layer",       100,    0.,     250000.0 ); 
  m_BLayerHitSize      = declare1D(m_name, "BLayerHitSize",    "B-Layer Hit: Size",       20,   0.,     20. ); 
  m_BLayerHitSizeZ     = declare1D(m_name, "BLayerHitSizeZ",   "B-Layer Hit: Size(z)",    8,    0.,     8. ); 
  m_BLayerHitSizePhi   = declare1D(m_name, "BLayerHitSizePhi", "B-Layer Hit: Size(#phi)", 10,   0.,     10. ); 
  m_BLayerlocalPhi     = declare1D(m_name, "BLayerlocalPhi",   "B-Layer local #phi",   50,  0.0, 1.0);
  m_BLayerlocalTheta   = declare1D(m_name, "BLayerlocalTheta", "B-Layer local #theta", 32, -1.6, 1.6);
  m_BLayerlocalPhi_sizePhi = declare2D(m_name, "BLayerlocalPhi_sizePhi", "B-Layer local #phi",   "B-Layer Size #phi", 50,  0.0, 1.0, 10, -0.5, 9.5);
  m_BLayerlocalTheta_sizeZ = declare2D(m_name, "BLayerlocalTheta_sizeZ", "B-Layer local #theta", "B-Layer Size z",    32, -1.6, 1.6, 10, -0.5, 9.5);

  m_BLayerbiasedResidualX   = declare1D(m_name, "BLayerbiasedResidualX", "BLayer biased residual x[mm]", 100, -0.2, 0.2);
  m_BLayerbiasedResidualY   = declare1D(m_name, "BLayerbiasedResidualY", "BLayer biased residual y[mm]", 100, -0.2, 0.2);
  m_BLayerunbiasedResidualX = declare1D(m_name, "BLayerunbiasedResidualX", "BLayer unbiased residual x[mm]", 100, -0.2, 0.2);
  m_BLayerunbiasedResidualY = declare1D(m_name, "BLayerunbiasedResidualY", "BLayer unbiased residual y[mm]", 100, -0.2, 0.2);
  m_BLayerResidualX         = declare1D(m_name, "BLayerResidualX", "BLayer residual x[mm]", 100, -0.2, 0.2);
  m_BLayerResidualY         = declare1D(m_name, "BLayerResidualY", "BLayer residual y[mm]", 100, -0.2, 0.2);

} // BookHists

void TrackHists::FillHists(const xAOD::TrackParticle* trk, float weight) const {

  if(m_usedBarcodes.size() > 100) { 
    std::cout << "WARNING::TrackHists " << m_name << " has " << m_usedBarcodes.size() 
      << " barcodes marked as used" << std::endl;
  }

  // only fill histogram 1 time per truth particle
//  if(  trk->isAvailable< int >("barcode")  ) {
//    // if used already - return!!
//    if( m_usedBarcodes.find((int)trk->auxdata< int >( "barcode" )) != m_usedBarcodes.end() ){
//      return;
//    }
//    m_usedBarcodes.insert( (int)trk->auxdata< int >( "barcode" ) );
//  }

  if (trk->isAvailable<float>("muTrk")) {
    float muValue = trk->auxdata<float>("muTrk");
    int   lumibk  = trk->auxdata<int>("lbTrk");
    int   pixClus = trk->auxdata<int>("pixClusTrk");
    m_muValue -> Fill(muValue,weight);
    m_lumibk  -> Fill(1.0*lumibk,weight);
    m_pixClus -> Fill(1.0*pixClus,weight);
  }


  uint8_t getInt(0);   // for accessing summary information
  float   getFlt(0.0); // for accessing summary information
  const xAOD::ParametersCovMatrix_t covTrk = trk->definingParametersCovMatrix();
  float d0err = sqrt(covTrk(0,0));
  float d0signif = trk->d0()/d0err;
  float z0err = sqrt(covTrk(1,1));

  /* track parameterization */
  m_qoverp       ->Fill(trk->qOverP(),weight); 
  m_pt           ->Fill(trk->pt()*1e-3,weight); 
  m_ptnarrow     ->Fill(trk->pt()*1e-3,weight); 
  m_eta          ->Fill(trk->eta(),weight); 
  m_abseta       ->Fill(TMath::Abs(trk->eta()),weight); 
  m_d0           ->Fill(trk->d0(),weight); 

  if (trk->vertex()) {
    m_d0_manual    ->Fill(TMath::Sign(
				      TMath::Sqrt(
						  TMath::Power(trk->vertex()->x(),2) +
						  TMath::Power(trk->vertex()->y(),2)),
				      trk->vertex()->x()*trk->pt()*TMath::Sin(trk->phi()) -
				      trk->vertex()->y()*trk->pt()*TMath::Cos(trk->phi())),
			  weight);
  } 

  m_z0           ->Fill(trk->z0(),weight); 
  m_z0_abseta    ->Fill(TMath::Abs(trk->eta()), trk->z0());
  m_phi          ->Fill(trk->phi(),weight); 

  float signQ = 0.0;
  if (trk->isAvailable<int>("signToQuark")) { signQ = trk->auxdata<int>("signToQuark"); }

  m_d0sign    ->Fill(trk->d0()*signQ,weight); 
  m_d0signsig ->Fill(trk->d0()*signQ/d0err,weight); 

  if (trk->d0()*signQ/d0err>0.0) {
    m_d0signsigP ->Fill(trk->d0()*signQ/d0err,weight); 
  }
  else {
    m_d0signsigN ->Fill(trk->d0()*signQ/d0err,weight); 
  }

  if (trk->isAvailable<float>("z0Corr1")) { m_z0Corr1 -> Fill(trk->auxdata<float>("z0Corr1"),weight); }
  if (trk->isAvailable<float>("z0Corr2")) { m_z0Corr2 -> Fill(trk->auxdata<float>("z0Corr2"),weight); }

  /* track parameterization errors */
  m_qoverp_err -> Fill(sqrt(covTrk(4,4)),weight); 
  m_eta_err    -> Fill(sqrt(covTrk(3,3)),weight);  // really theta - FIXME
  m_d0_err     -> Fill(d0err,weight);
  m_d0signif   -> Fill(d0signif,weight);
  m_z0_err     -> Fill(sqrt(covTrk(1,1)),weight); 
  m_phi_err    -> Fill(sqrt(covTrk(2,2)),weight); 

  //vx, vy, vz
  //parameterX, parameterY, parameterZ
  //parameterPX, parameterPY, parameterPZ
  //trackParameterCovarianceMatrices
  //parameterPosition
  //chiSquared  -> function
  //numberDoF   -> function
  //chiSqPerDoF
  //trackFitter
  //particleHypothesis
  //trackProperties
  //patternRecoInfo

  /* track properties */
  PrintMessage("Track Properties...");
  trk->summaryValue(getInt,xAOD::numberOfPixelHits); int nPixHits = getInt;
  trk->summaryValue(getInt,xAOD::numberOfSCTHits);   int nSCTHits = getInt;
  int nSiHits  = xAOD::TrackHelper::numberOfSiHits(trk);
  m_nPixHits -> Fill(1.0*nPixHits);
  m_nSCTHits -> Fill(1.0*nSCTHits); 
  m_nSiHits  -> Fill(1.0*nSiHits);
  
  PrintMessage("Eta Properties...");
  m_eta_nPixHits -> Fill(TMath::Abs(trk->eta()),1.0*nPixHits);
  m_eta_nSCTHits -> Fill(TMath::Abs(trk->eta()),1.0*nSCTHits); 
  m_eta_nSiHits  -> Fill(TMath::Abs(trk->eta()),1.0*nSiHits);

  trk->summaryValue(getInt,xAOD::numberOfGangedPixels);       m_nGangedPix        -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfGangedFlaggedFakes); m_nGangedPixFF      -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfContribPixelLayers); m_nPixLay           -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelSharedHits);    m_nPixSharedHits    -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelSplitHits);     m_nPixSplitHits     -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelSpoiltHits);    m_nPixSpoiltHits    -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelOutliers);      m_nPixOutliers      -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelHoles);         m_nPixHoles         -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfPixelDeadSensors);   m_nPixelDeadSensors -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfSCTSharedHits);      m_nSCTSharedHits    -> Fill(getInt,weight); 
  trk->summaryValue(getInt,xAOD::numberOfSCTSpoiltHits);      m_nSCTSpoiltHits    -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfSCTOutliers);        m_nSCTOutliers      -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfSCTHoles);           m_nSCTHoles         -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfSCTDoubleHoles);     m_nSCTDoubleHoles   -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfSCTDeadSensors);     m_nSCTDeadSensors   -> Fill(getInt,weight);

  PrintMessage("Si Dead Sensors Properties...");
  m_nSiDeadSensors -> Fill(xAOD::TrackHelper::numberOfSiDeadSensors(trk),weight);

  trk->summaryValue(getInt,xAOD::numberOfTRTHits);              m_nTRTHits     -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::numberOfTRTOutliers);          m_nTRTOutliers -> Fill(getInt,weight); 
  trk->summaryValue(getInt,xAOD::numberOfTRTHoles);             m_nTRTHoles    -> Fill(getInt,weight); 
  trk->summaryValue(getInt,xAOD::numberOfTRTHighThresholdHits); m_nTRTHTHits   -> Fill(getInt,weight); 
  //numberOfTRTHighThresholdOutliers
  //numberOfTRTDeadStraws
  //numberOfTRTTubeHits
  //numberOfTRTXenonHits
  //numberOfPrecisionLayers
  //numberOfPrecisionHoleLayers
  //numberOfPhiLayers
  //numberOfPhiHoleLayers
  //numberOfTriggerEtaLayers
  //numberOfTriggerEtaHoleLayers

  PrintMessage("Chisq/dof Properties...");
  m_chiSqPerDof -> Fill(xAOD::TrackHelper::chiSqPerDoF(trk),weight);
  m_Dof         -> Fill(trk->numberDoF(), weight);

  trk->summaryValue(getInt,xAOD::numberOfOutliersOnTrack);   m_nOutliers    -> Fill(getInt,weight);
  trk->summaryValue(getInt,xAOD::standardDeviationOfChi2OS); m_stdDevChi2OS -> Fill(getInt,weight);

  trk->summaryValue(getFlt,xAOD::eProbabilityComb); m_eProbComb -> Fill(getFlt,weight);
  trk->summaryValue(getFlt,xAOD::eProbabilityHT);   m_eProbHT   -> Fill(getFlt,weight);
  trk->summaryValue(getFlt,xAOD::eProbabilityToT);  m_eProbToT  -> Fill(getFlt,weight);
  trk->summaryValue(getFlt,xAOD::eProbabilityBrem); m_eProbBrem -> Fill(getFlt,weight);

  //===========================
  // Efficiency and fake study
  //===========================
  PrintMessage("Efficiency Properties...");
  const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(trk);

  if (truthParticle) {

    if (truthParticle->hasProdVtx()) {
      m_prodR -> Fill(truthParticle->prodVtx()->perp(),weight);
      m_prodZ -> Fill(truthParticle->prodVtx()->z(),weight);
    }

    float truthd0 = 0.0;
    if (truthParticle->isAvailable<float>("d0")) {
      truthd0 = truthParticle->auxdata<float>("d0");
    }
    // else if(truthParticle->hasProdVtx()) {
    //   const double p = TMath::Sqrt(TMath::Power(truthParticle->px(),2) +
    // 				   TMath::Power(truthParticle->py(),2) +
    // 				   TMath::Power(truthParticle->pz(),2));

    //   /*truthd0 = (truthParticle->prodVtx()->x()*truthParticle->py() - 
    // 	truthParticle->prodVtx()->y()*truthParticle->px())/p; */ //TODO: check same definition in ParticleAnalysis
      
    //   truthd0 = TMath::Sign(truthParticle->prodVtx()->perp(),
    // 			    truthParticle->prodVtx()->x()*truthParticle->py() - 
    // 			    truthParticle->prodVtx()->y()*truthParticle->px());
      
    //   if (p > 1E-10)
    // 	truthd0 = (truthParticle->prodVtx()->x()*truthParticle->py() - 
    // 		   truthParticle->prodVtx()->y()*truthParticle->px())/p;
    
    float truthz0 = 0.0;
    if (truthParticle->isAvailable<float>("z0")) {
      truthz0 = truthParticle->auxdata<float>("z0");
    }
    else if(truthParticle->hasProdVtx()) {
      truthz0 = truthParticle->prodVtx()->z();
    }

    m_z0_truthz0    -> Fill(trk->z0(), trk->z0() - truthz0);
    m_d0_truthd0    -> Fill(trk->d0(), trk->d0() - truthd0);

    PrintMessage("Truth Properties...");
    m_truthPt       -> Fill(truthParticle->pt()*1e-3,weight);
    m_truthPtnarrow -> Fill(truthParticle->pt()*1e-3,weight);
    m_truthEta      -> Fill(truthParticle->eta(),weight);
    m_truthAbseta   -> Fill(TMath::Abs(truthParticle->eta()),weight);
    m_truthPhi      -> Fill(truthParticle->phi(),weight);
    m_truthD0       -> Fill(truthd0,weight);
    m_truthZ0       -> Fill(truthz0,weight);
    
    float sigPt  = trk->pt()-truthParticle->pt();
    float sigQPt = (trk->charge()/trk->pt()-truthParticle->charge()/truthParticle->pt())*truthParticle->pt();
    float sigEta = trk->eta()-truthParticle->eta();
    float sigPhi = TVector2::Phi_mpi_pi(trk->phi()-truthParticle->phi());
    float sigD0  = trk->d0()-truthd0;
    float sigZ0  = trk->z0()-truthz0;

    PrintMessage("Bias Properties...");
    m_biasPt  -> Fill(sigPt,weight);
    m_biasQPt -> Fill(sigQPt,weight);
    m_biasEta -> Fill(sigEta,weight);
    m_biasPhi -> Fill(sigPhi,weight);
    m_biasD0  -> Fill(sigD0,weight);
    m_biasZ0  -> Fill(sigZ0,weight);
    
    //Bin number of X_vs_eta plots in which this track will fall into
    const int etaBinId = TMath::FloorNint(TMath::Abs(truthParticle->eta())/m_etaMax * double(m_etaVectorSize));

    if (etaBinId < m_etaVectorSize && etaBinId >= 0) { //TODO: check
      Assert("TrackHists::FillHists()\tetaBinId out of bounds", etaBinId < m_etaVectorSize && etaBinId >= 0);
      m_biasPt_abseta[etaBinId]  -> Fill(sigPt,  weight);
      m_biasQPt_abseta[etaBinId] -> Fill(sigQPt, weight);
      m_biasEta_abseta[etaBinId] -> Fill(sigEta/TMath::Abs(truthParticle->eta()), weight); //TODO: sigmaEta/Eta
      m_biasPhi_abseta[etaBinId] -> Fill(sigPhi, weight);
      m_biasD0_abseta[etaBinId]  -> Fill(sigD0,  weight);
      m_biasZ0_abseta[etaBinId]  -> Fill(sigZ0,  weight);

      if (truthParticle->eta() < 0) {
	m_biasPt_negeta[etaBinId]  -> Fill(sigPt,  weight);
	m_biasQPt_negeta[etaBinId]  -> Fill(sigQPt,  weight);
	m_biasEta_negeta[etaBinId]  -> Fill(sigEta,  weight);
	m_biasPhi_negeta[etaBinId]  -> Fill(sigPhi,  weight);
	m_biasD0_negeta[etaBinId]  -> Fill(sigD0,  weight);
	m_biasZ0_negeta[etaBinId]  -> Fill(sigZ0,  weight);
      }
      else {
	m_biasPt_poseta[etaBinId]  -> Fill(sigPt,  weight);
	m_biasQPt_poseta[etaBinId]  -> Fill(sigQPt,  weight);
	m_biasEta_poseta[etaBinId]  -> Fill(sigEta,  weight);
	m_biasPhi_poseta[etaBinId]  -> Fill(sigPhi,  weight);
	m_biasD0_poseta[etaBinId]  -> Fill(sigD0,  weight);
	m_biasZ0_poseta[etaBinId]  -> Fill(sigZ0,  weight);      
      }
    }

    if (trk->isAvailable<float>("matchedDR")) {
      bool passCut = false;
      
      float mindR = trk->auxdata<float>("matchedDR");
      float prob = xAOD::TrackHelper::truthMatchProb(trk);
      m_truthMatchProb -> Fill(prob,weight);
      
      PrintMessage("Before passCut...");
      //  passCut=true;
      //    if (nSiHits>8) { passCut=true; }
      //    if (nSiHits>9) { passCut=true; }
      if (nSiHits>10) { passCut=true; }
      
      if (passCut) {
	if (prob>0.95) { // !!!!!!!!! TODO: check !!!!!!!!!!!!!!
	  m_matchingDR -> Fill(mindR,weight);
	  if (mindR<0.01) { m_truthAbseta_dr001 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.02) { m_truthAbseta_dr002 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.03) { m_truthAbseta_dr003 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.04) { m_truthAbseta_dr004 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.05) { m_truthAbseta_dr005 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.06) { m_truthAbseta_dr006 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.07) { m_truthAbseta_dr007 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.08) { m_truthAbseta_dr008 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	  if (mindR<0.09) { m_truthAbseta_dr009 -> Fill(TMath::Abs(truthParticle->eta()),weight); }
	}
      }
    }
  }

  else {
    std::cout << "TrackHists::FillHists()\t No truth particle retrieved, not filling this track" << std::endl;
  }
  // Cluster study
  // m_IBLCharge           -> Fill(trk->auxdata<float>("IBLCharge"),weight); 
  // m_IBLHitSize          -> Fill(trk->auxdata<int>("IBLHitSize"),weight); 
  // m_IBLHitSizeZ         -> Fill(trk->auxdata<int>("IBLHitSizeZ"),weight); 
  // m_IBLHitSizePhi       -> Fill(trk->auxdata<int>("IBLHitSizePhi"),weight); 
  // m_IBLlocalPhi         -> Fill(trk->auxdata<float>("IBLlocalPhi"),weight);
  // m_IBLlocalTheta       -> Fill(trk->auxdata<float>("IBLlocalTheta"),weight);
  // m_IBLlocalPhi_sizePhi -> Fill(trk->auxdata<float>("IBLlocalPhi"),trk->auxdata<int>("IBLHitSizePhi"),weight);
  // m_IBLlocalTheta_sizeZ -> Fill(trk->auxdata<float>("IBLlocalTheta"),trk->auxdata<int>("IBLHitSizeZ"),weight);

  // m_IBLbiasedResidualX    -> Fill(trk->auxdata<float>("IBLbiasedResidualX"),weight);
  // m_IBLbiasedResidualY    -> Fill(trk->auxdata<float>("IBLbiasedResidualY"),weight);
  // m_IBLunbiasedResidualX  -> Fill(trk->auxdata<float>("IBLunbiasedResidualX"),weight);
  // m_IBLunbiasedResidualY  -> Fill(trk->auxdata<float>("IBLunbiasedResidualY"),weight);
  // m_IBLResidualX          -> Fill(trk->auxdata<float>("IBLResidualX"),weight);
  // m_IBLResidualY          -> Fill(trk->auxdata<float>("IBLResidualY"),weight);

  // m_BLayerCharge           -> Fill(trk->auxdata<float>("BLayerCharge"),weight); 
  // m_BLayerHitSize          -> Fill(trk->auxdata<int>("BLayerHitSize"),weight); 
  // m_BLayerHitSizeZ         -> Fill(trk->auxdata<int>("BLayerHitSizeZ"),weight); 
  // m_BLayerHitSizePhi       -> Fill(trk->auxdata<int>("BLayerHitSizePhi"),weight); 
  // m_BLayerlocalPhi         -> Fill(trk->auxdata<float>("BLayerlocalPhi"),weight);
  // m_BLayerlocalTheta       -> Fill(trk->auxdata<float>("BLayerlocalTheta"),weight);
  // m_BLayerlocalPhi_sizePhi -> Fill(trk->auxdata<float>("BLayerlocalPhi"),trk->auxdata<int>("BLayerHitSizePhi"),weight);
  // m_BLayerlocalTheta_sizeZ -> Fill(trk->auxdata<float>("BLayerlocalTheta"),trk->auxdata<int>("BLayerHitSizeZ"),weight);

  // m_BLayerbiasedResidualX    -> Fill(trk->auxdata<float>("BLayerbiasedResidualX"),weight);
  // m_BLayerbiasedResidualY    -> Fill(trk->auxdata<float>("BLayerbiasedResidualY"),weight);
  // m_BLayerunbiasedResidualX  -> Fill(trk->auxdata<float>("BLayerunbiasedResidualX"),weight);
  // m_BLayerunbiasedResidualY  -> Fill(trk->auxdata<float>("BLayerunbiasedResidualY"),weight);
  // m_BLayerResidualX          -> Fill(trk->auxdata<float>("BLayerResidualX"),weight);
  // m_BLayerResidualY          -> Fill(trk->auxdata<float>("BLayerResidualY"),weight);

  PrintMessage("Closing...");

} // FillHists

