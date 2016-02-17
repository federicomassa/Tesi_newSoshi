#ifndef ParticleAnalysis_TwoTrkEffStudy_H
#define ParticleAnalysis_TwoTrkEffStudy_H

// c++ include(s):
#include <iterator>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <EventLoop/Algorithm.h>
#include "AthLinks/ElementLink.h"

// package include(s):
//#include <ParticleAnalysis/Helix.h>
#include <HistManager/TrackHistManager.h>
#include <HistManager/TruthHists.h>
#include <HistManager/EventHists.h>
#include <HistManager/ClusterHistManager.h>
#include <HistManager/TwoTrackHists.h>

#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"

#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackMeasurementValidationContainer.h"
#endif // not __MAKECINT__

class TwoTrkEffStudy:public EL::Algorithm {
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  xAOD::TEvent *m_event;  //!
  int m_eventCounter;     //!

  // trigger tools member variables
  Trig::TrigDecisionTool *m_trigDecTool; //!
  TrigConf::xAODConfigTool *m_configTool; //!

  GoodRunsListSelectionTool *m_grl; //!
  static char* pileupWeightFile; //!

  ToolHandle<CP::IPileupReweightingTool> m_tool; //!

  // track hists
  TwoTrackHists* twotrkHist_reco_pair;  //!
  TwoTrackHists* twotrkHist_minX;  //!
  TwoTrackHists* twotrkHist_minY;  //!
  TrackHists* trkHist_reco_pair;  //!
  TrackHists* trkHist_reco_pair_leading;  //!
  TrackHists* trkHist_reco_pair_subleading;  //!
  ClusterHistManager* clusHist_reco_pair_leading;  //!
  ClusterHistManager* clusHist_reco_pair_subleading;  //!
  ClusterHists* clusHist_leading_low; //!
  ClusterHists* clusHist_leading_high; //!

  EventHists* eventHist_all; //!

// Everything in the header file  that refers to the xAOD edm needs 
#ifndef __MAKECINT__

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
	
  int  buildPixelModuleID(const xAOD::TrackMeasurementValidation* ) const; 
	
  void addTruthDistancePerLayer(   const xAOD::TrackParticle* tp, std::multimap<int, const xAOD::TrackMeasurementValidation*>& clustersOnAModule) const;  
 
  void countContributingParticles( xAOD::TrackMeasurementValidation* prd );
  void dressLayerHitContent( xAOD::TrackParticleContainer* trks, xAOD::TrackMeasurementValidationContainer* pixClusters, xAOD::TrackMeasurementValidationContainer* sctClusters );

  bool checkQuality(const xAOD::TrackParticle* trk, int outerLayer);

  double deltaPhi(const xAOD::TrackParticle* trk1, const xAOD::TrackParticle* trk2) {
    double dphi = trk1->phi() - trk2->phi();
    if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
    if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
    return dphi;
  }

  double deltaPhi(const xAOD::TrackParticle* trk, const xAOD::Jet* jet) {
    double dphi = trk->phi() - jet->phi();
    if (dphi >= TMath::Pi()) { dphi = 2.0*TMath::Pi()-dphi; }
    if (dphi < -TMath::Pi()) { dphi = 2.0*TMath::Pi()+dphi; }
    return dphi;
  }

  double deltaR(const xAOD::TrackParticle* trk1, const xAOD::TrackParticle* trk2) {
    double dphi = deltaPhi(trk1, trk2);
    double deta = trk1->eta()-trk2->eta();
    double dR = TMath::Sqrt(dphi*dphi+deta*deta);
    return dR;
  }

  double deltaR(const xAOD::TrackParticle* trk, const xAOD::Jet* jet) {
    double dphi = deltaPhi(trk, jet);
    double deta = trk->eta()-jet->eta();
    double dR = TMath::Sqrt(dphi*dphi+deta*deta);
    return dR;
  }

  void findPixelClusterOnLayer(const xAOD::TrackParticle* trk, 
                               xAOD::TrackMeasurementValidationContainer* pixClusters, 
                               int targetLayer, 
                               int *nPenetrateLayer, 
                               const xAOD::TrackStateValidation *msosPenetrate[4], 
                               int *nClusterLayer, 
                               const xAOD::TrackMeasurementValidation *pixPenetrateClus[4]);

//  void clusterPatternClassification(int nRefMSOSLayer, const xAOD::TrackStateValidation *msosRef[4], int nRefClusLayer, const xAOD::TrackMeasurementValidation *pixRefClus[4],
//                                    int nTarMSOSLayer, const xAOD::TrackStateValidation *msosTar[4], int nTarClusLayer, const xAOD::TrackMeasurementValidation *pixTarClus[4],
//                                    int *xclusPatternIndex, double *xminMSOSLocalDistance, double *xminClusLocalDistance, double *xminClusGlobalDistance);


  void clusterPatternClassification(int nRefMSOSLayer, const xAOD::TrackStateValidation *msosRef[4], int nRefClusLayer, const xAOD::TrackMeasurementValidation *pixRefClus[4],
                                    int nTarMSOSLayer, const xAOD::TrackStateValidation *msosTar[4], int nTarClusLayer, const xAOD::TrackMeasurementValidation *pixTarClus[4],
                                    int *xclusPatternIndex, double xminMSOSLocalDistance[2], double xminClusLocalDistance[2], double xminClusGlobalDistance[3]);


#endif // not __MAKECINT__


  // this is a standard constructor
  TwoTrkEffStudy();

  static void SetPileupWeightFile(char* apileupWeightFile) {
    pileupWeightFile = apileupWeightFile;
    return;
  }

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(TwoTrkEffStudy, 1);
};

#endif
