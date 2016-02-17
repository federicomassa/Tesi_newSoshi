#ifndef ParticleAnalysis_TrkEffStudy_H
#define ParticleAnalysis_TrkEffStudy_H

// c++ include(s):
#include <iterator>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <EventLoop/Algorithm.h>
#include "AthLinks/ElementLink.h"

// package include(s):
#include <HistManager/TrackHistManager.h>
#include <HistManager/TruthHists.h>
#include <HistManager/EventHists.h>
#include <HistManager/ClusterHistManager.h>


#ifndef __MAKECINT__
#include "xAODTruth/TruthParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackStateValidationContainer.h"
#include "xAODTracking/TrackMeasurementValidationContainer.h"
#endif // not __MAKECINT__


class TrkEffStudy : public EL::Algorithm
{
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

  // switches
  static bool s_doClusterHists;
  static bool s_doBeforeAfterBLayer;  

  // truth hists
  TruthHists* truHist_mother; //!
  TruthHists* truHist_stable; //!

  // track hists
  TrackHistManager* trkHist_reco;        //!
  TrackHistManager* trkHist_reco_bforLay0;        //!
  TrackHistManager* trkHist_reco_aftrLay0;        //!
  TrackHistManager* trkHist_sisp;        //!
  TrackHistManager* trkHist_sisp_reco;       //!
  TrackHistManager* trkHist_sisp_lost ;       //!
  TrackHistManager* trkHist_sisp_lost0SSCT;        //!  
  TrackHistManager* trkHist_sisp_lostgt0SSCT;        //!   

  TrackHists* trkHist_pseudo_all;             //!
  TrackHists* trkHist_pseudo_SiSp;            //!
  TrackHists* trkHist_pseudo_SiSp_bhit;       //!
  TrackHists* trkHist_pseudo_all_bhit;        //!
  TrackHists* trkHist_pseudo_all_bhitMerge;   //!
  TrackHists* trkHist_pseudo_all_max2ShareSCT;   //!
  TrackHists* trkHist_pseudo_all_max4ShareSCT;   //!
  TrackHists* trkHist_pseudo_reco;            //!
  TrackHists* trkHist_pseudo_reco_has_bhit;   //!
  TrackHists* trkHist_pseudo_reco_has_bhitShare;   //!
  TrackHists* trkHist_pseudo_reco_has_bhitSplit;   //!
  TrackHists* trkHist_pseudo_reco_has_bhitMerge;   //!
  TrackHists* trkHist_pseudo_reco_Lay0Hits;        //!
  TrackHists* trkHist_pseudo_reco_Lay0HitsFound;   //!
  TrackHists* trkHist_pseudo_reco_Lay0HitsExtra;   //!
  TrackHists* trkHist_pseudo_reco_Lay1Hits;        //!
  TrackHists* trkHist_pseudo_reco_Lay1HitsFound;   //!
  TrackHists* trkHist_pseudo_reco_Lay1HitsExtra;   //!
  TrackHists* trkHist_pseudo_reco_Lay2Hits;        //!
  TrackHists* trkHist_pseudo_reco_Lay2HitsFound;   //!
  TrackHists* trkHist_pseudo_reco_Lay2HitsExtra;   //!
  TrackHists* trkHist_pseudo_reco_Lay3Hits;        //!
  TrackHists* trkHist_pseudo_reco_Lay3HitsFound;   //!
  TrackHists* trkHist_pseudo_reco_Lay3HitsExtra;   //!
  TrackHists* trkHist_pseudo_lost;            //!
  TrackHists* trkHist_pseudo_lost_YES_SiSpMatch;   //!      
  TrackHists* trkHist_pseudo_lost_NO_SiSpMatch;    //!  
  TrackHists* trkHist_pseudo_lost0SSCT;            //!       
  TrackHists* trkHist_pseudo_lost0SSCT_YES_SiSpMatch;  //! 
  TrackHists* trkHist_pseudo_lost0SSCT_NO_SiSpMatch;   //!  

  ClusterHistManager* clusHist_lay0;     //!
  ClusterHistManager* clusHist_lay0Found;     //!
  // cluster hists
  /*
  ClusterHistManager* clusHist_all;      //!
  ClusterHistManager* clusHist_lay0;     //!
  ClusterHistManager* clusHist_lay1;     //!
  ClusterHistManager* clusHist_lay2;     //!
  ClusterHistManager* clusHist_lay3;     //!
  */
  std::vector <ClusterHistManager*> pixClusHistsOnPseudoTrack_all; //!
  std::vector <ClusterHistManager*> pixClusHistsOnPseudoTrack_lost; //!  
  std::vector <ClusterHistManager*> pixClusHistsOnPseudoTrack_lost0SSCT; //!
  std::vector <ClusterHistManager*> pixClusHistsOnPseudoTrack_lost0SSCT_YES_SiSpMatch; //!  
  std::vector <ClusterHistManager*> pixClusHistsOnPseudoTrack_lost0SSCT_NO_SiSpMatch; //!  

  EventHists* eventHist_all; //!
  EventHists* eventHist_all_reco_pstrk; //!
  EventHists* eventHist_lost_pstrk; //!  
  EventHists* eventHist_lost0SSCT_pstrk; //! 
  EventHists* eventHist_lostgt0SSCT_pstrk ; //! 

// Everything in the header file  that refers to the xAOD edm needs 
#ifndef __MAKECINT__

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
	
  void filldetIdtoMsosMap( const xAOD::TrackParticleContainer* tc, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const;

  void MinDistanceOtherTrack( const xAOD::TrackParticle*, const xAOD::TrackParticleContainer*, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const;
  
  int  buildPixelModuleID(const xAOD::TrackMeasurementValidation* ) const; 
	
  void addTruthDistancePerLayer(   const xAOD::TrackParticle* tp, std::multimap<int, const xAOD::TrackMeasurementValidation*>& clustersOnAModule) const;  
 
  void fillPixelClusterHists ( const xAOD::TrackStateValidation* msos, xAOD::TrackMeasurementValidationContainer* pixClusters, std::vector<ClusterHistManager*> &clusHists);

  void countContributingParticles( xAOD::TrackMeasurementValidation* prd );
  void dressLayerHitContent( xAOD::TrackParticleContainer* trks, xAOD::TrackMeasurementValidationContainer* pixClusters, xAOD::TrackMeasurementValidationContainer* sctClusters );
  const xAOD::TrackParticle* getHighestProbTrack( std::multimap<int, const xAOD::TrackParticle* >& bcodeMap, int bcode );
  std::vector< const xAOD::TrackMeasurementValidation* > getAllPixelClusters( const std::vector<ElementLink< xAOD::TrackStateValidationContainer > >& mot );


#endif // not __MAKECINT__


  // this is a standard constructor
  TrkEffStudy ();

  // activate/deactivate switches
  static void doClusterHists(bool do_it) {
    s_doClusterHists = do_it;
    return;
  }
  static void doBeforeAfterBLayer(bool do_it) {
    s_doBeforeAfterBLayer = do_it;
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
  ClassDef(TrkEffStudy, 1);
};

#endif
