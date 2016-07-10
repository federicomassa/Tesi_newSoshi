/******************************************
 *
 * Code to analysize simple samples for 
 * study tracking efficiency
 *
 * G. Facini
 * Wed Sep 10 07:21:58 CEST 2014
 *
 * Migrate into rel20. -- Mar  2 11:44:17 CET 2015 S. Tsuno
 *
 ******************************************/

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ParticleAnalysis/TrkEffStudy.h>
#include <HistManager/EventFeaturesInterface.h>
#include <Utility/TrackHelper.h>

// EDM includes: - if move to header file will not compile!
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "TFile.h"

#include <map>
#include <limits>
#include <string>
#include <sstream>

// this is needed to distribute the algorithm to the workers
ClassImp(TrkEffStudy)

// switches - made static because need to be set (in job steering macro) before calling constructor
bool TrkEffStudy :: s_doClusterHists = false;
bool TrkEffStudy :: s_doBeforeAfterBLayer = false;

TrkEffStudy :: TrkEffStudy ()
{

  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  
  truHist_mother = new TruthHists("mother");
  truHist_stable = new TruthHists("stable");

  // In TrackHistManager constructor:
  // 1st bool is to fill plots for tracks matched to primary particles
  // 2nd bool is to fill plots for tracks matched to secondary particles
  // 3rd bool is to fill plots for fake tracks 
  // 4th bool is if to fill duplicates (need 1st and 2nd be true)
  // 5th bool is if to fill plots based on hit content on track (need 1st be true)
  
  // recoTracks histograms
  trkHist_reco   = new TrackHistManager("reco",false,false,false,false,false);
  if( s_doBeforeAfterBLayer ) {
    trkHist_reco_bforLay0   = new TrackHistManager("reco_bforLay0",true,false,false,false,true);
    trkHist_reco_aftrLay0   = new TrackHistManager("reco_aftrLay0",true,false,false,false,true);
  }
  // sisp seeds histograms
  trkHist_sisp              = new TrackHistManager("sisp",false,false,false,false,false);
  trkHist_sisp_reco         = new TrackHistManager("sisp_reco",false,false,false,false,false); 
  trkHist_sisp_lost         = new TrackHistManager("sisp_lost",false,false,false,false,false);
  trkHist_sisp_lost0SSCT    = new TrackHistManager("sisp_lost0SSCT",false,false,false,false,false);
  trkHist_sisp_lostgt0SSCT  = new TrackHistManager("sisp_lostgt0SSCT",false,false,false,false,false);
  // pseudoTracks histograms
  trkHist_pseudo_all           		 = new TrackHists("pseudo_all");
  trkHist_pseudo_all_bhit      		 = new TrackHists("pseudo_all_bhit");
  trkHist_pseudo_all_bhitMerge 		 = new TrackHists("pseudo_all_bhitMerge");
  // trkHist_pseudo_all_max2ShareSCT = new TrackHists("pseudo_all_max2ShareSCT");
  // trkHist_pseudo_all_max4ShareSCT = new TrackHists("pseudo_all_max4ShareSCT");
  trkHist_pseudo_SiSp          		 = new TrackHists("pseudo_SiSp");
  trkHist_pseudo_SiSp_bhit     		 = new TrackHists("pseudo_SiSp_bhit");
  trkHist_pseudo_reco          		 = new TrackHists("pseudo_reco");
  // trkHist_pseudo_reco_has_bhit = new TrackHists("pseudo_reco_has_bhit");
  // trkHist_pseudo_reco_has_bhitShare = new TrackHists("pseudo_reco_has_bhitShare");
  // trkHist_pseudo_reco_has_bhitSplit = new TrackHists("pseudo_reco_has_bhitSplit");
  // trkHist_pseudo_reco_has_bhitMerge = new TrackHists("pseudo_reco_has_bhitMerge");
  trkHist_pseudo_reco_Lay0Hits        = new TrackHists("pseudo_reco_Lay0Hits");
  trkHist_pseudo_reco_Lay0HitsFound   = new TrackHists("pseudo_reco_Lay0HitsFound");
  trkHist_pseudo_reco_Lay0HitsExtra   = new TrackHists("pseudo_reco_Lay0HitsExtra");
  trkHist_pseudo_reco_Lay1Hits        = new TrackHists("pseudo_reco_Lay1Hits");
  trkHist_pseudo_reco_Lay1HitsFound   = new TrackHists("pseudo_reco_Lay1HitsFound");
  trkHist_pseudo_reco_Lay1HitsExtra   = new TrackHists("pseudo_reco_Lay1HitsExtra");
  trkHist_pseudo_reco_Lay2Hits        = new TrackHists("pseudo_reco_Lay2Hits");
  trkHist_pseudo_reco_Lay2HitsFound   = new TrackHists("pseudo_reco_Lay2HitsFound");
  trkHist_pseudo_reco_Lay2HitsExtra   = new TrackHists("pseudo_reco_Lay2HitsExtra");
  trkHist_pseudo_reco_Lay3Hits        = new TrackHists("pseudo_reco_Lay3Hits");
  trkHist_pseudo_reco_Lay3HitsFound   = new TrackHists("pseudo_reco_Lay3HitsFound");
  trkHist_pseudo_reco_Lay3HitsExtra   = new TrackHists("pseudo_reco_Lay3HitsExtra");
  trkHist_pseudo_lost          		 = new TrackHists("pseudo_lost");  
  trkHist_pseudo_lost_YES_SiSpMatch      = new TrackHists("pseudo_lost_YES_SiSpMatch"); 
  trkHist_pseudo_lost_NO_SiSpMatch       = new TrackHists("pseudo_lost_NO_SiSpMatch");  
  trkHist_pseudo_lost0SSCT               = new TrackHists("pseudo_lost0SSCT");  
  trkHist_pseudo_lost0SSCT_YES_SiSpMatch = new TrackHists("pseudo_lost0SSCT_YES_SiSpMatch");  
  trkHist_pseudo_lost0SSCT_NO_SiSpMatch  = new TrackHists("pseudo_lost0SSCT_NO_SiSpMatch");

  clusHist_lay0 = new ClusterHistManager("lay0", 0, true);
  clusHist_lay0Found = new ClusterHistManager("lay0Found", 0, true);

  if(s_doClusterHists) {
    /*
       clusHist_all  = new ClusterHistManager("all", -1);
       clusHist_lay0 = new ClusterHistManager("lay0", 0);
       clusHist_lay1 = new ClusterHistManager("lay1", 1);
       clusHist_lay2 = new ClusterHistManager("lay2", 2);
       clusHist_lay3 = new ClusterHistManager("lay3", 3);
       */
    pixClusHistsOnPseudoTrack_all.reserve(4);
    pixClusHistsOnPseudoTrack_lost0SSCT.reserve(4);
    for (int iLay = 0; iLay < 5; iLay++){
      std::stringstream ssi; ssi << iLay-1;  // first element in vector has label '-1' (meaning "all layers"). Second elemt has label '0' (meaning layer 0) and so on..
      std:: string iLayString = ssi.str();

      ClusterHistManager* pixclushistonpseudotrack_all = new ClusterHistManager("Pixel_OnPseudoTrack_all_lay"+iLayString , iLay-1);
      pixClusHistsOnPseudoTrack_all.push_back(pixclushistonpseudotrack_all);

      ClusterHistManager*  pixclushistonpseudotrack_lost = new ClusterHistManager("Pixel_OnPseudoTrack_lost_lay"+iLayString , iLay-1);
      pixClusHistsOnPseudoTrack_lost.push_back(pixclushistonpseudotrack_lost);      

      ClusterHistManager*  pixclushistonpseudotrack_lost0SSCT = new ClusterHistManager("Pixel_OnPseudoTrack_lost0SSCT_lay"+iLayString , iLay-1);
      pixClusHistsOnPseudoTrack_lost0SSCT.push_back(pixclushistonpseudotrack_lost0SSCT);

      ClusterHistManager*  pixclushistonpseudotrack_lost0SSCT_YES_SiSpMatch = new ClusterHistManager("Pixel_OnPseudoTrack_lost0SSCT_YES_SiSpMatch_lay"+iLayString , iLay-1);
      pixClusHistsOnPseudoTrack_lost0SSCT_YES_SiSpMatch.push_back(pixclushistonpseudotrack_lost0SSCT_YES_SiSpMatch);      
      ClusterHistManager*  pixclushistonpseudotrack_lost0SSCT_NO_SiSpMatch = new ClusterHistManager("Pixel_OnPseudoTrack_lost0SSCT_NO_SiSpMatch_lay"+iLayString , iLay-1);
      pixClusHistsOnPseudoTrack_lost0SSCT_NO_SiSpMatch.push_back(pixclushistonpseudotrack_lost0SSCT_NO_SiSpMatch);      

    }  
  }

  // event histograms
  eventHist_all               = new EventHists("all");
  eventHist_all_reco_pstrk    = new EventHists("all_reco_pstrk");
  eventHist_lost_pstrk        = new EventHists("lost_pstrk");  
  eventHist_lost0SSCT_pstrk   = new EventHists("lost0SSCT_pstrk"); 
  eventHist_lostgt0SSCT_pstrk = new EventHists("lostgt0SSCT_pstrk");
}


EL::StatusCode TrkEffStudy :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD ();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "PartAnaxAOD" ).ignore(); // call before opening first file

  // count number of events
  m_eventCounter = 0;
  std::cout << "Set counter to 0 " << m_eventCounter << std::endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  
  truHist_mother->Init( wk() );
  truHist_stable->Init( wk() );

  trkHist_reco  ->Init( wk() );
  if(s_doBeforeAfterBLayer ) { 
    trkHist_reco_bforLay0  ->Init( wk() );
    trkHist_reco_aftrLay0  ->Init( wk() );
  }
  trkHist_sisp             ->Init( wk() );
  trkHist_sisp_reco        ->Init( wk() );
  trkHist_sisp_lost        ->Init( wk() );
  trkHist_sisp_lost0SSCT   ->Init( wk() );
  trkHist_sisp_lostgt0SSCT ->Init( wk() );

  trkHist_pseudo_all           ->Init( wk() );
  trkHist_pseudo_all_bhit      ->Init( wk() );
  trkHist_pseudo_all_bhitMerge ->Init( wk() ); 
  // trkHist_pseudo_all_max2ShareSCT ->Init( wk() );
  // trkHist_pseudo_all_max4ShareSCT ->Init( wk() );   
  trkHist_pseudo_SiSp         ->Init( wk() );
  trkHist_pseudo_SiSp_bhit    ->Init( wk() );
  trkHist_pseudo_reco         ->Init( wk() );
  // trkHist_pseudo_reco_has_bhit      ->Init( wk() );
  // trkHist_pseudo_reco_has_bhitShare ->Init( wk() );
  // trkHist_pseudo_reco_has_bhitSplit ->Init( wk() );
  // trkHist_pseudo_reco_has_bhitMerge ->Init( wk() );
  trkHist_pseudo_reco_Lay0Hits         ->Init( wk() );
  trkHist_pseudo_reco_Lay0HitsFound    ->Init( wk() );
  trkHist_pseudo_reco_Lay0HitsExtra    ->Init( wk() );
  trkHist_pseudo_reco_Lay1Hits         ->Init( wk() );
  trkHist_pseudo_reco_Lay1HitsFound    ->Init( wk() );
  trkHist_pseudo_reco_Lay1HitsExtra    ->Init( wk() );
  trkHist_pseudo_reco_Lay2Hits         ->Init( wk() );
  trkHist_pseudo_reco_Lay2HitsFound    ->Init( wk() );
  trkHist_pseudo_reco_Lay2HitsExtra    ->Init( wk() );
  trkHist_pseudo_reco_Lay3Hits         ->Init( wk() );
  trkHist_pseudo_reco_Lay3HitsFound    ->Init( wk() );
  trkHist_pseudo_reco_Lay3HitsExtra    ->Init( wk() );
  trkHist_pseudo_lost                    ->Init( wk() );
  trkHist_pseudo_lost_YES_SiSpMatch	 ->Init( wk() );
  trkHist_pseudo_lost_NO_SiSpMatch	 ->Init( wk() );
  trkHist_pseudo_lost0SSCT        	 ->Init( wk() );
  trkHist_pseudo_lost0SSCT_YES_SiSpMatch ->Init( wk() );
  trkHist_pseudo_lost0SSCT_NO_SiSpMatch  ->Init( wk() );

  clusHist_lay0->Init( wk() );
  clusHist_lay0Found->Init( wk() );

  if(s_doClusterHists) {
    /*
    clusHist_all ->Init( wk() );
    clusHist_lay0->Init( wk() );
    clusHist_lay1->Init( wk() );
    clusHist_lay2->Init( wk() );
    clusHist_lay3->Init( wk() );
    */
    for (auto histIt = pixClusHistsOnPseudoTrack_all.begin(); histIt != pixClusHistsOnPseudoTrack_all.end(); ++histIt){
      (*histIt)->Init( wk() );
    }    
    for (auto histIt = pixClusHistsOnPseudoTrack_lost.begin(); histIt != pixClusHistsOnPseudoTrack_lost.end(); ++histIt){
      (*histIt)->Init( wk() );
    }      
    for (auto histIt = pixClusHistsOnPseudoTrack_lost0SSCT.begin(); histIt != pixClusHistsOnPseudoTrack_lost0SSCT.end(); ++histIt){
      (*histIt)->Init( wk() );
    }     
    for (auto histIt = pixClusHistsOnPseudoTrack_lost0SSCT_YES_SiSpMatch.begin(); histIt != pixClusHistsOnPseudoTrack_lost0SSCT_YES_SiSpMatch.end(); ++histIt){
      (*histIt)->Init( wk() );
    }      
    for (auto histIt = pixClusHistsOnPseudoTrack_lost0SSCT_NO_SiSpMatch.begin(); histIt != pixClusHistsOnPseudoTrack_lost0SSCT_NO_SiSpMatch.end(); ++histIt){
      (*histIt)->Init( wk() );
    }        
  }
  
  eventHist_all               ->Init( wk() );
  eventHist_all_reco_pstrk    ->Init( wk() );
  eventHist_lost_pstrk        ->Init( wk() );
  eventHist_lost0SSCT_pstrk   ->Init( wk() ); 
  eventHist_lostgt0SSCT_pstrk ->Init( wk() );

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TrkEffStudy :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: changeInput (bool /* firstFile */)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  TFile* infile = wk()->inputFile ();
  std::cout << infile->GetName() << std::endl;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_event = wk()->xaodEvent();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int
	m_eventCounter=0;
	
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: execute ()
{
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  
  // print every 100 events, so we know where we are:
  if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
  //std::cout << "Event count " << m_eventCounter << std::endl;
  m_eventCounter++;
//  if(m_eventCounter > 1e5) { 
//    if(m_eventCounter == 1e5) { 
//      std::cout << "STOPPING after event " << m_eventCounter << std::endl;
//      std::cout << "Look in here to change -> TrkEffStudy :: execute()" << std::endl;
//    }
//    return EL::StatusCode::SUCCESS; 
//  }

  EventFeatures eventFeatures;

  //---------------------------
  // Reset barcode map in histogram class
  //---------------------------
  trkHist_reco->resetBarcodes();
  if( s_doBeforeAfterBLayer ) {
    trkHist_reco_bforLay0->resetBarcodes();
    trkHist_reco_aftrLay0->resetBarcodes();
  }
  trkHist_sisp->resetBarcodes();
  trkHist_sisp_reco->resetBarcodes();
  trkHist_sisp_lost->resetBarcodes();
  trkHist_sisp_lost0SSCT->resetBarcodes();
  trkHist_sisp_lostgt0SSCT->resetBarcodes();

  trkHist_pseudo_all               ->resetBarcodes();
  trkHist_pseudo_all_bhit          ->resetBarcodes();
  trkHist_pseudo_all_bhitMerge     ->resetBarcodes();
  // trkHist_pseudo_all_max2ShareSCT  ->resetBarcodes();
  // trkHist_pseudo_all_max4ShareSCT  ->resetBarcodes();
  trkHist_pseudo_SiSp              ->resetBarcodes();
  trkHist_pseudo_SiSp_bhit         ->resetBarcodes();
  trkHist_pseudo_reco              ->resetBarcodes();
  // trkHist_pseudo_reco_has_bhit     ->resetBarcodes();
  // trkHist_pseudo_reco_has_bhitShare->resetBarcodes();
  // trkHist_pseudo_reco_has_bhitSplit->resetBarcodes();
  // trkHist_pseudo_reco_has_bhitMerge->resetBarcodes();
  trkHist_pseudo_reco_Lay0Hits         ->resetBarcodes();
  trkHist_pseudo_reco_Lay0HitsFound    ->resetBarcodes();
  trkHist_pseudo_reco_Lay1Hits         ->resetBarcodes();
  trkHist_pseudo_reco_Lay1HitsFound    ->resetBarcodes();
  trkHist_pseudo_reco_Lay2Hits         ->resetBarcodes();
  trkHist_pseudo_reco_Lay2HitsFound    ->resetBarcodes();
  trkHist_pseudo_reco_Lay3Hits         ->resetBarcodes();
  trkHist_pseudo_reco_Lay3HitsFound    ->resetBarcodes();
  trkHist_pseudo_lost                    ->resetBarcodes();
  trkHist_pseudo_lost_YES_SiSpMatch	 ->resetBarcodes();
  trkHist_pseudo_lost_NO_SiSpMatch	 ->resetBarcodes();
  trkHist_pseudo_lost0SSCT	         ->resetBarcodes();
  trkHist_pseudo_lost0SSCT_YES_SiSpMatch ->resetBarcodes();
  trkHist_pseudo_lost0SSCT_NO_SiSpMatch  ->resetBarcodes();

  //---------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  xAOD::TStore store; // needed for shallow copies

  bool debugInitial(false);
  bool debug(false);
  bool debugMarco(false);
  if(debug) { Info("execute()", "Start Event"); }

//  // check if the event is data or MC
//  // (many tools are applied either to data or MC)
//  bool isMC = false;
//  // check if the event is MC
//  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
//    isMC = true; // can do something with this later
//  } 

  //---------------------
  // TRUTH PARTICLES
  //---------------------
  // get track container of interest
  const xAOD::TruthParticleContainer* truthPartsOrig = 0;
  if ( !m_event->retrieve( truthPartsOrig, "TruthParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Truth Particle container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //---------------------
  // PSEUDO TRACKS
  //---------------------
  // get track container of interest
  const xAOD::TrackParticleContainer* pseudoTracksOrig = 0;
  if ( !m_event->retrieve( pseudoTracksOrig, "InDetPseudoTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Pseduo Track container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventFeatures.nPseudoTracks = pseudoTracksOrig->size(); 

  // loop over the tracks in the container 
  if( debugInitial ) {
    for(xAOD::TrackParticleContainer::const_iterator pseudoTrk_itr = pseudoTracksOrig->begin(); 
        pseudoTrk_itr != pseudoTracksOrig->end(); pseudoTrk_itr++) {
      Info("execute()", "  original pseudo track pT = %.2f GeV", ((*pseudoTrk_itr)->pt() * 0.001));
      const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(*pseudoTrk_itr);
      if( truthParticle ) {
        Info("execute()", "   barcode %i", truthParticle->barcode());
      }
    } // end for loop over pseudo tracks
  } // debugInitial


  //---------------------
  // RECONSTRUCTED TRACKS
  //---------------------
  // get track container of interest
  const xAOD::TrackParticleContainer* recoTracksOrig = 0;
  if ( !m_event->retrieve( recoTracksOrig, "InDetTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Reconstructed Track container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventFeatures.nRecoTracks = recoTracksOrig->size(); 
 
  // loop over the tracks in the container
  if( debugInitial ) {
    for(xAOD::TrackParticleContainer::const_iterator recoTrk_itr = recoTracksOrig->begin();
        recoTrk_itr != recoTracksOrig->end(); recoTrk_itr++) {
      Info("execute()", "  original reco track pT = %.2f GeV", ((*recoTrk_itr)->pt() * 0.001));
      const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(*recoTrk_itr);
      if( truthParticle ) {
        Info("execute()", "   barcode %i", truthParticle->barcode());
      }
    } // end for loop over reconstructed tracks
  } // debugInitial

  // make a shallow copy
  std::pair< xAOD::TrackParticleContainer*, xAOD::ShallowAuxContainer* > recoTracksShallowCopyPair = xAOD::shallowCopyContainer( *recoTracksOrig );
  xAOD::TrackParticleContainer* recoTracks = recoTracksShallowCopyPair.first;

  //---------------------
  // SiSpSeededTracks TRACKS
  //---------------------
  typedef std::vector<ElementLink< xAOD::TrackParticleContainer > > SiSpSeededTracks;
  typedef std::vector<ElementLink< xAOD::TrackParticleContainer > >::const_iterator SiSpSeededTracksIter;
  // get track container of interest
  const xAOD::TrackParticleContainer* sispTrksOrig = 0;
  if ( !m_event->retrieve( sispTrksOrig, "SiSpSeededTrackParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve SiSpSeededTrack container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventFeatures.nSiSpSeeds = sispTrksOrig->size();

  // make a shallow copy
  std::pair< xAOD::TrackParticleContainer*, xAOD::ShallowAuxContainer* > sispTrksShallowCopyPair = xAOD::shallowCopyContainer( *sispTrksOrig );
  xAOD::TrackParticleContainer* sispTracks = sispTrksShallowCopyPair.first;

  //---------------------
  // Pixel Clusters
  //---------------------
  // get cluster container of interest
  const xAOD::TrackMeasurementValidationContainer* pixClustersOrig = 0;
  if ( !m_event->retrieve( pixClustersOrig, "PixelClusters" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve Pixel Cluster container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventFeatures.nPixelClusters = pixClustersOrig->size();
  
  // loop over the clusters in the container 
  if(debugInitial) {
    xAOD::TrackMeasurementValidationContainer::const_iterator pixCluster_itr = pixClustersOrig->begin();
    xAOD::TrackMeasurementValidationContainer::const_iterator pixCluster_end = pixClustersOrig->end();
    for( ; pixCluster_itr != pixCluster_end; ++pixCluster_itr ) {
      Info("execute()", "  pixel cluster sizeZ   = %.2i ", ((*pixCluster_itr)->auxdataConst< int >( "sizeZ" )));
      Info("execute()", "  pixel cluster sizePhi = %.2i ", ((*pixCluster_itr)->auxdataConst< int >( "sizePhi" )));
      Info("execute()", "  pixel cluster size    = %.2i ", ((*pixCluster_itr)->auxdataConst< int >( "size" )));      
    } // end for loop over clusters
  } // initial debug

  // make a shallow copy
  std::pair< xAOD::TrackMeasurementValidationContainer*, xAOD::ShallowAuxContainer* > pixClustersPair = xAOD::shallowCopyContainer( *pixClustersOrig );
  xAOD::TrackMeasurementValidationContainer* pixClusters =  pixClustersPair.first;

  //---------------------
  // SCT Clusters
  //---------------------
  // get cluster container of interest
  const xAOD::TrackMeasurementValidationContainer* sctClustersOrig = 0;
  if ( !m_event->retrieve( sctClustersOrig, "SCT_Clusters" ).isSuccess() ){ // retrieve arguments: container type, container key
    Error("execute()", "Failed to retrieve SCT Cluster container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventFeatures.nSCTClusters = sctClustersOrig->size();
  
  // now count total number of Si clusters in event
  eventFeatures.nSiClusters = pixClustersOrig->size() + sctClustersOrig->size();

  // make a shallow copy
  std::pair< xAOD::TrackMeasurementValidationContainer*, xAOD::ShallowAuxContainer* > sctClustersPair = xAOD::shallowCopyContainer( *sctClustersOrig );
  xAOD::TrackMeasurementValidationContainer* sctClusters =  sctClustersPair.first;
  
  //---------------------
  // PUT IN WORK
  //---------------------

  //---------------------
  // Build a map of every cluster on a module
  //---------------------  
  std::multimap<int, const xAOD::TrackMeasurementValidation*> clusterPerModule;
  for( auto cluster: *pixClusters ) {
    clusterPerModule.insert( std::pair< int, const xAOD::TrackMeasurementValidation*>( buildPixelModuleID(cluster), cluster));
  } // end for loop over clusters

  //---------------------
  // FIND THE MOTHER
  //---------------------
  // Find the mother and fill a container with stable primary particles
  // keep only primary stable particles in this container
  const xAOD::TruthParticle* mother = 0;
  float motherPt(-1), motherR(-1);
  //Create the new container and its auxiliary store for only primary truth particles 
  //need this since link to mother and children is not yet in the AR...
  ConstDataVector< xAOD::TruthParticleContainer > primaryTruth(SG::VIEW_ELEMENTS);
  //xAOD::TruthParticleContainer* truthPriParts = new xAOD::TruthParticleContainer();
  //xAOD::TruthParticleAuxContainer* truthPriPartsAux = new xAOD::TruthParticleAuxContainer();
  //truthPriParts->setStore( truthPriPartsAux ); //< Connect the two
  for(xAOD::TruthParticleContainer::const_iterator truthPart_itr = truthPartsOrig->begin(); 
      truthPart_itr != truthPartsOrig->end(); truthPart_itr++) {
    if((*truthPart_itr)->barcode() <= 10e3 || (*truthPart_itr)->barcode() >= 200e3) { continue; }
    if((*truthPart_itr)->status() == 999) { // the mother
      motherPt = (*truthPart_itr)->pt();
      if( (*truthPart_itr)->hasDecayVtx() ) { motherR = (*truthPart_itr)->decayVtx()->perp(); }
      mother = (*truthPart_itr);
      continue;
    }
    
    if( !(xAOD::TrackHelper::isPrimary( *truthPart_itr)) ) continue;

    // kill particles that decay to themselves
    // for particles that oscillate - need to have an extra line
    //   if daughter is same pdgid with different sign kill particle
    //   if parent is same pdgid with different sign - go

    // Kill particle if has same pdgid as parent ( no abs )
    int pdgid((*truthPart_itr)->pdgId());
    bool selfDecay(false);
    if( (*truthPart_itr)->hasProdVtx() ) { 
      size_t nParents = (*truthPart_itr)->prodVtx()->nIncomingParticles();
      for( size_t iparent=0; iparent < nParents; iparent++ ) {
        const xAOD::TruthParticle* parent = (*truthPart_itr)->prodVtx()->incomingParticle(iparent);
        if( parent->pdgId() == pdgid ) { selfDecay = true; break; }  // do not check sign
      }
    } // has prodVtx
    if(selfDecay) { 
      std::cout << " self decay " << (*truthPart_itr)->barcode() << "\t" << (*truthPart_itr)->pt() << std::endl;
      continue; 
    }

    // Kill particle if it decays to itself but with a different sign - oscillation!
    // The daughter will not be discarded when it comes through
//    bool oscillation(false);
//    int absPdgId( abs(pdgid) );
//    for( size_t ichild=0; ichild < (*truthPart_itr)->nChildren(); ichild++ ) {
//      const xAOD::TruthParticle* child = (*truthPart_itr)->child(ichild);
//      if( child->absPdgId() == absPdgId && child->pdgId() != pdgid ) { 
//        oscillation = true; 
//        break;
//      }
//    }
    // daughters are not that clean....especially for pions...maybe okay for b's
//    if(oscillation) { 
//      std::cout << " oscillation " << (*truthPart_itr)->barcode() << "\t" << (*truthPart_itr)->pt() << std::endl;
//      continue; 
//    }

    // Copy this particle to the new container:
    //xAOD::TruthParticle* tp = new xAOD::TruthParticle();
    //tp->makePrivateStore( **truthPart_itr );
    //truthPriParts->push_back( tp );
    primaryTruth.push_back( *truthPart_itr );
  }
  const xAOD::TruthParticleContainer* truthPriParts = primaryTruth.asDataVector();

  if( !mother ) {
    Info("execute()", "NO MOTHER");
    return EL::StatusCode::SUCCESS;
  }
 
  //---------------------
  // Prune Pseudotracks
  //---------------------

  // Only keep the pseudotracks from stable primary particles
  //Create the new container and its auxiliary store for pseudotracks
  xAOD::TrackParticleContainer* pseudoTracks = new xAOD::TrackParticleContainer();
  xAOD::TrackParticleAuxContainer* pseudoTracksAux = new xAOD::TrackParticleAuxContainer();
  pseudoTracks->setStore( pseudoTracksAux ); //< Connect the two
  unsigned int n_PrimaryPseudoTracks(0), n_SecondaryPseudoTracks(0);
  for(xAOD::TrackParticleContainer::const_iterator pseudoTrk_itr = pseudoTracksOrig->begin(); 
    pseudoTrk_itr != pseudoTracksOrig->end(); pseudoTrk_itr++) {
    const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(*pseudoTrk_itr);
    if( !truthParticle ) { 
      if(debug) { 
        Info("execute()", "No truth link for pseudotrack?"); 
        std::cout << "  truthMatchProbability " << (*pseudoTrk_itr)->auxdata< float >( "truthMatchProbability" ) << std::endl;
      }
      continue;  // kill pseudotracks without a link ( and thus no barcode? )
    }

    if (xAOD::TrackHelper::isPrimary(truthParticle))        {n_PrimaryPseudoTracks++;}
    else if (xAOD::TrackHelper::isSecondary(truthParticle)) {n_SecondaryPseudoTracks++;}

    if (xAOD::TrackHelper::isPrimary(truthParticle))
    {
      // Copy this track to the new container:
      xAOD::TrackParticle* tp = new xAOD::TrackParticle();
      tp->makePrivateStore( **pseudoTrk_itr);
      pseudoTracks->push_back( tp );
      addTruthDistancePerLayer( tp, clusterPerModule );
    }
  }

  eventFeatures.nPrimaryPseudoTracks   = n_PrimaryPseudoTracks; 
  eventFeatures.nSecondaryPseudoTracks = n_SecondaryPseudoTracks;

  // only look at events where each stable particle has a corresponding pseudotrack
  // if ( pseudoTracks->size() < truthPriParts->size() ) { return EL::StatusCode::SUCCESS; }

  // get cluster on pseudo tracks -> do not want noise clusters and all clusters
  // from primary or secondary particles will be on a pseudotrack
  // kill clusters not containing hits from a primary particle?
  const xAOD::TrackStateValidationContainer* pseudoPixelMSOS = 0;
  if ( !m_event->retrieve( pseudoPixelMSOS, "PseudoPixelMSOSs" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Pseudotrack Pixel TrackStateValidationContainer. Exiting." );
    return EL::StatusCode::FAILURE;
  }

//std::cout << "executing 7" << std::endl;

  // dress tracks with number of particles contributing to hits on each layer
  dressLayerHitContent( pseudoTracks, pixClusters, sctClusters );

//std::cout << "executing 7-0" << std::endl;

  dressLayerHitContent( recoTracks, pixClusters, sctClusters); 

//std::cout << "executing 7-1" << std::endl;

  dressLayerHitContent( sispTracks, pixClusters, sctClusters); 

//std::cout << "executing 8" << std::endl; 

  // plot clusters but only plot 1 sibling for split clusters
  for( xAOD::TrackStateValidationContainer::const_iterator msos_itr = pseudoPixelMSOS->begin();
          msos_itr != pseudoPixelMSOS->end(); msos_itr++) {

     // The detector type enum is defined in TrkEventPrimitives/TrackStateDefs.h
    if( (*msos_itr)->detType() != 1 ) { continue; }  // confirm this is a pixel hit
     // The measurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
    if( (*msos_itr)->type()    != 0 ) { continue; }  // confirm this is a measurement type

    // check if already plotting this MSOS's PRD's sibling
      // siblings are copies of the parent with a new local position
      // the identifier and positions can be different but the layer, size, and charge information will be the same
    bool skip(false);
    for( xAOD::TrackStateValidationContainer::const_iterator msos_jtr = pseudoPixelMSOS->begin(); msos_jtr != msos_itr; ++msos_jtr) {
      if ( xAOD::TrackHelper::areSiblings( *msos_itr, *msos_jtr) ) { skip = true; break; }
    }
    if(skip) { continue; }

    // fill pixel cluster plots 
    fillPixelClusterHists( *msos_itr, pixClusters, pixClusHistsOnPseudoTrack_all);
    
  } // loop over pseudo track pixel TrackStateValidationContainer 
  
  std::multimap<uint64_t, const xAOD::TrackStateValidation*> detIdtoMsosMap_pseudo;
  filldetIdtoMsosMap(pseudoTracks, detIdtoMsosMap_pseudo);
  std::multimap<uint64_t, const xAOD::TrackStateValidation*> detIdtoMsosMap_sisp;
  filldetIdtoMsosMap(sispTracks, detIdtoMsosMap_sisp);

//std::cout << "executing 9" << std::endl; 

  //---------------------
  // FILL PLOTS!!!
  //---------------------

  //eventFeatures.dumpInfo();
  eventHist_all->FillHists( eventFeatures, 1.0 ); 

  truHist_mother->FillHists( mother, 1.0 );

  xAOD::TruthParticleContainer::const_iterator truth_itr = truthPriParts->begin();
  xAOD::TruthParticleContainer::const_iterator truth_end = truthPriParts->end();
  if(debug) { std::cout << "Truth : mother pT and R " << motherPt*0.001 << " " << motherR << std::endl; }
  for( ; truth_itr != truth_end; ++truth_itr ) {
    (*truth_itr)->auxdecor< float >("motherPt") = motherPt;
    (*truth_itr)->auxdecor< float >("motherR")  = motherR;
    truHist_stable->FillHists( (*truth_itr), 1.0 );
    if(debug) { std::cout << "Truth : " << (*truth_itr)->barcode() << " " << (*truth_itr)->pt()*0.001 << std::endl; }
  }

  typedef std::multimap<int, const xAOD::TrackParticle* > barcodeMap;
  typedef std::multimap<int, const xAOD::TrackParticle* >::iterator barcodeMapIter;
  barcodeMap recoTracksBarcodeMap;
  barcodeMap sispTracksBarcodeMap;

  // reco tracks
  xAOD::TrackParticleContainer::iterator reco_itr = (recoTracks)->begin();
  xAOD::TrackParticleContainer::iterator reco_end = (recoTracks)->end();
  for( ; reco_itr != reco_end; ++reco_itr ) {
    
    addTruthDistancePerLayer( *reco_itr, clusterPerModule );
  
    const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(*reco_itr);

    if( truthParticle ) {
      recoTracksBarcodeMap.insert( std::pair<int, const xAOD::TrackParticle*> (truthParticle->barcode(), *reco_itr) ); 
      // make and easier to access barcode and set mother properties
      (*reco_itr)->auxdata< int >("barcode") = truthParticle->barcode();
    } else {
      if(debug) { Info("execute()", "no barcode for reco track"); }
      (*reco_itr)->auxdata< int >("barcode") = -888;
    }
    (*reco_itr)->auxdata< float >("motherPt") = motherPt;
    (*reco_itr)->auxdata< float >("motherR")  = motherR;
    
    trkHist_reco->FillHists( (*reco_itr), 1.0 ); 
    if( s_doBeforeAfterBLayer && truthParticle ) {
      if( truthParticle->hasProdVtx() ) {
        if( truthParticle->prodVtx()->perp() < 49 ) {
          trkHist_reco_bforLay0->FillHists( (*reco_itr), 1.0 ); 
        } else if( truthParticle->prodVtx()->perp() > 52 ) {
          trkHist_reco_aftrLay0->FillHists( (*reco_itr), 1.0 ); 
        }
      }
    }
  } // loop over reco tracks

  // SiSp tracks
  for( auto sispTrack: *sispTracks ) {
  
    addTruthDistancePerLayer( sispTrack, clusterPerModule );

    const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(sispTrack);

    if( truthParticle ) {
      // make  easier to access barcode and set mother properties
      sispTracksBarcodeMap.insert( std::pair<int, const xAOD::TrackParticle*> (truthParticle->barcode(), sispTrack) ); 
      sispTrack->auxdata< int >("barcode") = truthParticle->barcode();
    } else {
      if(debug) { Info("execute()", "no barcode for reco track"); }
      sispTrack->auxdata< int >("barcode") = -888;
    }
    sispTrack->auxdata< float >("motherPt") = motherPt;
    sispTrack->auxdata< float >("motherR")  = motherR;
    // set min dist for all pix layers inside the auxdata, then will retrieve it in TrackHist   

    MinDistanceOtherTrack(sispTrack, sispTracks, detIdtoMsosMap_sisp); 

    trkHist_sisp->FillHists( sispTrack, 1.0 );

  } // loop over sisp tracks

  // match pseudos to reco, sisp
  xAOD::TrackParticleContainer::iterator pseudo_itr = pseudoTracks->begin();
  xAOD::TrackParticleContainer::iterator pseudo_end = pseudoTracks->end();
  
  bool evtHasLost = false;
  bool evtHasLost0SSCT = false;
  
  for( ; pseudo_itr != pseudo_end; ++pseudo_itr ) {

    // make and easier to access barcode and set mother properties
    const xAOD::TruthParticle *truthParticle = xAOD::TrackHelper::truthParticle(*pseudo_itr);
 
    if( truthParticle ) {
      (*pseudo_itr)->auxdata< int >("barcode") = truthParticle->barcode();
    } else {
      if(debug) { Info("execute()", "no barcode for pseudo track %.2f", ((*pseudo_itr)->pt() * 0.001)); }
      (*pseudo_itr)->auxdata< int >("barcode") = -999;
    }
    (*pseudo_itr)->auxdata< float >("motherPt") = motherPt;
    (*pseudo_itr)->auxdata< float >("motherR")  = motherR;
    
    // set min dist for all pix layers inside the auxdata, then will retrieve it in TrackHist
    MinDistanceOtherTrack(*pseudo_itr, pseudoTracks, detIdtoMsosMap_pseudo); 

    // count matches - hopefully only one
    int nMatchReco = recoTracksBarcodeMap.count( (*pseudo_itr)->auxdata< int >("barcode") );
    int nMatchSisp = sispTracksBarcodeMap.count( (*pseudo_itr)->auxdata< int >("barcode") );
    (*pseudo_itr)->auxdata< int >( "nMatchReco" ) = nMatchReco;
    (*pseudo_itr)->auxdata< int >( "nMatchSisp" ) = nMatchSisp;

     // get the match in sisp container with the highest truth prob score
    const xAOD::TrackParticle* theSiSpMatch = getHighestProbTrack( sispTracksBarcodeMap, (*pseudo_itr)->auxdata< int >("barcode") );
    if(debug && theSiSpMatch){
       std::cout << "\t pseudoTrack with barcode: " << (*pseudo_itr)->auxdata< int >("barcode") << 
                    " is best-matched with siSpSeeded Track with barcode: " << theSiSpMatch->auxdata< int >("barcode") << 
		    " and truth match probability: " << theSiSpMatch->auxdata< float >( "truthMatchProbability" )  
		    << std::endl;
    }
    if (theSiSpMatch) xAOD::TrackHelper::countMissingAndWrongHitsOnTrack( theSiSpMatch, (*pseudo_itr) );

    int nMatchSispBHit = 0;
    auto range = sispTracksBarcodeMap.equal_range( (*pseudo_itr)->auxdata< int >("barcode") );
    uint8_t getInt(0);   // for accessing summary information
    for( auto sispIter=range.first; sispIter!=range.second; ++sispIter ){
      //has a b-layer hit?
      ((*sispIter).second)->summaryValue( getInt, xAOD::numberOfBLayerHits );
      if( getInt > 0 ) { nMatchSispBHit++; }
      break;
    }
    (*pseudo_itr)->auxdata< int >( "nMatchSispBHit" ) = nMatchSispBHit;

    //fill hists!
    trkHist_pseudo_all->FillHists( (*pseudo_itr), 1.0 );
    if( nMatchSisp > 0 ) { 
      trkHist_pseudo_SiSp->FillHists( (*pseudo_itr), 1.0 );
      if( nMatchSispBHit > 0 ) { 
        trkHist_pseudo_SiSp_bhit->FillHists( (*pseudo_itr), 1.0 );
      }
    }
    // has a b-layer hit?
    (*pseudo_itr)->summaryValue( getInt, xAOD::numberOfBLayerHits ); 
    
    if( getInt > 0 ) { // has one! 
      trkHist_pseudo_all_bhit     ->FillHists( (*pseudo_itr), 1.0 );
    }
    // is b-layer hit merged? Need to go to cluster and count!!
    // already done above for all hits on track!
    if( ((*pseudo_itr)->auxdata< std::vector<int> >( "pixClusPriPartContent" )).at(0) > 1 ) {
      trkHist_pseudo_all_bhitMerge->FillHists( (*pseudo_itr), 1.0 );
    }

//    (*pseudo_itr)->summaryValue( getInt, xAOD::numberOfSCTSharedHits ); 
//    if( getInt < 3 ) { // max 2 shared SCT Hits
//      trkHist_pseudo_all_max2ShareSCT->FillHists( (*pseudo_itr), 1.0 );
//    }
//    if( getInt < 5 ) { // max 2 shared SCT Hits
//      trkHist_pseudo_all_max4ShareSCT->FillHists( (*pseudo_itr), 1.0 );
//    }

    //if( nMatchSisp == 0 ) { }
    /*
    if( nMatchReco == 0 && nMatchSisp == 0) {
      trkHist_pseudo_lost_NO_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
      (*pseudo_itr)->summaryValue( getInt, xAOD::numberOfSCTSharedHits ); 
      if( getInt == 0 ) {
        trkHist_pseudo_lost0SSCT_NO_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
      }
    }
    */
    // check lost pseudotracks
    if( nMatchReco == 0 ) {
      
      trkHist_pseudo_lost->FillHists( (*pseudo_itr), 1.0 );
      evtHasLost = true;
      
      // fill histograms for best-matching sisp track
      if(theSiSpMatch) trkHist_sisp_lost ->FillHists( theSiSpMatch, 1.0 );
      
      (*pseudo_itr)->summaryValue( getInt, xAOD::numberOfSCTSharedHits );      
       
      if( getInt == 0 ) {
        trkHist_pseudo_lost0SSCT->FillHists( (*pseudo_itr), 1.0 );
	evtHasLost0SSCT = true;
      }
      
      if(nMatchSisp == 0){  // lost with no match in SiSp
        trkHist_pseudo_lost_NO_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
        if( getInt == 0 ) {
          trkHist_pseudo_lost0SSCT_NO_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
        }
      } else {  // lost with at least one match in SiSp
        trkHist_pseudo_lost_YES_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
        if( getInt == 0 ) {

	  trkHist_pseudo_lost0SSCT_YES_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
          
	  // fill histograms for best-matching sisp track
          if(theSiSpMatch) trkHist_sisp_lost0SSCT ->FillHists( theSiSpMatch, 1.0 );
          
//
//          trkHist_pseudo_lost0SSCT_YES_SiSpMatch->FillHists( (*pseudo_itr), 1.0 );
//          if(debug) { Info("execute()", "Iterating on sispTracksBarcodeMap. \n Pseudotrack has barcode %i ", (*pseudo_itr)->auxdata< int >("barcode") ); }

//          // get all sisp tracks with the same barcode
//          auto range = sispTracksBarcodeMap.equal_range( (*pseudo_itr)->auxdata< int >("barcode") );
//          for( auto sispIter=range.first; sispIter!=range.second; ++sispIter ){
//            trkHist_sisp_lost0SSCT ->FillHists( (*sispIter).second, 1.0 );
//          }
//
        } else {	  
	  // fill histograms for best-matching sisp track
          if(theSiSpMatch) trkHist_sisp_lostgt0SSCT ->FillHists( theSiSpMatch, 1.0 );
        }
      }

    } else { // nMatchReco != 1

      trkHist_pseudo_reco->FillHists( (*pseudo_itr), 1.0 );

      if(theSiSpMatch) trkHist_sisp_reco ->FillHists( theSiSpMatch, 1.0 );

      // get the match with the highest truth prob score
//      const xAOD::TrackParticle* theMatch = getHighestProbTrack( recoTracksBarcodeMap, (*pseudo_itr)->auxdata< int >("barcode") );
//      // check if matched reco track has a b-layer hit
//      if(theMatch) theMatch->summaryValue( getInt, xAOD::numberOfBLayerHits ); 
//      if( getInt > 0 ) { 
//        trkHist_pseudo_reco_has_bhit     ->FillHists( (*pseudo_itr), 1.0 );
//        // has a shared b-layer hit
//        theMatch->summaryValue( getInt, xAOD::numberOfBLayerSharedHits); 
//        if( getInt > 0 ) { trkHist_pseudo_reco_has_bhitShare->FillHists( *pseudo_itr, 1.0 ); }
//        // has a split b-layer hit
//        theMatch->summaryValue( getInt, xAOD::numberOfBLayerSplitHits); 
//        if( getInt > 0 ) { trkHist_pseudo_reco_has_bhitSplit->FillHists( *pseudo_itr, 1.0 ); }
//        // is b-layer hit merged? Need to go to cluster and count!!
//        // already done above for all hits on track!
//        if( (theMatch->auxdata< std::vector<int> >( "pixClusPriPartContent" )).at(0) > 1 ) {
//          trkHist_pseudo_reco_has_bhitMerge->FillHists( *pseudo_itr, 1.0 );
//        }
//      }

      static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change
      static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack(measurementNames);

      const xAOD::TrackParticle* theMatch = getHighestProbTrack( recoTracksBarcodeMap, (*pseudo_itr)->auxdata< int >("barcode") );


      // compare hit by hit
      if( (*pseudo_itr)->isAvailable< MeasurementsOnTrack >( measurementNames ) 
          && theMatch 
          && (theMatch)->isAvailable< MeasurementsOnTrack >( measurementNames )
          ) { 
        const std::vector<ElementLink< xAOD::TrackStateValidationContainer > >& psdoMOT = acc_MeasurementsOnTrack( *(*pseudo_itr) );
        const std::vector<ElementLink< xAOD::TrackStateValidationContainer > >& recoMOT = acc_MeasurementsOnTrack( *(theMatch) );
        std::vector< const xAOD::TrackMeasurementValidation* > psdoClusters = getAllPixelClusters( psdoMOT );
        std::vector< const xAOD::TrackMeasurementValidation* > recoClusters = getAllPixelClusters( recoMOT );

        // make map of IDs and layers to search
        std::map < int, int > rIDs;
        for( auto cl : recoClusters ) { 
          rIDs[ (cl)->identifier() ] = cl->auxdata<int>("layer");
        }

        std::map < int, int > pIDs;
        bool tskip(false);
        bool tstop(false); // lazy
        for( auto cl : psdoClusters ) { 
          pIDs[ (cl)->identifier() ] = cl->auxdata<int>("layer");
          // fill hists for each cluster
          switch ( cl->auxdata<int>("layer") )
          { 
            case 0 :
              trkHist_pseudo_reco_Lay0Hits -> FillHists( (*pseudo_itr), 1.0 );
              trkHist_pseudo_reco_Lay0Hits ->resetBarcodes();
              tskip = false;
              tstop = false;
              // loop over clusters already filled and make sure not filling a sibling
              for( auto clFilled : psdoClusters ) { 
                if( clFilled == cl ) { tstop = true; }
                if(tstop) { continue; }
                if( xAOD::TrackHelper::areSiblings( cl, clFilled ) ) { tskip = true; tstop = true; }
              }
              if( !tskip ) { 
                tstop = false;
                xAOD::TrackMeasurementValidationContainer::iterator clus_itr = (pixClusters)->begin();
                xAOD::TrackMeasurementValidationContainer::iterator clus_end = (pixClusters)->end();
                for( ; clus_itr != clus_end; ++clus_itr ) {
                  if(tstop) { continue; }
                  if( (*clus_itr)->identifier() != (cl)->identifier() ) { continue; }
                  if( (*clus_itr)->auxdata< float >( "charge" ) != (cl)->auxdata< float >( "charge" ) ) { continue; }
                  clusHist_lay0 -> FillHists( *clus_itr, 1.0 );
                  if( rIDs.count( (cl)->identifier() ) > 0 ) {
                    clusHist_lay0Found -> FillHists( *clus_itr, 1.0 );
                  }
                  tstop = true;
                }
              }
              break;
            case 1 :
              trkHist_pseudo_reco_Lay1Hits -> FillHists( (*pseudo_itr), 1.0 );
              trkHist_pseudo_reco_Lay1Hits ->resetBarcodes();
              break;
            case 2 :
              trkHist_pseudo_reco_Lay2Hits -> FillHists( (*pseudo_itr), 1.0 );
              trkHist_pseudo_reco_Lay2Hits ->resetBarcodes();
              break;
            case 3 :
              trkHist_pseudo_reco_Lay3Hits -> FillHists( (*pseudo_itr), 1.0 );
              trkHist_pseudo_reco_Lay3Hits ->resetBarcodes();
              break;
            default :
              break;
          }
        }


        // loop over psuedo and look for match in reco collection
        for( auto pC : pIDs ) {
          if( rIDs.count( pC.first ) == 0 ) {
            // fill hist for missing clusters
//            switch ( pC.second ) // layer
//            { 
//              case 0 :
//                trkHist_pseudo_reco_Lay0HitsLost -> FillHists( (*pseudo_itr), 1.0 );
//                trkHist_pseudo_reco_Lay0HitsLost ->resetBarcodes();
//                break;
//              case 1 :
//                trkHist_pseudo_reco_Lay1HitsLost -> FillHists( (*pseudo_itr), 1.0 );
//                trkHist_pseudo_reco_Lay1HitsLost ->resetBarcodes();
//                break;
//              case 2:
//                trkHist_pseudo_reco_Lay2HitsLost -> FillHists( (*pseudo_itr), 1.0 );
//                trkHist_pseudo_reco_Lay2HitsLost ->resetBarcodes();
//                break;
//              case 3:
//                trkHist_pseudo_reco_Lay3HitsLost -> FillHists( (*pseudo_itr), 1.0 );
//                trkHist_pseudo_reco_Lay3HitsLost ->resetBarcodes();
//                break;
//              default :
//                break;
//            }

          } else {

            switch ( pC.second ) // layer
            { 
              case 0 :
                trkHist_pseudo_reco_Lay0HitsFound -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay0HitsFound -> resetBarcodes();
                break;
              case 1 :
                trkHist_pseudo_reco_Lay1HitsFound -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay1HitsFound -> resetBarcodes();
                break;
              case 2:
                trkHist_pseudo_reco_Lay2HitsFound -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay2HitsFound -> resetBarcodes();
                break;
              case 3:
                trkHist_pseudo_reco_Lay3HitsFound -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay3HitsFound -> resetBarcodes();
                break;
              default :
                break;
            }
          }
        } // loop over hits on pseudotrack

        // loop over reco and look for match in pseudo collection
        for( auto rC : rIDs ) {
          if( pIDs.count( rC.first ) == 0 ) {
            // fill hist for spurious clusters
            switch ( rC.second ) // layer
            { 
              case 0 :
                trkHist_pseudo_reco_Lay0HitsExtra -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay0HitsExtra -> resetBarcodes();
                break;
              case 1 :
                trkHist_pseudo_reco_Lay1HitsExtra -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay1HitsExtra -> resetBarcodes();
                break;
              case 2:
                trkHist_pseudo_reco_Lay2HitsExtra -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay2HitsExtra -> resetBarcodes();
                break;
              case 3:
                trkHist_pseudo_reco_Lay3HitsExtra -> FillHists( (*pseudo_itr), 1.0 );
                trkHist_pseudo_reco_Lay3HitsExtra -> resetBarcodes();
                break;
              default :
                break;
            }
          }
        }
      } // if the match exists and so do MSOS


//      // comparing numbers of hits - this is too simplistic
//      if(theMatch) {
//        // for each layer if pseudotrack has a hit check if reco track has a hit on same layer
//        std::vector<int> pseudoHits  = (*pseudo_itr) ->auxdata< std::vector<int> >( "nHitsOnPixelLayer" );
//        std::vector<int> recoHits    = (theMatch)    ->auxdata< std::vector<int> >( "nHitsOnPixelLayer" );
//        if( pseudoHits.at(0) > 0 ) {
//          for( int i = 0; i < pseudoHits.at(0); i++ ) {
//            trkHist_pseudo_reco_Lay0Hits -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay0Hits ->resetBarcodes();
//          }
//          for( int i = 0; i < recoHits.at(0); i++ ) {
//            trkHist_pseudo_reco_Lay0HitsFound -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay0HitsFound ->resetBarcodes();
//          }
//        }
//        if( pseudoHits.at(1) > 0 ) {
//          for( int i = 0; i < pseudoHits.at(1); i++ ) {
//            trkHist_pseudo_reco_Lay1Hits -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay1Hits ->resetBarcodes();
//          }
//          for( int i = 0; i < recoHits.at(1); i++ ) {
//            trkHist_pseudo_reco_Lay1HitsFound -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay1HitsFound ->resetBarcodes();
//          }
//        }
//        if( pseudoHits.at(2) > 0 ) {
//          for( int i = 0; i < pseudoHits.at(2); i++ ) {
//            trkHist_pseudo_reco_Lay2Hits -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay2Hits ->resetBarcodes();
//          }
//          for( int i = 0; i < recoHits.at(2); i++ ) {
//            trkHist_pseudo_reco_Lay2HitsFound -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay2HitsFound ->resetBarcodes();
//          }
//        }
//        if( pseudoHits.at(3) > 0 ) {
//          for( int i = 0; i < pseudoHits.at(3); i++ ) {
//            trkHist_pseudo_reco_Lay3Hits -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay3Hits ->resetBarcodes();
//          }
//          for( int i = 0; i < recoHits.at(3); i++ ) {
//            trkHist_pseudo_reco_Lay3HitsFound -> FillHists( (*pseudo_itr), 1.0 );
//            trkHist_pseudo_reco_Lay3HitsFound ->resetBarcodes();
//          }
//        }
//      } // theMatch exists

    } // nMatchReco 0
  } // loop over pseudo tracks  
  
  if(!evtHasLost){ 
  eventHist_all_reco_pstrk  ->FillHists( eventFeatures, 1.0 );
  } else {
     eventHist_lost_pstrk -> FillHists( eventFeatures, 1.0 );
     if (evtHasLost0SSCT) {eventHist_lost0SSCT_pstrk   ->FillHists( eventFeatures, 1.0 );}
     else                 {eventHist_lostgt0SSCT_pstrk ->FillHists( eventFeatures, 1.0 );}
  }

  ///*  Get clusters on pseudotrack
  // Use the shallow copy of pseudoTrack container
  for(pseudo_itr = pseudoTracks->begin() ; pseudo_itr != pseudo_end; ++pseudo_itr ) {    
     
     // look only at lost pseudotracks with no shared SCT hits 
     
     uint8_t getInt(0);   // for accessing summary information
     if( (*pseudo_itr)->auxdata< int >( "nMatchReco" ) != 0 ) { continue; }
     (*pseudo_itr)->summaryValue( getInt, xAOD::numberOfSCTSharedHits ); 
     //if( getInt != 0 ) { continue; } 
    
     // ok, now fill histograms for clusters on track
     // some aliases
     typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
     typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
     static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change
     static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack(measurementNames);

     if(debugMarco) {std::cout << "\t lost pseudoTrack has barcode: " << (*pseudo_itr)->auxdata< int >("barcode") << " - nr. shared SCT hits: " << static_cast<int>(getInt) <<  " - Has pixel MSOSs? " << (*pseudo_itr)->isAvailable< MeasurementsOnTrack >( measurementNames ) << std::endl;}

     if( (*pseudo_itr)->isAvailable< MeasurementsOnTrack >( measurementNames ) ) 
     {  
       const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack( *(*pseudo_itr) );
       // now loop over this psudotrack TSOSs
       for( MeasurementsOnTrackIter msos_iter = measurementsOnTrack.begin(); msos_iter != measurementsOnTrack.end(); ++msos_iter )
       {
	 
	 //Check if the element link is valid
	 if( ! (*msos_iter).isValid() ) { continue; }	 
	 
	 const  xAOD::TrackStateValidation *msos = *(*msos_iter);
	 
	 if( !(msos->detType()==1) ) { continue; }
         // The measurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
         if( !(msos->type()==0) ) { continue; }  // confirm this is a measurement type

         // check if already plotting this MSOS's PRD's sibling
           // siblings are copies of the parent with a new local position
           // the identifier and positions can be different but the layer, size, and charge information will be the same
         bool skip(false);
         for( MeasurementsOnTrackIter  msos_newiter = measurementsOnTrack.begin(); msos_newiter != msos_iter; ++msos_newiter) {
	    const  xAOD::TrackStateValidation *this_msos = *(*msos_newiter);
            if ( xAOD::TrackHelper::areSiblings( msos, this_msos) ) { skip = true; break; }
         }
         if(skip) { continue; }
	 
         // fill cluster plots 
         fillPixelClusterHists( msos, pixClusters, pixClusHistsOnPseudoTrack_lost);
         if (getInt == 0){
           fillPixelClusterHists( msos, pixClusters, pixClusHistsOnPseudoTrack_lost0SSCT);
           if( (*pseudo_itr)->auxdata< int >( "nMatchSisp" ) != 0 ) {
             fillPixelClusterHists( msos, pixClusters, pixClusHistsOnPseudoTrack_lost0SSCT_YES_SiSpMatch);
           } else {
             fillPixelClusterHists( msos, pixClusters, pixClusHistsOnPseudoTrack_lost0SSCT_NO_SiSpMatch);
           }
         } 
       }
     }

  }
  //*/

  //---------------------
  // CLEAN UP
  //---------------------
 //delete truthPriParts;
 //delete truthPriPartsAux;
 delete pseudoTracks;
 delete pseudoTracksAux;
 delete recoTracksShallowCopyPair.first;
 delete recoTracksShallowCopyPair.second;
 delete sispTrksShallowCopyPair.first;
 delete sispTrksShallowCopyPair.second;
 delete pixClustersPair.first;
 delete pixClustersPair.second;
 delete sctClustersPair.first;
 delete sctClustersPair.second;

 return EL::StatusCode::SUCCESS;
}


EL::StatusCode TrkEffStudy :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrkEffStudy :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}

int TrkEffStudy :: buildPixelModuleID(const xAOD::TrackMeasurementValidation* prd) const
{
  int bec        = prd->auxdataConst< int >("bec");   // -4,4
  int layer      = prd->auxdataConst< int >("layer");      // 0,4
  int eta_module = prd->auxdataConst< int >("eta_module"); // -20,20
  int phi_module = prd->auxdataConst< int >("phi_module"); // 0-50

   return bec * 1000000 + layer * 100000 + phi_module * 100 + eta_module + 50; 

}


void TrkEffStudy :: addTruthDistancePerLayer( const xAOD::TrackParticle* tp, std::multimap<int, const xAOD::TrackMeasurementValidation*>& clustersOnAModule) const
{
  //Calculates the distance to the neaeast true paricle
  //Note the track needs to have a cluster on the layer in question 
  //If no other cluster is on that module then you will get the default value
  std::vector<float> resultsR(4,-9999);
  std::vector<float> resultsX(4,-9999);
  std::vector<float> resultsY(4,-9999);

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
  static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change

  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack( measurementNames );

  // Check if there are MSOS attached
  if( ! tp->isAvailable< MeasurementsOnTrack >( measurementNames ) ) {
    return;
  }
  // Get the MSOS's
  const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack( *tp );;

  // Loop over track TrackStateValidation's
  
  const xAOD::TruthParticle*  truth = xAOD::TrackHelper::truthParticle(tp);
  if(!truth)
  {
    //Warning("truthDistancePerLayer", "No Truth Partilce associated to the track -- return dummy vector");
    return;
  }
  
  for( MeasurementsOnTrackIter msos_iter = measurementsOnTrack.begin();
       msos_iter != measurementsOnTrack.end(); ++msos_iter)
  {  
    //Check if the element link is valid
    if( ! (*msos_iter).isValid() ) {
      continue;
    }
    
    const xAOD::TrackStateValidation* msos = *(*msos_iter); 
   
    // The detector type enum is defined in TrkEventPrimitives/TrackStateDefs.h
    // The mesurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
    // We will port the enums to xAOD world in the future
    if( msos->detType() != 1 || msos->type() ==  6    )  //its   pixel && and not a hole
      continue;
    
    //DEBUG("execute()",  "->MSOS is a pixel with the following ID %li", msos->detElementId()); 
    
    //Get pixel cluster
    if(  msos->trackMeasurementValidationLink().isValid() && *(msos->trackMeasurementValidationLink()) ){
      const xAOD::TrackMeasurementValidation* pixelCluster =  *(msos->trackMeasurementValidationLink());        
              
      // if it is not in the barrel
      if( pixelCluster->auxdata<int>("bec") != 0  )
        continue;
      
      int layer  = pixelCluster->auxdata<int>("layer");
      
      static SG::AuxElement::ConstAccessor< std::vector<int> >     acc_sihit_barcode("sihit_barcode");
      static SG::AuxElement::ConstAccessor< std::vector<float> >   acc_sihit_startPosX("sihit_startPosX");
      static SG::AuxElement::ConstAccessor< std::vector<float> >   acc_sihit_startPosY("sihit_startPosY");
      static SG::AuxElement::ConstAccessor< std::vector<float> >   acc_sihit_endPosX("sihit_endPosX");
      static SG::AuxElement::ConstAccessor< std::vector<float> >   acc_sihit_endPosY("sihit_endPosY");


      const std::vector<int>   &  sihit_barcode   = acc_sihit_barcode(*pixelCluster);
      const std::vector<float> &  sihit_startPosX = acc_sihit_startPosX(*pixelCluster);
      const std::vector<float> &  sihit_startPosY = acc_sihit_startPosY(*pixelCluster);
      const std::vector<float> &  sihit_endPosX   = acc_sihit_endPosX(*pixelCluster);
      const std::vector<float> &  sihit_endPosY   = acc_sihit_endPosY(*pixelCluster);


      int myTP_Location(-1); 
      for( unsigned int i(0); i < sihit_barcode.size(); ++i ){
        if( sihit_barcode[i] == truth->barcode() ){
          myTP_Location = i;
          break;
        }
      }
      
      //Cluster is incorrectly assigned to the track
      if(myTP_Location == -1){
        continue;
      }

      float myhitX = (sihit_endPosX[myTP_Location] + sihit_startPosX[myTP_Location]) *0.5;
      float myhitY = (sihit_endPosY[myTP_Location] + sihit_startPosY[myTP_Location]) *0.5;

      std::vector<int> otherParticles;  
      for( unsigned int i(0); i < sihit_barcode.size(); ++i ){
        if( (int)i == myTP_Location){
          continue;
        }
        //only consider primary particles
        if( sihit_barcode[i] >0  && sihit_barcode[i]< 200000 ){
          otherParticles.push_back(i);
        }
      }

      //Multiple paricles in the cluster which is the closest one 
      if (otherParticles.size() > 0 ){
        float minDeltaR2(100000);
        int closest(0);
        for( auto i : otherParticles ){
          float dx = myhitX - (sihit_endPosX[i] + sihit_startPosX[i]) *0.5;
          float dy = myhitY - (sihit_endPosY[i] + sihit_startPosY[i]) *0.5;
          float deltaR2 = dx*dx + dy*dy;
          if(deltaR2 < minDeltaR2){
            closest =  (int)i;
            minDeltaR2 = deltaR2;
          }
        }
        resultsR[layer] = minDeltaR2;
        resultsX[layer] = myhitX - (sihit_endPosX[closest] + sihit_startPosX[closest]) *0.5;
        resultsY[layer] = myhitY - (sihit_endPosY[closest] + sihit_startPosY[closest]) *0.5;
        continue;
      }


      //Get all of the clusters on that module
      int clusterModuleID = buildPixelModuleID( pixelCluster ); 
      // equal_range(**) returns pair<iterator,iterator> representing the range
      // of element with key **
      
      auto pair_it =  clustersOnAModule.equal_range(clusterModuleID);
    
      //std::cout << "Layer :" <<  layer << " Cluster does not contain more than one truth particle -- look for other clusters on the module" << std::endl;
      //std::cout << "Module ID: " << clusterModuleID << " other clusters " <<  clustersOnAModule.count(clusterModuleID) << "  "<< clustersOnAModule.size() << std::endl;

      //Find the closest cluster in the same module;
      const xAOD::TrackMeasurementValidation* closestCluster = 0;
      float minDeltaR2(1.e20);
      
      //Copy the barcode vector so we can use it a comparison later
      const std::vector<int>     sihit_barcode_Copy = sihit_barcode;
      // Loop through range of maps of key **
      for (auto it = pair_it.first;
          it != pair_it.second;
          ++it)
      {
        // Note you can not use the pointers directly as they may change 
        if(it->second->identifier() == pixelCluster->identifier() )
          continue;
        
        // Double check that if the identifier differ but if they contain exactly the same truth particles continue -- split cluster 
        const std::vector<int> &    sihit_barcode_Other   = acc_sihit_barcode(*it->second);
        if( sihit_barcode_Copy == sihit_barcode_Other )
          continue;
        
        float dx = pixelCluster->localX() - it->second->localX();  
        float dy = pixelCluster->localY() - it->second->localY();
        float deltaR2 = dx*dx + dy*dy;
        if(deltaR2 < minDeltaR2){
          closestCluster = it->second;
          minDeltaR2 = deltaR2;
        }    
      }
       
      // NO other cluster on the Module
      if (!closestCluster){
        //std::cout << "No cluster found on module" << std::endl;
        continue;
      } 

      const std::vector<int> &    sihit_barcode_OC   = acc_sihit_barcode(*closestCluster);
      const std::vector<float> &  sihit_startPosX_OC = acc_sihit_startPosX(*closestCluster);
      const std::vector<float> &  sihit_startPosY_OC = acc_sihit_startPosY(*closestCluster);
      const std::vector<float> &  sihit_endPosX_OC   = acc_sihit_endPosX(*closestCluster);
      const std::vector<float> &  sihit_endPosY_OC   = acc_sihit_endPosY(*closestCluster);
      
      otherParticles.clear();  
      for( unsigned int i(0); i < sihit_barcode_OC.size(); ++i ){
        if( sihit_barcode_OC[i] == truth->barcode() ){
          //std::cout << "Truth particle has same barcode as the one i am interested in" << std::endl;  

          //float dx = myhitX - (sihit_endPosX_OC[i] + sihit_startPosX_OC[i]) *0.5;
          //float dy = myhitY - (sihit_endPosY_OC[i] + sihit_startPosY_OC[i]) *0.5;
          //std::cout << "dx " << dx << " dy " <<  dy << std::endl; 

          continue;
        }
        //only consider primary particles
        if( sihit_barcode_OC[i] >0  && sihit_barcode_OC[i]< 200000 ){
          otherParticles.push_back(i);
        }
      }

      //Multiple paricles in the cluster which is the closeest one 
      if (otherParticles.size() > 0 ){
        float minDeltaR2(100000);
        int closest(0);
        for( auto i : otherParticles ){
          float dx = myhitX - (sihit_endPosX_OC[i] + sihit_startPosX_OC[i]) *0.5;
          float dy = myhitY - (sihit_endPosY_OC[i] + sihit_startPosY_OC[i]) *0.5;
          float deltaR2 = dx*dx + dy*dy;
          if(deltaR2 < minDeltaR2){
            closest =  (int)i;
            minDeltaR2 = deltaR2;
          }
        }
        //std::cout << "Found another truth particle in another cluster" << std::endl;  
        
        resultsR[layer] = minDeltaR2;
        resultsX[layer] = myhitX - (sihit_endPosX_OC[closest] + sihit_startPosX_OC[closest]) *0.5;
        resultsY[layer] = myhitY - (sihit_endPosY_OC[closest] + sihit_startPosY_OC[closest]) *0.5;
        continue;
      } else {
        //std::cout << "Everything look like it comes from the same truth particle" << std::endl;  
      }
    } // End Check on Element link     
  } // End for loop over MSOS's
  
  //dress the TrackParticle
 
  static SG::AuxElement::Decorator< float > truedist1pixlay("truedist1pixlay");
  static SG::AuxElement::Decorator< float > truedist2pixlay("truedist2pixlay");
  static SG::AuxElement::Decorator< float > truedist3pixlay("truedist3pixlay");
  static SG::AuxElement::Decorator< float > truedist4pixlay("truedist4pixlay");

  static SG::AuxElement::Decorator< float > truedistEta1pixlay("truedistEta1pixlay");
  static SG::AuxElement::Decorator< float > truedistEta2pixlay("truedistEta2pixlay");
  static SG::AuxElement::Decorator< float > truedistEta3pixlay("truedistEta3pixlay");
  static SG::AuxElement::Decorator< float > truedistEta4pixlay("truedistEta4pixlay");

  static SG::AuxElement::Decorator< float > truedistPhi1pixlay("truedistPhi1pixlay");
  static SG::AuxElement::Decorator< float > truedistPhi2pixlay("truedistPhi2pixlay");
  static SG::AuxElement::Decorator< float > truedistPhi3pixlay("truedistPhi3pixlay");
  static SG::AuxElement::Decorator< float > truedistPhi4pixlay("truedistPhi4pixlay");
 
  truedist1pixlay(*tp) = sqrt(resultsR[0]);
  truedist2pixlay(*tp) = sqrt(resultsR[1]);
  truedist3pixlay(*tp) = sqrt(resultsR[2]);
  truedist4pixlay(*tp) = sqrt(resultsR[3]);

  truedistEta1pixlay(*tp) = resultsY[0];
  truedistEta2pixlay(*tp) = resultsY[1];
  truedistEta3pixlay(*tp) = resultsY[2];
  truedistEta4pixlay(*tp) = resultsY[3];

  truedistPhi1pixlay(*tp) = resultsX[0];
  truedistPhi2pixlay(*tp) = resultsX[1];
  truedistPhi3pixlay(*tp) = resultsX[2];
  truedistPhi4pixlay(*tp) = resultsX[3];

  
}

void TrkEffStudy :: filldetIdtoMsosMap( const xAOD::TrackParticleContainer* tc, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const
{
	for(auto trk : (*tc))
	{
		auto measureOnTrack = xAOD::TrackHelper::measureOnTrack(trk);
		for( auto msos_iter : measureOnTrack)
		{
			//Check if the element link is valid
			if( !(msos_iter).isValid() )	continue;
		
			auto msos = *(msos_iter);
			
			if( (msos->detType() != 1 &&	msos->detType() != 2) ||  msos->type() != 0 )  //its a pixel or SCT && and a measurement on track
				continue;
		
			detIdtoMsosMap.insert(std::pair<uint64_t,const xAOD::TrackStateValidation*>(msos->detElementId(),msos));
		} // End for loop over MSOS's
	}
}

void TrkEffStudy :: MinDistanceOtherTrack( const xAOD::TrackParticle* tp, const xAOD::TrackParticleContainer* tc, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const
{

  float mindistPix[4] = {1e5,1e5,1e5,1e5};
	float mindistSCT[9] = {1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5};

  int layer = -1;
	auto measureOnTrack = xAOD::TrackHelper::measureOnTrack(tp);
	for (auto msos_iter : measureOnTrack) {
		//Check if the element link is valid
		if (!(msos_iter).isValid())	{ continue; }

		auto thisMsos = *(msos_iter);
		if ((thisMsos->detType()!= 1 &&	thisMsos->detType()!= 2) || thisMsos->type()!=0) { continue; }  //its a pixel or SCT && and a measurement on track
		
		if (thisMsos->trackMeasurementValidationLink().isValid() && *(thisMsos->trackMeasurementValidationLink())) {
			auto cluster =  *(thisMsos->trackMeasurementValidationLink());       
			layer =	cluster->auxdataConst<int>("layer");
		}
		auto range = detIdtoMsosMap.equal_range(thisMsos->detElementId());
		for (auto otherMsos=range.first; otherMsos!=range.second; ++otherMsos) {
			if (otherMsos->second==thisMsos) continue;
			float loc_dx = otherMsos->second->localX()-thisMsos->localX();
			float loc_dy = otherMsos->second->localY()-thisMsos->localY();
			float dist = sqrt(loc_dx*loc_dx+loc_dy*loc_dy); 
			if (thisMsos->detType()==1) { 
        if (dist<mindistPix[layer]) { mindistPix[layer] = dist; }
      }
      else { 
        if (dist<mindistSCT[layer]) { mindistSCT[layer] = dist; }
      }
		}
	} // End for loop over MSOS's

  static SG::AuxElement::Decorator< float > dist1pixlay("dist1pixlay");
  static SG::AuxElement::Decorator< float > dist2pixlay("dist2pixlay");
  static SG::AuxElement::Decorator< float > dist3pixlay("dist3pixlay");
  static SG::AuxElement::Decorator< float > dist4pixlay("dist4pixlay");
  dist1pixlay(*tp) = mindistPix[0];
  dist2pixlay(*tp) = mindistPix[1];
  dist3pixlay(*tp) = mindistPix[2];
  dist4pixlay(*tp) = mindistPix[3];

	static SG::AuxElement::Decorator< float > dist1sctlay("dist1sctlay");
	static SG::AuxElement::Decorator< float > dist2sctlay("dist2sctlay");
	static SG::AuxElement::Decorator< float > dist3sctlay("dist3sctlay");
	static SG::AuxElement::Decorator< float > dist4sctlay("dist4sctlay");
	static SG::AuxElement::Decorator< float > dist5sctlay("dist5sctlay");
	static SG::AuxElement::Decorator< float > dist6sctlay("dist6sctlay");
	static SG::AuxElement::Decorator< float > dist7sctlay("dist7sctlay");
	static SG::AuxElement::Decorator< float > dist8sctlay("dist8sctlay");
	static SG::AuxElement::Decorator< float > dist9sctlay("dist9sctlay");
  dist1sctlay(*tp) = mindistSCT[0];
  dist2sctlay(*tp) = mindistSCT[1];
  dist3sctlay(*tp) = mindistSCT[2];
  dist4sctlay(*tp) = mindistSCT[3];
  dist5sctlay(*tp) = mindistSCT[4];
  dist6sctlay(*tp) = mindistSCT[5];
  dist7sctlay(*tp) = mindistSCT[6];
  dist8sctlay(*tp) = mindistSCT[7];
  dist9sctlay(*tp) = mindistSCT[8];


/*
  std::vector<float> mindistPix(4,std::numeric_limits<float>::max());
	std::vector<float> mindistSCT(4,std::numeric_limits<float>::max());
	std::vector<float>* mindist = 0;
	int layer = -1;
	auto measureOnTrack = xAOD::TrackHelper::measureOnTrack(tp);
	for( auto msos_iter : measureOnTrack )
	{
		//Check if the element link is valid
		if( !(msos_iter).isValid() )	continue;
	
		auto thisMsos = *(msos_iter);
		
		if( (thisMsos->detType() != 1 &&	thisMsos->detType() != 2) ||  thisMsos->type() != 0 )  //its a pixel or SCT && and a measurement on track
			continue;
		
		if(  thisMsos->trackMeasurementValidationLink().isValid() && *(thisMsos->trackMeasurementValidationLink()) )
		{
			auto cluster =  *(thisMsos->trackMeasurementValidationLink());       
			layer =	cluster->auxdataConst<int>("layer");
			if( thisMsos->detType()==1) mindist = &mindistPix;
			else mindist = &mindistSCT;
		}
		
		auto range = detIdtoMsosMap.equal_range(thisMsos->detElementId());
		for( auto otherMsos=range.first; otherMsos!=range.second; ++otherMsos ){
			if(otherMsos->second==thisMsos) continue;
			float loc_dx =	otherMsos->second->localX()	-	thisMsos->localX();
			float loc_dy =	otherMsos->second->localY()	-	thisMsos->localY();
			float dist = sqrt( loc_dx*loc_dx + loc_dy*loc_dy ); 
			if( dist < (*mindist)[layer] )	(*mindist)[layer] = dist;
		}
	} // End for loop over MSOS's
	
	static SG::AuxElement::Decorator< std::vector<float> > distpixlayer("distpixlayer");
	for( auto minDistLayer : mindistPix )		distpixlayer(*tp).push_back(minDistLayer);
	static SG::AuxElement::Decorator< std::vector<float> > distsctlayer("distsctlayer");
	for( auto minDistLayer : mindistSCT )		distsctlayer(*tp).push_back(minDistLayer);
*/


}

void TrkEffStudy :: fillPixelClusterHists ( const xAOD::TrackStateValidation* msos, xAOD::TrackMeasurementValidationContainer* pixClusters, std::vector<ClusterHistManager*> &clusHists) {

  // find the PRD shallow copy so can add number of properties i.e. number of contributing particles
  // Pixel 
  xAOD::TrackMeasurementValidationContainer::iterator pixclus_itr = pixClusters->begin();
  xAOD::TrackMeasurementValidationContainer::iterator pixclus_end = pixClusters->end();
  xAOD::TrackMeasurementValidation* pixprd = 0;
  for( ; pixclus_itr != pixclus_end; ++pixclus_itr ) {
    if( (*pixclus_itr)->identifier() != (*(msos->trackMeasurementValidationLink()))->identifier() ) { continue; }
    if( (*pixclus_itr)->auxdata< float >( "charge" ) != (*(msos->trackMeasurementValidationLink()))->auxdata< float >( "charge" ) ) { continue; }
    pixprd = (*pixclus_itr);
    break;
  }
  if(!pixprd) { 
    std::cout << "Cannot find MSOS PRD in Pixel Cluster container" << std::endl;
    return; 
  }

  //if( pixprd->auxdata< int >( "nPrimaryParticles" ) == 0 ) { return; }
  if(!s_doClusterHists) { return; }

  for (auto histIt = clusHists.begin(); histIt != clusHists.end(); ++histIt){
    (*histIt)->FillHists( pixprd, 1.0 );
    // fill different components of histogram vector based on prd detType
  }

} // fillPixelClusterHists

void TrkEffStudy :: countContributingParticles( xAOD::TrackMeasurementValidation* prd ) {

    std::vector<int> bcodes;
    int nPrimaryParticles(0);
    const std::vector< std::vector<int> > sdo_barcodes = prd->auxdata< std::vector< std::vector<int> > >( "sdo_depositsBarcode" );
    for(unsigned int i = 0; i < sdo_barcodes.size(); i++) { 
      for(unsigned int j = 0; j < sdo_barcodes.at(i).size(); j++) { 
        int bc( sdo_barcodes.at(i).at(j) );
        std::vector<int>::iterator it = find(bcodes.begin(), bcodes.end(), bc);
        if( it == bcodes.end() ) { // if at end then new barcode
          bcodes.push_back( bc );
          if( bc <= 10e3 || bc >= 200e3 ) { continue; }
          nPrimaryParticles++;
        }
      }
    } // loop over sdo barcodes
    prd->auxdata< int >( "nParticles" ) = bcodes.size();
    prd->auxdata< int >( "nPrimaryParticles" ) = nPrimaryParticles;
    prd->auxdata< std::vector<int> >( "SDOBarcodes") = bcodes;

}


void TrkEffStudy :: dressLayerHitContent( xAOD::TrackParticleContainer* trks, xAOD::TrackMeasurementValidationContainer* pixClusters, xAOD::TrackMeasurementValidationContainer* sctClusters) {

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
  static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change

  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack(measurementNames);

  // for each track, fill in vector for how many particles are contributing to each hit 
  xAOD::TrackParticleContainer::iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::iterator trk_end = trks->end();
  for( ; trk_itr != trk_end; ++trk_itr ) {
    std::vector<int> nPartPix    = { -1, -1, -1, -1 };
    std::vector<int> nPriPartPix = { -1, -1, -1, -1 };
    ///*
    std::vector<int> nPartSCT    = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> nPriPartSCT = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    //*/
    /*
    std::vector<int> nPartSCT	 = { -1, -1, -1, -1 };
    std::vector<int> nPriPartSCT = { -1, -1, -1, -1 };
    */
    std::vector<int> nHitsOnPixelLayer = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<int> nHitsOnSCTLayer   = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };	

    (*trk_itr)->auxdata< float >( "BLayerCharge"  )  = -1;
    (*trk_itr)->auxdata< int >( "BLayerHitSize" )    = -1;
    (*trk_itr)->auxdata< int >( "BLayerHitSizeZ" )   = -1;
    (*trk_itr)->auxdata< int >( "BLayerHitSizePhi" ) = -1;
    (*trk_itr)->auxdata< int >( "nBLayerMergedHits" ) = 0;

//std::cout << "dressing 0" << std::endl;

    // Get the MSOS's
    if( ! (*trk_itr)->isAvailable< MeasurementsOnTrack >( measurementNames ) ) { 
      (*trk_itr)->auxdata< std::vector<int> >( "pixClusPartContent" )    = nPartPix;
      (*trk_itr)->auxdata< std::vector<int> >( "pixClusPriPartContent" ) = nPriPartPix;
      (*trk_itr)->auxdata< std::vector<int> >( "sctClusPartContent" )    = nPartSCT;
      (*trk_itr)->auxdata< std::vector<int> >( "sctClusPriPartContent" ) = nPriPartSCT;
      
      (*trk_itr)->auxdata< std::vector<int> >( "nHitsOnPixelLayer" ) = nHitsOnPixelLayer;
      (*trk_itr)->auxdata< std::vector<int> >( "nHitsOnSCTLayer" )   = nHitsOnSCTLayer;
      continue; 
    }

//std::cout << "dressing 1" << std::endl;

    const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack( *(*trk_itr) );
    // Loop over track TrackStateValidation's
    for( MeasurementsOnTrackIter msos_iter = measurementsOnTrack.begin();
        msos_iter != measurementsOnTrack.end(); ++msos_iter)
    {  

//std::cout << "\t dressing 1-0" << std::endl;

      //Check if the element link is valid
      if( ! (*msos_iter).isValid() ) {
        continue;
      }

      const xAOD::TrackStateValidation* msos = *(*msos_iter); 
      if(  !msos->trackMeasurementValidationLink().isValid() ) { continue; }
      if(  !(*(msos->trackMeasurementValidationLink())) )      { continue; }

//std::cout << "\t dressing 1-1" << std::endl;


      // The detector type enum is defined in TrkEventPrimitives/TrackStateDefs.h
      // The mesurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
      // We will port the enums to xAOD world in the future
      if( msos->detType() == 1 ) { // its a pixel 
        if( msos->type() ==  6 ) { continue; } // not a hole

        //Get pixel cluster then match to clusters in container
        const xAOD::TrackMeasurementValidation* msosClus =  *(msos->trackMeasurementValidationLink());        

//std::cout << "\t is pixel - dressing 1-2" << std::endl;

        // find this cluster in the shallow copy container
        xAOD::TrackMeasurementValidationContainer::iterator clus_itr = (pixClusters)->begin();
        xAOD::TrackMeasurementValidationContainer::iterator clus_end = (pixClusters)->end();
        for( ; clus_itr != clus_end; ++clus_itr ) {
          if( (*clus_itr)->identifier() != (msosClus)->identifier() ) { continue; }
          if( (*clus_itr)->auxdata< float >( "charge" ) != (msosClus)->auxdata< float >( "charge" ) ) { continue; }

//std::cout << "\t\t  is pixel - dressing 1-1-0" << std::endl;

 
          countContributingParticles((*clus_itr));

//std::cout << "\t\t  is pixel - dressing 1-1-1" << std::endl;

          // multiple hits on one layer? FIXME - this should fix it
	  if( nHitsOnPixelLayer.at((*clus_itr)->auxdata< int >( "layer" )) < 2 ){
            nPartPix.at( (*clus_itr)->auxdata< int >( "layer" ) ) = (*clus_itr)->auxdata< int >( "nParticles" );
            nPriPartPix.at( (*clus_itr)->auxdata< int >( "layer" ) ) = (*clus_itr)->auxdata< int >( "nPrimaryParticles" );
          }
	  // for the track under exam, count the number of hits on each layer
          switch( (*clus_itr)->auxdata< int >( "layer" ) )
          {
            case 0 :
              nHitsOnPixelLayer[0]++;
              break;
            case 1 :
              nHitsOnPixelLayer[1]++;
              break;
            case 2 :
              nHitsOnPixelLayer[2]++;
              break;            
            case 3 :
              nHitsOnPixelLayer[3]++;
              break;	    
            default:
              std::cout << "Could not find pixel layer" << std::endl;
              break;
          }

          if( (*clus_itr)->auxdata< int >( "layer" ) == 0 && (*trk_itr)->auxdata< float >( "BLayerCharge"  ) == -1 ) {
            (*trk_itr)->auxdata< float >( "BLayerCharge"  )  = (*clus_itr)->auxdata< float >( "charge" );
            (*trk_itr)->auxdata< int >( "BLayerHitSize" )    = (*clus_itr)->auxdata< int >( "size" );
            (*trk_itr)->auxdata< int >( "BLayerHitSizeZ" )   = (*clus_itr)->auxdata< int >( "sizeZ" );
            (*trk_itr)->auxdata< int >( "BLayerHitSizePhi" ) = (*clus_itr)->auxdata< int >( "sizePhi" );
            if( (*clus_itr)->auxdata< int >( "nPrimaryParticles" ) > 1 ) {
              (*trk_itr)->auxdata< int >( "nBLayerMergedHits" ) = (*trk_itr)->auxdata< int >( "nBLayerMergedHits" ) + 1;
            }
          }
          
          break; // once is enough! Move to next MSOS...
        } // loop over cluster container looking for a match
      } // pixels
      else if( msos->detType() == 2 ) { // SCT
        if( msos->type() ==  6 ) { continue; } // not a hole
      
        //Get sct cluster then match to clusters in container
        const xAOD::TrackMeasurementValidation* msosClus =  *(msos->trackMeasurementValidationLink());
	
//std::cout << "\t is sct - dressing 1-2" << std::endl; 
	        
        // find this cluster in the shallow copy container
        xAOD::TrackMeasurementValidationContainer::iterator sct_itr = (sctClusters)->begin();
        xAOD::TrackMeasurementValidationContainer::iterator sct_end = (sctClusters)->end();
        for( ; sct_itr != sct_end; ++sct_itr ) {
          if( (*sct_itr)->identifier() != (msosClus)->identifier() ) { continue; }
          //if( (*sct_itr)->auxdata< float >( "charge" ) != (msosClus)->auxdata< float >( "charge" ) ) { continue; }

          // what is layer number of 1st SCT layer ? 3 for Run I 4 for Run II...
          int sct_layer = (*sct_itr)->auxdata< int >( "layer" )/* - 3*/;

// ****  how is sct_layer indexed?

//std::cout << " sct_layer idx: " << sct_layer << std::endl;


//std::cout << "\t\t  is sct - dressing 1-1-0" << std::endl;

          countContributingParticles((*sct_itr));

//std::cout << "\t\t  is sct - dressing 1-1-1" << std::endl;

	  // consider only first measurement on layer
	  if( nHitsOnSCTLayer.at( sct_layer ) < 2 ){	  
            nPartSCT.at( sct_layer )    = (*sct_itr)->auxdata< int >( "nParticles" );
            nPriPartSCT.at( sct_layer ) = (*sct_itr)->auxdata< int >( "nPrimaryParticles" );
          }
//std::cout << "\t\t  is sct - dressing 1-1-2" << std::endl;	
  
	  // for the track under exam, count the number of hits on each layer - NB: 4 (double) SCT layers!
          switch( (*sct_itr)->auxdata< int >( "layer" ) )
          {
            case 0 :
              nHitsOnSCTLayer[0]++;
              break;
            case 1 :
              nHitsOnSCTLayer[1]++;
              break;
            case 2 :
              nHitsOnSCTLayer[2]++;
              break;            
            case 3 :
              nHitsOnSCTLayer[3]++;
              break;	
            case 4 :
              nHitsOnSCTLayer[4]++;
              break;
            case 5 :
              nHitsOnSCTLayer[5]++;
              break;
            case 6 :
              nHitsOnSCTLayer[6]++;
              break;            
            case 7 :
              nHitsOnSCTLayer[7]++;
              break;	
            case 8 :
              nHitsOnSCTLayer[8]++;
              break;			
            default:
              std::cout << "Could not find SCT layer" << std::endl;
              break;
          }

          // would be good to do for 1st SCT layer
	  //if( (*clus_itr)->auxdata< int >( "layer" ) == 0 && (*trk_itr)->auxdata< float >( "BLayerCharge"  ) == -1 ) {
          //  (*trk_itr)->auxdata< float >( "BLayerCharge"  )  = (*clus_itr)->auxdata< float >( "charge" );
          //  (*trk_itr)->auxdata< int >( "BLayerHitSize" )    = (*clus_itr)->auxdata< int >( "size" );
          //  (*trk_itr)->auxdata< int >( "BLayerHitSizeZ" )   = (*clus_itr)->auxdata< int >( "sizePhi" );
          //  (*trk_itr)->auxdata< int >( "BLayerHitSizePhi" ) = (*clus_itr)->auxdata< int >( "sizeZ" );
          //  if( (*clus_itr)->auxdata< int >( "nPrimaryParticles" ) > 1 ) {
          //    (*trk_itr)->auxdata< int >( "nBLayerMergedHits" ) = (*trk_itr)->auxdata< int >( "nBLayerMergedHits" ) + 1;
          //  }
          //}	
	  
          break; // once is enough! Move to next MSOS...
        } // sct
      } // pixel or sct msos  
    } // track's MSOS

//std::cout << "dressing 2" << std::endl;
    
    (*trk_itr)->auxdata< std::vector<int> >( "pixClusPartContent" )    = nPartPix;
    (*trk_itr)->auxdata< std::vector<int> >( "pixClusPriPartContent" ) = nPriPartPix;
    (*trk_itr)->auxdata< std::vector<int> >( "sctClusPartContent" )    = nPartSCT;
    (*trk_itr)->auxdata< std::vector<int> >( "sctClusPriPartContent" ) = nPriPartSCT;

    (*trk_itr)->auxdata< std::vector<int> >( "nHitsOnPixelLayer" ) = nHitsOnPixelLayer;
    (*trk_itr)->auxdata< std::vector<int> >( "nHitsOnSCTLayer" )   = nHitsOnSCTLayer;

  } // loop over tracks
} // dressLayerHitContent

const xAOD::TrackParticle* TrkEffStudy :: getHighestProbTrack( std::multimap<int, const xAOD::TrackParticle* >& bcodeMap, int bcode ) {

  typedef std::multimap<int, const xAOD::TrackParticle* >::iterator barcodeMapIter;
  std::pair <barcodeMapIter, barcodeMapIter> matches = bcodeMap.equal_range( bcode );

  float bestScore = 0;
  const xAOD::TrackParticle*  bestMatch = 0;
  for( barcodeMapIter it = matches.first; it!=matches.second; ++it ) { 
    if( (*it).second->isAvailable< float >(  "truthMatchProbability" ) &&  
        (*it).second->auxdata< float >( "truthMatchProbability" ) > bestScore ) { 
      bestScore = (*it).second->auxdata< float >( "truthMatchProbability" );
      bestMatch = (*it).second ;
    }
  }
  return bestMatch;

} // getHighestProbTrack

std::vector< const xAOD::TrackMeasurementValidation* > TrkEffStudy :: getAllPixelClusters( const std::vector<ElementLink< xAOD::TrackStateValidationContainer > >& mot ) {
  std::vector< const xAOD::TrackMeasurementValidation* > pixClus;
  for( std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator msos_iter = mot.begin(); msos_iter != mot.end(); ++msos_iter)
  {
    if( ! (*msos_iter).isValid() ) { continue; }
    const  xAOD::TrackStateValidation *msos = *(*msos_iter);
    if(  !msos->trackMeasurementValidationLink().isValid() ) { continue; }
    if(  !(*(msos->trackMeasurementValidationLink())) )      { continue; }
    if( msos->detType() != 1 ) { continue; } // not a pixel
    if( msos->type() ==  6 )   { continue; } // is a hole
    //Get pixel cluster 
    const xAOD::TrackMeasurementValidation* msosClus =  *(msos->trackMeasurementValidationLink());        
    if( msosClus ) { pixClus.push_back( msosClus ); }
  }
  return pixClus;
}


