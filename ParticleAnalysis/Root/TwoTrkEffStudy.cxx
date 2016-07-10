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
#include <ParticleAnalysis/TwoTrkEffStudy.h>
#include <HistManager/EventFeaturesInterface.h>
#include <Utility/TrackHelper.h>

// EDM includes: - if move to header file will not compile!
#include "AthContainers/ConstDataVector.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTracking/TrackParticleAuxContainer.h"

#include "xAODTracking/VertexContainer.h"

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "TSystem.h"
#include "TFile.h"

#include <map>
#include <limits>
#include <string>
#include <sstream>

// this is needed to distribute the algorithm to the workers
ClassImp(TwoTrkEffStudy)

// switches - made static because need to be set (in job steering macro) before calling constructor
char* TwoTrkEffStudy::pileupWeightFile = "";

TwoTrkEffStudy::TwoTrkEffStudy() {

  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}


EL::StatusCode TwoTrkEffStudy::setupJob(EL::Job& job) {
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

EL::StatusCode TwoTrkEffStudy::histInitialize() {
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  // In TrackHistManager constructor:
  // 1st bool is to fill plots for tracks matched to primary particles
  // 2nd bool is to fill plots for tracks matched to secondary particles
  // 3rd bool is to fill plots for fake tracks 
  // 4th bool is if to fill duplicates (need 1st and 2nd be true)
  // 5th bool is if to fill plots based on hit content on track (need 1st be true)
  
  // NN-cluster splitting study using data/mc
  twotrkHist_reco_pair         = new TwoTrackHists("reco_pair");
  twotrkHist_minX              = new TwoTrackHists("minX");
  twotrkHist_minY              = new TwoTrackHists("minY");
  trkHist_reco_pair            = new TrackHists("reco_pair");
  trkHist_reco_pair_leading    = new TrackHists("reco_pair_leading");
  trkHist_reco_pair_subleading = new TrackHists("reco_pair_subleading");

  clusHist_reco_pair_leading    = new ClusterHistManager("reco_pair_leading", 0, true);
  clusHist_reco_pair_subleading = new ClusterHistManager("reco_pair_subleading", 0, true);

//  clusHist_reco_pair_leading    = new ClusterHistManager("reco_pair_leading", 1, true);
//  clusHist_reco_pair_subleading = new ClusterHistManager("reco_pair_subleading", 1, true);

  clusHist_leading_low  = new ClusterHists("leading_low");
  clusHist_leading_high = new ClusterHists("leading_high");

  // event histograms
  eventHist_all               = new EventHists("all");

  twotrkHist_reco_pair         -> Init( wk() );
  twotrkHist_minX              -> Init( wk() );
  twotrkHist_minY              -> Init( wk() );
  trkHist_reco_pair            -> Init( wk() );
  trkHist_reco_pair_leading    -> Init( wk() );
  trkHist_reco_pair_subleading -> Init( wk() );

  clusHist_reco_pair_leading    -> Init( wk() );
  clusHist_reco_pair_subleading -> Init( wk() );

  clusHist_leading_low  -> Init( wk() );
  clusHist_leading_high -> Init( wk() );

  eventHist_all               ->Init( wk() );

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TwoTrkEffStudy::fileExecute() {
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TwoTrkEffStudy::changeInput (bool /* firstFile */) {
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  TFile* infile = wk()->inputFile ();
  std::cout << infile->GetName() << std::endl;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TwoTrkEffStudy::initialize() {
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

  // GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  const char* fullGRLFilePath = gSystem->ExpandPathName("data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml");
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(fullGRLFilePath);
  m_grl->setProperty( "GoodRunsListVec", vecStringGRL);
  m_grl->setProperty("PassThrough", false);   // if true (default) will ignore result of GRL and will just pass all events
  m_grl->initialize();

  // Trigger decision tool
  m_configTool = new TrigConf::xAODConfigTool("xAODConfigTool");
  ToolHandle<TrigConf::ITrigConfigTool> configHandle(m_configTool);
  if (configHandle->initialize().isFailure()) {
    Error("initialize()","Failed to initialize trigger config handle");
    return EL::StatusCode::FAILURE;
  }

  m_trigDecTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  m_trigDecTool->setProperty("ConfigTool"     , configHandle  );
  m_trigDecTool->setProperty("TrigDecisionKey","xTrigDecision");

  if (m_trigDecTool->initialize().isFailure()) {
    Error("initialize()","Failed to initialize trigger decision tool");
    return EL::StatusCode::FAILURE;
  }

  // pileup re-weighting tool
  m_tool = new CP::PileupReweightingTool("prw");
  std::vector<std::string> confFiles;
  std::vector<std::string> lcalcFiles;
//  confFiles.push_back("JZ3W.pileup_weight.root");
  confFiles.push_back(pileupWeightFile);
//  lcalcFiles.push_back("ilumicalc_histograms_None_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j100_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j110_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j150_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j175_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j200_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j260_267073-271744.root");
  lcalcFiles.push_back("ilumicalc_histograms_HLT_j300_267073-271744.root");

  dynamic_cast<CP::PileupReweightingTool&>(*m_tool).setProperty( "ConfigFiles", confFiles);
  dynamic_cast<CP::PileupReweightingTool&>(*m_tool).setProperty( "LumiCalcFiles", lcalcFiles);
  m_tool->initialize();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TwoTrkEffStudy::execute() {
  // Here you do everything that needs to be done on every single
  // event, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // print every 100 events, so we know where we are:
  if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
  //std::cout << "Event count " << m_eventCounter << std::endl;
  m_eventCounter++;
//  if(m_eventCounter > 1e5) { 
//    if(m_eventCounter == 1e5) { 
//      std::cout << "STOPPING after event " << m_eventCounter << std::endl;
//      std::cout << "Look in here to change -> TwoTrkEffStudy :: execute()" << std::endl;
//    }
//    return EL::StatusCode::SUCCESS; 
//  }

  EventFeatures eventFeatures;

  //---------------------------
  // Reset barcode map in histogram class
  //---------------------------
  trkHist_reco_pair            -> resetBarcodes();
  trkHist_reco_pair_leading    -> resetBarcodes();
  trkHist_reco_pair_subleading -> resetBarcodes();

  //---------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  if (!m_event->retrieve(eventInfo,"EventInfo").isSuccess()) {
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  xAOD::TStore store; // needed for shallow copies

  bool debugInitial(false);
  bool debug(false);
  if(debug) { Info("execute()", "Start Event"); }

  // check if the event is data or MC
  // (many tools are applied either to data or MC)
  bool isMC = false;
  // check if the event is MC
  if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
    isMC = true; // can do something with this later
  } 

  //---------------------
  // TRUTH PARTICLES
  //---------------------
  // get track container of interest
  const xAOD::TruthParticleContainer* truthPartsOrig = 0;
  if (isMC) {
    if ( !m_event->retrieve( truthPartsOrig, "TruthParticles" ).isSuccess() ){ // retrieve arguments: container type, container key
      Error("execute()", "Failed to retrieve Truth Particle container. Exiting." );
      return EL::StatusCode::FAILURE;
    }
  }

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
  if (debugInitial) {
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

  //===============
  // Jet container
  //===============
  const xAOD::JetContainer* jetAntiKt4EMTopo = 0;
  if (!m_event->retrieve(jetAntiKt4EMTopo,"AntiKt4EMTopoJets").isSuccess()) {
    Error("execute()", "Failed to retrieve AntiKt4EMTopoJets container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  //==================
  // Vertex container
  //==================
  const xAOD::VertexContainer* primaryVertex = 0;
  if (!m_event->retrieve(primaryVertex,"PrimaryVertices").isSuccess()) {
    Error("execute()", "Failed to retrieve PrimaryVertices container. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  //---------------------
  // PUT IN WORK
  //---------------------
  // if data check if event passes GRL
  if (!isMC) {
    if (!m_grl->passRunLB(*eventInfo)){ return EL::StatusCode::SUCCESS; }
  }

  // event cleaning
  if (!isMC) {
    if (eventInfo->errorState(xAOD::EventInfo::EventFlagSubDet::Tile)==xAOD::EventInfo::Error) { return EL::StatusCode::SUCCESS; }
    if (eventInfo->errorState(xAOD::EventInfo::EventFlagSubDet::LAr)==xAOD::EventInfo::Error)  { return EL::StatusCode::SUCCESS; }
    if (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18))                               { return EL::StatusCode::SUCCESS; }
  }

  // solve jet bug fix. 
  bool flag_jet850 = false;
  for (xAOD::JetContainer::const_iterator jet_itr=jetAntiKt4EMTopo->begin(); jet_itr!=jetAntiKt4EMTopo->end(); jet_itr++) {
//    std::cout << "jet Pt " << (*jet_itr)->pt()/1000.0 << std::endl;
    if ((*jet_itr)->pt()/1000.0>850.0) { flag_jet850 = true; }
  }
  if (flag_jet850) { return EL::StatusCode::SUCCESS; }

  // store primary vertex z-position
  float vtxZ = 0.0;
  for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
    if ((*vtx_itr)->vertexType()==1) {  // primary vertex
      vtxZ = (*vtx_itr)->z();
    }
  }


  // pileup weight
  m_tool->apply(*eventInfo);
  float pileupWeight = 1.0;
  float mcEvWeight = 1.0;
  float xsecWeight = 1.0;
  float muValue = 0.0;
  float muValueCorr = 0.0;

  if (isMC) {
    pileupWeight = eventInfo->auxdata<float>("PileupWeight");
    muValue     = eventInfo->auxdata<float>("averageInteractionsPerCrossing");
    muValueCorr = eventInfo->auxdata<float>("corrected_averageInteractionsPerCrossing");
    mcEvWeight = eventInfo->mcEventWeight(0);

    int mcDSID = eventInfo->mcChannelNumber();
    if      (mcDSID==361020) { xsecWeight=78.420e6*1.0/400000.0; }         // JZ0W : 361020 : 400000 : 78.420mb : 1.00
    else if (mcDSID==361021) { xsecWeight=78.420e6*6.7198e-4/399500.0; }   // JZ1W : 361021 : 400000 : 78.420mb : 6.7198E-04
    else if (mcDSID==361022) { xsecWeight= 2.433e6*3.3264e-4/286000.0; }   // JZ2W : 361022 : 400000 :  2.433mb : 3.3264E-04  2433200 2.4332e6
    else if (mcDSID==361023) { xsecWeight=26.454e3*3.1953e-4/400000.0; }   // JZ3W : 361023 : 400000 : 26.454ub : 3.1953E-04    26454 26.454e3
    else if (mcDSID==361024) { xsecWeight=  254.63*5.3009e-4/399398.0; }   // JZ4W : 361024 : 399400 : 254.63nb : 5.3009E-04      254.63
    else if (mcDSID==361025) { xsecWeight=   4.553*9.2325e-4/399300.0; }   // JZ5W : 361025 : 399300 :  4.553nb : 9.2325E-04
  }
  else {
    muValue     = eventInfo->auxdata<float>("averageInteractionsPerCrossing");
    muValueCorr = eventInfo->auxdata<float>("corrected_averageInteractionsPerCrossing");
  }
  eventFeatures.muValue     = muValue;
  eventFeatures.muValueCorr = muValueCorr;

  // Trigger decision tool
  float triggerWeight = 1.0;
  if (m_trigDecTool) {

    bool pass_HLT_j100 = m_trigDecTool->isPassed("HLT_j100");
    bool pass_HLT_j110 = m_trigDecTool->isPassed("HLT_j110");
    bool pass_HLT_j150 = m_trigDecTool->isPassed("HLT_j150");
    bool pass_HLT_j175 = m_trigDecTool->isPassed("HLT_j175");
    bool pass_HLT_j200 = m_trigDecTool->isPassed("HLT_j200");
    bool pass_HLT_j260 = m_trigDecTool->isPassed("HLT_j260");
    bool pass_HLT_j300 = m_trigDecTool->isPassed("HLT_j300");
    bool pass_HLT_j320 = m_trigDecTool->isPassed("HLT_j320");
    bool pass_HLT_j360 = m_trigDecTool->isPassed("HLT_j360");
    bool pass_HLT_j380 = m_trigDecTool->isPassed("HLT_j380");
    bool pass_HLT_j400 = m_trigDecTool->isPassed("HLT_j400");
    bool pass_HLT_j420 = m_trigDecTool->isPassed("HLT_j420");
    bool pass_HLT_j460 = m_trigDecTool->isPassed("HLT_j460");

    bool pass_desd_jetA = (pass_HLT_j100 || pass_HLT_j110 || pass_HLT_j150 || pass_HLT_j175 || pass_HLT_j200 || pass_HLT_j260 || pass_HLT_j300);
    bool pass_desd_jetB = (pass_HLT_j320);
    bool pass_desd_jetC = (pass_HLT_j360 || pass_HLT_j380 || pass_HLT_j400);
    bool pass_desd_jetD = (pass_HLT_j420 && !pass_HLT_j460);
    bool pass_desd_jetE = (pass_HLT_j460);

    if      (pass_desd_jetE) { eventFeatures.trigCounter = 5; }
    else if (pass_desd_jetD) { eventFeatures.trigCounter = 4; }
    else if (pass_desd_jetC) { eventFeatures.trigCounter = 3; }
    else if (pass_desd_jetB) { eventFeatures.trigCounter = 2; }
    else if (pass_desd_jetA) { eventFeatures.trigCounter = 1; }
    else                     { eventFeatures.trigCounter = 0; }

    if (isMC) {
      int mcDSID = eventInfo->mcChannelNumber();
      if (mcDSID==361023) {
        if      (pass_desd_jetE) { triggerWeight = 0.0; }
        else if (pass_desd_jetD) { triggerWeight = 0.0; }
        else if (pass_desd_jetC) { triggerWeight = 178246.0/301.960; }
        else if (pass_desd_jetB) { triggerWeight = 178246.0/301.960; }
        else if (pass_desd_jetA) { triggerWeight = 178246.0/301.960; }
        else                     { triggerWeight = 0.0; }
      }
      if (mcDSID==361024) {
        if      (pass_desd_jetE) { triggerWeight = 7444.7; }
        else if (pass_desd_jetD) { triggerWeight = 7919.7; }
        else if (pass_desd_jetC) { triggerWeight = 14121.7; }
        else if (pass_desd_jetB) { triggerWeight = 73832.5; }
        else if (pass_desd_jetA) { triggerWeight = 0.0; }
        else                     { triggerWeight = 0.0; }
      }
    }

//    std::cout << pass_HLT_j100 << pass_HLT_j110 << pass_HLT_j150 << pass_HLT_j175 << pass_HLT_j200 << pass_HLT_j260 << pass_HLT_j300 << pass_HLT_j320 << pass_HLT_j360 << pass_HLT_j380 << pass_HLT_j400 << pass_HLT_j420 << pass_HLT_j460 << std::endl;

//sel_jet600 = 'AntiKt4EMTopoJets.JetConstitScaleMomentum_pt >= 600.*GeV'
//sel_jet800 = 'AntiKt4EMTopoJets.JetConstitScaleMomentum_pt >= 800.*GeV'
//sel_jet1000 = 'AntiKt4EMTopoJets.JetConstitScaleMomentum_pt >= 1000.*GeV'

//desd_jetA = '(HLT_j100 || HLT_j110 || HLT_j150 || HLT_j175 || HLT_j200 || HLT_j260 || HLT_j300)'
//desd_jetB = '( HLT_j320 )'
//desd_jetC = '( HLT_j360 ||  HLT_j380 || HLT_j400 )'
//desd_jetD = '( HLT_j420 && !HLT_j460 )'
//desd_jetE = '( HLT_j460 )'
//desd_jetF = '( HLT_j460 && count('+sel_jet600+')>0 && count('+sel_jet800+')==0 )'
//desd_jetG = '( HLT_j460 && count('+sel_jet800+')>0 && count('+sel_jet1000+')==0 )'
//desd_jetH = '( HLT_j460 && count('+sel_jet1000+')>0 )'


//    if (!(m_trigDecTool->isPassed("HLT_noalg_mb_L1RD3_FILLED"))) { return EL::StatusCode::SUCCESS; }
//    if (!(m_trigDecTool->isPassed("HLT_noalg_mb_L1RD3_FILLED"))) { return EL::StatusCode::SUCCESS; }

//    std::cout << "trig " << m_trigDecTool->isPassed("HLT_mb_mbts_L1MBTS_1") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("HLT_noalg_mbts_MBTS_1") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("HLT_noalg_mb_L1MBTS_1") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("HLT_noalg_mb_L1RD0_FILLED") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("HLT_noalg_mb_L1RD3_FILLED") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("L1_MBTS_1_1") << std::endl;
//    std::cout << "trig " << m_trigDecTool->isPassed("L1_RD0_FILLED") << std::endl;

//    auto chainGroup = m_trigDecTool->getChainGroup("HLT*");
//    for (auto &trig:chainGroup->getListOfTriggers()) {
//      auto cg = m_trigDecTool->getChainGroup(trig);
//      std::string thisTrig = trig;
//      std::cout << thisTrig << std::endl;
//      Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale (L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
//    }

  }


  // Trigger selection
//  if (eventFeatures.trigCounter==0) { return EL::StatusCode::SUCCESS; }


  // weight for histogram
  float histWeight = pileupWeight*xsecWeight*mcEvWeight*triggerWeight;


  //==========================================
  // Build a map of every cluster on a module
  //==========================================
  std::multimap<int, const xAOD::TrackMeasurementValidation*> clusterPerModule;
  for (auto cluster:*pixClusters) {
    clusterPerModule.insert(std::pair<int, const xAOD::TrackMeasurementValidation*>(buildPixelModuleID(cluster), cluster));
  } // end for loop over clusters

  //=========================
  // Dress clusters on track
  //=========================
  dressLayerHitContent( recoTracks, pixClusters, sctClusters); 

  //=========================================================================================================================
  // Build track-pair class for reconstructed tracks
  //=========================================================================================================================
  /*
   nHitsOnPixelLayer[0] ... IBL
   nHitsOnPixelLayer[1] ... b-layer
   nHitsOnPixelLayer[2] ... pixel-2 (barrel)
   nHitsOnPixelLayer[3] ... pixel-3 (barrel)
   nHitsOnPixelLayer[4] ... pixel-0 (endcap)
   nHitsOnPixelLayer[5] ... pixel-1 (endcap)
   nHitsOnPixelLayer[6] ... pixel-2 (endcap)

   Radius : 33.0 [mm]   ... IBL at Stave position
            50.5 [mm]   ... b-layer
            88.5 [mm]   ... pixel-2
           122.5 [mm]   ... pixel-3
  */

//  int cut_ntrackIsolation = 0;                       // 0 track except target track.
  int cut_ntrackIsolation = 100;                       // 0 track except target track.
  double cut_leadingTrkPt = 1.0;                    // [GeV]
  double cut_subleadingTrkPt = 1.0;                  // [GeV]
  double cut_maxTrackSeparation = 0.1;               // [--]

  xAOD::TrackParticleContainer* twoRecoTracks = new xAOD::TrackParticleContainer();
  xAOD::TrackParticleAuxContainer* twoRecoTracksAux = new xAOD::TrackParticleAuxContainer();
  twoRecoTracks->setStore(twoRecoTracksAux); //< Connect the two

  // IBL - Barrel
  int targetLayer = 0;
//  int targetLayer = 1;
  float checkPt = -1.0; float checkEta = -100.0; float checkPhi = -100.0;
  for (xAOD::TrackParticleContainer::const_iterator twoRecoTrk_itr=recoTracks->begin(); twoRecoTrk_itr!=recoTracks->end(); twoRecoTrk_itr++) {

    if ((*twoRecoTrk_itr)->pt()/1000.0<cut_leadingTrkPt) { continue; } // pt cut
    if (!checkQuality(*twoRecoTrk_itr,targetLayer+1)) { continue; }    // Quality cut

    //===================================
    // Find the closest track on surface
    //===================================
    int nonIsoTrk = 0;
    double mindR = 10000.0;
    xAOD::TrackParticleContainer::const_iterator twoRecoSiblingTrk_itr;
    for (xAOD::TrackParticleContainer::const_iterator recoTrk_itr=recoTracks->begin(); recoTrk_itr!=recoTracks->end(); recoTrk_itr++) {
      // skip same track
      if (fabs((*recoTrk_itr)->pt()-(*twoRecoTrk_itr)->pt())<1e-4 && fabs((*recoTrk_itr)->phi()-(*twoRecoTrk_itr)->phi())<1e-4 && fabs((*recoTrk_itr)->eta()-(*twoRecoTrk_itr)->eta())<1e-4) { continue; }
      if ((*recoTrk_itr)->pt()/1000.0<cut_subleadingTrkPt) { continue; } // pt cut
      if (!checkQuality(*recoTrk_itr,targetLayer+1)) { continue; }       // Quality cut

      double dR = deltaR(*twoRecoTrk_itr,*recoTrk_itr);
      if (dR<mindR) { mindR = dR; twoRecoSiblingTrk_itr = recoTrk_itr; }
      if (dR<0.1) { nonIsoTrk++; }
    }

    if (mindR>0.0 && mindR<cut_maxTrackSeparation && nonIsoTrk<=1+cut_ntrackIsolation) {
      // avoid double counting tracks.
      if (fabs(checkPt-(*twoRecoSiblingTrk_itr)->pt())>1e-4 && fabs(checkEta-(*twoRecoSiblingTrk_itr)->eta())>1e-4 && fabs(checkPhi-(*twoRecoSiblingTrk_itr)->phi())) {
        // Find cluster on surface
        int nRefMSOSLayer = 0;
        int nRefClusLayer = 0;
        const xAOD::TrackStateValidation *msosRef[4];
        const xAOD::TrackMeasurementValidation *pixRefClus[4];
        findPixelClusterOnLayer(*twoRecoTrk_itr, pixClusters, targetLayer, &nRefMSOSLayer, &msosRef[0], &nRefClusLayer, &pixRefClus[0]);

        int nTarMSOSLayer = 0;
        int nTarClusLayer = 0;
        const xAOD::TrackStateValidation *msosTar[4];
        const xAOD::TrackMeasurementValidation *pixTarClus[4];
        findPixelClusterOnLayer(*twoRecoSiblingTrk_itr, pixClusters, targetLayer, &nTarMSOSLayer, &msosTar[0], &nTarClusLayer, &pixTarClus[0]);

        int TargetLayer = targetLayer;
        int clusPatternIndex = -1;
        double minMSOSLocalDistance[2]  = {-1.0,-1.0};
        double minClusLocalDistance[2]  = {-1.0,-1.0};
        double minClusGlobalDistance[3] = {-1.0,-1.0,-1.0};

        if (nRefMSOSLayer==0 || nTarMSOSLayer==0 || nRefMSOSLayer==3 || nTarMSOSLayer==3) {
//          std::cout << (*twoRecoSiblingTrk_itr)->phi() << " " << (*twoRecoSiblingTrk_itr)->eta() << std::endl;
        }
        else {
          // Fill cluster pattern
          clusterPatternClassification(nRefMSOSLayer,msosRef,nRefClusLayer,pixRefClus,
                                       nTarMSOSLayer,msosTar,nTarClusLayer,pixTarClus,
                                       &clusPatternIndex,&minMSOSLocalDistance[0],&minClusLocalDistance[0],&minClusGlobalDistance[0]);

          // Copy this track to the new container:
          xAOD::TrackParticle* tp = new xAOD::TrackParticle();
          tp->makePrivateStore( **twoRecoTrk_itr);

          static SG::AuxElement::Decorator<xAOD::TrackParticleContainer::const_iterator> twoRecoSiblingTrkIter("twoRecoSiblingTrkIter");
          twoRecoSiblingTrkIter(*tp) = twoRecoSiblingTrk_itr;

          static SG::AuxElement::Decorator<int> nRefHitAtLayerOfInterest("nRefHitAtLayerOfInterest");
          static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> RefMSOS1("RefMSOS1");
          static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> RefMSOS2("RefMSOS2");
          static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelRefClus1("PixelRefClus1");
          static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelRefClus2("PixelRefClus2");
          nRefHitAtLayerOfInterest(*tp) = nRefClusLayer;
          RefMSOS1(*tp) = msosRef[0];
          RefMSOS2(*tp) = msosRef[1];
          PixelRefClus1(*tp) = pixRefClus[0];
          PixelRefClus2(*tp) = pixRefClus[1];

          static SG::AuxElement::Decorator<int> nTarHitAtLayerOfInterest("nTarHitAtLayerOfInterest");
          static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> TarMSOS1("TarMSOS1");
          static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> TarMSOS2("TarMSOS2");
          static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelTarClus1("PixelTarClus1");
          static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelTarClus2("PixelTarClus2");
          nTarHitAtLayerOfInterest(*tp) = nTarClusLayer;
          TarMSOS1(*tp) = msosTar[0];
          TarMSOS2(*tp) = msosTar[1];
          PixelTarClus1(*tp) = pixTarClus[0];
          PixelTarClus2(*tp) = pixTarClus[1];

          static SG::AuxElement::Decorator<int>   TargetLayerOfInterest("TargetLayerOfInterest");
          static SG::AuxElement::Decorator<int>   ClusterPatternIndex("ClusterPatternIndex");
          static SG::AuxElement::Decorator<float> minMSOSLocalDistanceAtLayer("minMSOSLocalDistanceAtLayer");
          static SG::AuxElement::Decorator<float> minClusLocalDistanceAtLayer("minClusLocalDistanceAtLayer");
          static SG::AuxElement::Decorator<float> minClusGlobalDistanceAtLayer("minClusGlobalDistanceAtLayer");
          TargetLayerOfInterest(*tp)        = TargetLayer;
          ClusterPatternIndex(*tp)          = clusPatternIndex;
          minMSOSLocalDistanceAtLayer(*tp)  = (float)TMath::Sqrt(minMSOSLocalDistance[0]*minMSOSLocalDistance[0]+minMSOSLocalDistance[1]*minMSOSLocalDistance[1]);
          minClusLocalDistanceAtLayer(*tp)  = (float)TMath::Sqrt(minClusLocalDistance[0]*minClusLocalDistance[0]+minClusLocalDistance[1]*minClusLocalDistance[1]);
          minClusGlobalDistanceAtLayer(*tp) = (float)TMath::Sqrt(minClusGlobalDistance[0]*minClusGlobalDistance[0]+minClusGlobalDistance[1]*minClusGlobalDistance[1]+minClusGlobalDistance[2]*minClusGlobalDistance[2]);

          static SG::AuxElement::Decorator<float> minMSOSLocalDistanceAtLayerX("minMSOSLocalDistanceAtLayerX");
          static SG::AuxElement::Decorator<float> minClusLocalDistanceAtLayerX("minClusLocalDistanceAtLayerX");
          static SG::AuxElement::Decorator<float> minClusGlobalDistanceAtLayerX("minClusGlobalDistanceAtLayerX");
          minMSOSLocalDistanceAtLayerX(*tp)  = (float)minMSOSLocalDistance[0];
          minClusLocalDistanceAtLayerX(*tp)  = (float)minClusLocalDistance[0];
          minClusGlobalDistanceAtLayerX(*tp) = (float)minClusGlobalDistance[0];

          static SG::AuxElement::Decorator<float> minMSOSLocalDistanceAtLayerY("minMSOSLocalDistanceAtLayerY");
          static SG::AuxElement::Decorator<float> minClusLocalDistanceAtLayerY("minClusLocalDistanceAtLayerY");
          static SG::AuxElement::Decorator<float> minClusGlobalDistanceAtLayerY("minClusGlobalDistanceAtLayerY");
          minMSOSLocalDistanceAtLayerY(*tp)  = (float)minMSOSLocalDistance[1];
          minClusLocalDistanceAtLayerY(*tp)  = (float)minClusLocalDistance[1];
          minClusGlobalDistanceAtLayerY(*tp) = (float)minClusGlobalDistance[1];

          static SG::AuxElement::Decorator<float>   z0Corr1("z0Corr1");
          static SG::AuxElement::Decorator<float>   z0Corr2("z0Corr2");
          z0Corr1(*tp) = (*twoRecoTrk_itr)->z0()+(*twoRecoTrk_itr)->vz()-vtxZ;
          z0Corr2(*tp) = (*twoRecoSiblingTrk_itr)->z0()+(*twoRecoSiblingTrk_itr)->vz()-vtxZ;

          // look for jet matches with tracks.
          double mnDRJ1 = 100.0; double ptJ1 = 0.0;
          double mnDRJ2 = 100.0; double ptJ2 = 0.0;
          for (xAOD::JetContainer::const_iterator jet_itr=jetAntiKt4EMTopo->begin(); jet_itr!=jetAntiKt4EMTopo->end(); jet_itr++) {
            if ((*jet_itr)->pt()/1000.0>10.0) {
              double dR1 = deltaR(*twoRecoTrk_itr,*jet_itr);
              double dR2 = deltaR(*twoRecoSiblingTrk_itr,*jet_itr);
              if (dR1<mnDRJ1) { mnDRJ1=dR1; ptJ1=(*jet_itr)->pt(); }
              if (dR2<mnDRJ2) { mnDRJ2=dR2; ptJ2=(*jet_itr)->pt(); }
            }
          }

          static SG::AuxElement::Decorator<float> mindRJet1("mindRJet1");
          static SG::AuxElement::Decorator<float> mindRJet2("mindRJet2");
          static SG::AuxElement::Decorator<float> pTJet1("pTJet1");
          static SG::AuxElement::Decorator<float> pTJet2("pTJet2");
          mindRJet1(*tp) = (float)mnDRJ1;
          mindRJet2(*tp) = (float)mnDRJ2;
          pTJet1(*tp)    = (float)ptJ1;
          pTJet2(*tp)    = (float)ptJ2;

          twoRecoTracks->push_back( tp );

          // avoid double counting
          checkPt = (*twoRecoTrk_itr)->pt(); checkEta = (*twoRecoTrk_itr)->eta(); checkPhi = (*twoRecoTrk_itr)->phi();
        }
      }
    }
  }
  //=========================================================================================================================

  //---------------------
  // FILL PLOTS!!!
  //---------------------
  //eventFeatures.dumpInfo();
  eventHist_all->FillHists(eventFeatures,histWeight); 

  typedef std::multimap<int, const xAOD::TrackParticle* > barcodeMap;

  //=========================================================================================================================
  // Fill histograms for track-pair for reconstructed tracks
  //=========================================================================================================================
  barcodeMap twoRecoTracksBarcodeMap;
  barcodeMap neighborBarcodeMap;
  for (xAOD::TrackParticleContainer::const_iterator twoRecoTrk_itr=twoRecoTracks->begin(); twoRecoTrk_itr!=twoRecoTracks->end(); twoRecoTrk_itr++) {

    int TargetLayer  = (*twoRecoTrk_itr)->auxdata<int>("TargetLayerOfInterest");

    xAOD::TrackParticleContainer::const_iterator neighbor_itr = (*twoRecoTrk_itr)->auxdata<xAOD::TrackParticleContainer::const_iterator>("twoRecoSiblingTrkIter");

    static SG::AuxElement::Decorator<xAOD::TrackParticleContainer::const_iterator> twoRecoSiblingTrkIter("twoRecoSiblingTrkIter");

    int nRefLayer = (*twoRecoTrk_itr)->auxdata<int>("nRefHitAtLayerOfInterest");
    const xAOD::TrackStateValidation *msosRef1 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("RefMSOS1");
    const xAOD::TrackStateValidation *msosRef2 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("RefMSOS2");
    const xAOD::TrackMeasurementValidation *pixRefClus1 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus1");
    const xAOD::TrackMeasurementValidation *pixRefClus2 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus2");

    int nTarLayer = (*twoRecoTrk_itr)->auxdata<int>("nTarHitAtLayerOfInterest");
    const xAOD::TrackStateValidation *msosTar1 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("TarMSOS1");
    const xAOD::TrackStateValidation *msosTar2 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("TarMSOS2");
    const xAOD::TrackMeasurementValidation *pixTarClus1 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelTarClus1");
    const xAOD::TrackMeasurementValidation *pixTarClus2 = (*twoRecoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelTarClus2");

    float minMSOSLocalDistance = (*twoRecoTrk_itr)->auxdata<float>("minMSOSLocalDistanceAtLayer");

    // Associate truth particle in the reference track
//STSTST      addTruthDistancePerLayer( *twoRecoTrk_itr, clusterPerModule );

    const xAOD::TruthParticle *twoRecoTracks_truth = xAOD::TrackHelper::truthParticle(*twoRecoTrk_itr);
    if( twoRecoTracks_truth ) {
      twoRecoTracksBarcodeMap.insert( std::pair<int, const xAOD::TrackParticle*> (twoRecoTracks_truth->barcode(), *twoRecoTrk_itr) ); 
      (*twoRecoTrk_itr)->auxdecor<int>("barcode") = twoRecoTracks_truth->barcode();
    } else {
      (*twoRecoTrk_itr)->auxdecor<int>("barcode") = -888;
    }
    (*twoRecoTrk_itr)->auxdecor<float>("motherPt") = -1.0;
    (*twoRecoTrk_itr)->auxdecor<float>("motherR")  = -1.0;

    // Associate truth particle in the neighbor track
//STSTST      addTruthDistancePerLayer( *neighbor_itr, clusterPerModule );
    const xAOD::TruthParticle *neighbor_truth = xAOD::TrackHelper::truthParticle(*neighbor_itr);
    if( neighbor_truth ) {
      neighborBarcodeMap.insert( std::pair<int, const xAOD::TrackParticle*> (neighbor_truth->barcode(), *neighbor_itr) ); 
      (*neighbor_itr)->auxdecor<int>("barcode") = neighbor_truth->barcode();
    } else {
      (*neighbor_itr)->auxdecor<int>("barcode") = -888;
    }
    (*neighbor_itr)->auxdecor<float>("motherPt") = -1.0;
    (*neighbor_itr)->auxdecor<float>("motherR")  = -1.0;

    // Now fill track histogram
    twotrkHist_reco_pair -> FillHists((*twoRecoTrk_itr),histWeight);
    trkHist_reco_pair -> FillHists((*twoRecoTrk_itr),histWeight);
    trkHist_reco_pair -> FillHists((*neighbor_itr),histWeight);
    if ((*twoRecoTrk_itr)->pt()>(*neighbor_itr)->pt()) {
      trkHist_reco_pair_leading    -> FillHists((*twoRecoTrk_itr),histWeight);
      if (nRefLayer>0) { 
        clusHist_reco_pair_leading -> FillHists(pixRefClus1,histWeight); 
        if (minMSOSLocalDistance<0.4) { clusHist_leading_low  -> FillHists(pixRefClus1,histWeight); }
        else                          { clusHist_leading_high -> FillHists(pixRefClus1,histWeight); }
      }
      if (nRefLayer>1) { 
        clusHist_reco_pair_leading -> FillHists(pixRefClus2,histWeight); 
        if (minMSOSLocalDistance<0.4) { clusHist_leading_low  -> FillHists(pixRefClus2,histWeight); }
        else                          { clusHist_leading_high -> FillHists(pixRefClus2,histWeight); }
      }

      trkHist_reco_pair_subleading -> FillHists((*neighbor_itr),histWeight);
      if (nTarLayer>0) { clusHist_reco_pair_subleading -> FillHists(pixTarClus1,histWeight); }
      if (nTarLayer>1) { clusHist_reco_pair_subleading -> FillHists(pixTarClus2,histWeight); }
    }
    else {
      trkHist_reco_pair_subleading -> FillHists((*twoRecoTrk_itr),histWeight);
      if (nRefLayer>0) { clusHist_reco_pair_subleading -> FillHists(pixRefClus1,histWeight); }
      if (nRefLayer>1) { clusHist_reco_pair_subleading -> FillHists(pixRefClus2,histWeight); }

      trkHist_reco_pair_leading    -> FillHists((*neighbor_itr),histWeight);
      if (nTarLayer>0) { 
        clusHist_reco_pair_leading -> FillHists(pixTarClus1,histWeight); 
        if (minMSOSLocalDistance<0.4) { clusHist_leading_low  -> FillHists(pixTarClus1,histWeight); }
        else                          { clusHist_leading_high -> FillHists(pixTarClus1,histWeight); }
      }
      if (nTarLayer>1) { 
        clusHist_reco_pair_leading -> FillHists(pixTarClus2,histWeight); 
        if (minMSOSLocalDistance<0.4) { clusHist_leading_low  -> FillHists(pixTarClus2,histWeight); }
        else                          { clusHist_leading_high -> FillHists(pixTarClus2,histWeight); }
      }
    }
  }


  for (xAOD::JetContainer::const_iterator jet_itr=jetAntiKt4EMTopo->begin(); jet_itr!=jetAntiKt4EMTopo->end(); jet_itr++) {
    if ((*jet_itr)->pt()/1000.0>10.0 && TMath::Abs((*jet_itr)->eta())<2.5) {
      double minMSOSLocalX = 100.0;
      double minMSOSLocalY = 100.0;
      xAOD::TrackParticleContainer::const_iterator twoRecoTrk_minX;
      xAOD::TrackParticleContainer::const_iterator twoRecoTrk_minY;
      for (xAOD::TrackParticleContainer::const_iterator twoRecoTrk_itr=twoRecoTracks->begin(); twoRecoTrk_itr!=twoRecoTracks->end(); twoRecoTrk_itr++) {
        float mnDRJ1 = (*twoRecoTrk_itr)->auxdata<float>("mindRJet1");
        float mnDRJ2 = (*twoRecoTrk_itr)->auxdata<float>("mindRJet2");
        float ptJ1   = (*twoRecoTrk_itr)->auxdata<float>("pTJet1");
        float ptJ2   = (*twoRecoTrk_itr)->auxdata<float>("pTJet2");

        if (TMath::Abs(ptJ1-(*jet_itr)->pt())<0.01 && TMath::Abs(ptJ2-(*jet_itr)->pt())<0.01) {
          if (TMath::Abs(ptJ1-ptJ2)<0.01 && mnDRJ1<0.4 && mnDRJ2<0.4) {
            float minMSOSLocalDistanceX = (*twoRecoTrk_itr)->auxdata<float>("minMSOSLocalDistanceAtLayerX");
            float minMSOSLocalDistanceY = (*twoRecoTrk_itr)->auxdata<float>("minMSOSLocalDistanceAtLayerY");
            if (minMSOSLocalDistanceX<minMSOSLocalX) {
              minMSOSLocalX = minMSOSLocalDistanceX;
              twoRecoTrk_minX = twoRecoTrk_itr;
            }
            if (minMSOSLocalDistanceX<minMSOSLocalY) {
              minMSOSLocalY = minMSOSLocalDistanceY;
              twoRecoTrk_minY = twoRecoTrk_itr;
            }
          }
        }
      }
      if (minMSOSLocalX<100.0) { twotrkHist_minX -> FillHists((*twoRecoTrk_minX),histWeight); }
      if (minMSOSLocalY<100.0) { twotrkHist_minY -> FillHists((*twoRecoTrk_minY),histWeight); }
    }
  }




  //---------------------
  // CLEAN UP
  //---------------------
  delete recoTracksShallowCopyPair.first;
  delete recoTracksShallowCopyPair.second;
  delete pixClustersPair.first;
  delete pixClustersPair.second;
  delete sctClustersPair.first;
  delete sctClustersPair.second;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TwoTrkEffStudy :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TwoTrkEffStudy :: finalize ()
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

  //clean up the memory
  if (m_trigDecTool) {
    delete m_trigDecTool;
    m_trigDecTool = 0;
  }

  if (m_configTool) {
    delete m_configTool;
    m_configTool = 0;
  }

  // GRL
  if (m_grl) {
    delete m_grl;
    m_grl = 0;
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TwoTrkEffStudy :: histFinalize ()
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

int TwoTrkEffStudy::buildPixelModuleID(const xAOD::TrackMeasurementValidation* prd) const {
  int bec        = prd->auxdataConst< int >("bec");   // -4,4
  int layer      = prd->auxdataConst< int >("layer");      // 0,4
  int eta_module = prd->auxdataConst< int >("eta_module"); // -20,20
  int phi_module = prd->auxdataConst< int >("phi_module"); // 0-50
  return bec * 1000000 + layer * 100000 + phi_module * 100 + eta_module + 50; 
}

void TwoTrkEffStudy::addTruthDistancePerLayer(const xAOD::TrackParticle* tp, std::multimap<int, const xAOD::TrackMeasurementValidation*>& clustersOnAModule) const {
  //Calculates the distance to the neaeast true paricle
  //Note the track needs to have a cluster on the layer in question 
  //If no other cluster is on that module then you will get the default value
  std::vector<float> resultsR(4,-9999);
  std::vector<float> resultsX(4,-9999);
  std::vector<float> resultsY(4,-9999);

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
//  static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change
//  static const char* measurementNames = "IDTIDE1msosLink";   //Note the prefix could change
  static const char* measurementNames = "msosLink";   //Note the prefix could change

  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack( measurementNames );

  // Check if there are MSOS attached
  if (!tp->isAvailable<MeasurementsOnTrack>(measurementNames)) { return; }
  // Get the MSOS's
  const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack(*tp);

  // Loop over track TrackStateValidation's
  const xAOD::TruthParticle*  truth = xAOD::TrackHelper::truthParticle(tp);
  if (!truth) { return; }
  
  for (MeasurementsOnTrackIter msos_iter=measurementsOnTrack.begin(); msos_iter!=measurementsOnTrack.end(); ++msos_iter) {  
    //Check if the element link is valid
    if (!(*msos_iter).isValid()) { continue; }
    
    const xAOD::TrackStateValidation* msos = *(*msos_iter); 
   
    // The detector type enum is defined in TrkEventPrimitives/TrackStateDefs.h
    // The mesurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
    // We will port the enums to xAOD world in the future
    if (msos->detType()!=1 || msos->type()==6) { continue; }  //its pixel and not a hole
    
    //Get pixel cluster
    if (msos->trackMeasurementValidationLink().isValid() && *(msos->trackMeasurementValidationLink())) {
      const xAOD::TrackMeasurementValidation* pixelCluster =  *(msos->trackMeasurementValidationLink());        
              
      // if it is not in the barrel
      if (pixelCluster->auxdata<int>("bec")!=0) { continue; }
      
      int layer = pixelCluster->auxdata<int>("layer");
      
      static SG::AuxElement::ConstAccessor<std::vector<int>>   acc_sihit_barcode("sihit_barcode");
      static SG::AuxElement::ConstAccessor<std::vector<float>> acc_sihit_startPosX("sihit_startPosX");
      static SG::AuxElement::ConstAccessor<std::vector<float>> acc_sihit_startPosY("sihit_startPosY");
      static SG::AuxElement::ConstAccessor<std::vector<float>> acc_sihit_endPosX("sihit_endPosX");
      static SG::AuxElement::ConstAccessor<std::vector<float>> acc_sihit_endPosY("sihit_endPosY");

      const std::vector<int>   &sihit_barcode   = acc_sihit_barcode(*pixelCluster);
      const std::vector<float> &sihit_startPosX = acc_sihit_startPosX(*pixelCluster);
      const std::vector<float> &sihit_startPosY = acc_sihit_startPosY(*pixelCluster);
      const std::vector<float> &sihit_endPosX   = acc_sihit_endPosX(*pixelCluster);
      const std::vector<float> &sihit_endPosY   = acc_sihit_endPosY(*pixelCluster);

      int myTP_Location(-1); 
      for (unsigned int i=0; i<sihit_barcode.size(); ++i) {
        if (sihit_barcode[i]==truth->barcode()) {
          myTP_Location = i;
          break;
        }
      }
      
      //Cluster is incorrectly assigned to the track
      if (myTP_Location==-1) { continue; }

      float myhitX = (sihit_endPosX[myTP_Location]+sihit_startPosX[myTP_Location])*0.5;
      float myhitY = (sihit_endPosY[myTP_Location]+sihit_startPosY[myTP_Location])*0.5;

      std::vector<int> otherParticles;
      for (unsigned int i=0; i<sihit_barcode.size(); ++i) {
        if ((int)i==myTP_Location) { continue; }
        //only consider primary particles
        if (sihit_barcode[i]>0 && sihit_barcode[i]<200000) {
          otherParticles.push_back(i);
        }
      }

      //Multiple paricles in the cluster which is the closest one 
      if (otherParticles.size()>0) {
        float minDeltaR2(100000);
        int closest(0);
        for (auto i:otherParticles) {
          float dx = myhitX-(sihit_endPosX[i]+sihit_startPosX[i])*0.5;
          float dy = myhitY-(sihit_endPosY[i]+sihit_startPosY[i])*0.5;
          float deltaR2 = dx*dx+dy*dy;
          if (deltaR2<minDeltaR2) {
            closest = (int)i;
            minDeltaR2 = deltaR2;
          }
        }
        resultsR[layer] = minDeltaR2;
        resultsX[layer] = myhitX-(sihit_endPosX[closest]+sihit_startPosX[closest])*0.5;
        resultsY[layer] = myhitY-(sihit_endPosY[closest]+sihit_startPosY[closest])*0.5;
        continue;
      }

      //Get all of the clusters on that module
      int clusterModuleID = buildPixelModuleID(pixelCluster); 
      // equal_range(**) returns pair<iterator,iterator> representing the range
      // of element with key **
      
      auto pair_it =  clustersOnAModule.equal_range(clusterModuleID);
    
      //std::cout << "Layer :" <<  layer << " Cluster does not contain more than one truth particle -- look for other clusters on the module" << std::endl;
      //std::cout << "Module ID: " << clusterModuleID << " other clusters " <<  clustersOnAModule.count(clusterModuleID) << "  "<< clustersOnAModule.size() << std::endl;

      //Find the closest cluster in the same module;
      const xAOD::TrackMeasurementValidation* closestCluster = 0;
      float minDeltaR2(1.e20);
      
      //Copy the barcode vector so we can use it a comparison later
      const std::vector<int> sihit_barcode_Copy = sihit_barcode;
      // Loop through range of maps of key **
      for (auto it = pair_it.first; it!=pair_it.second; ++it) {
        // Note you can not use the pointers directly as they may change 
        if (it->second->identifier()==pixelCluster->identifier()) { continue; }
        
        // Double check that if the identifier differ but if they contain exactly the same truth particles continue -- split cluster 
        const std::vector<int> &sihit_barcode_Other = acc_sihit_barcode(*it->second);
        if (sihit_barcode_Copy==sihit_barcode_Other) { continue; }
        
        float dx = pixelCluster->localX()-it->second->localX();  
        float dy = pixelCluster->localY()-it->second->localY();
        float deltaR2 = dx*dx+dy*dy;
        if (deltaR2<minDeltaR2) {
          closestCluster = it->second;
          minDeltaR2 = deltaR2;
        }    
      }
       
      // NO other cluster on the Module
      if (!closestCluster) {
        //std::cout << "No cluster found on module" << std::endl;
        continue;
      } 

      const std::vector<int>   &sihit_barcode_OC   = acc_sihit_barcode(*closestCluster);
      const std::vector<float> &sihit_startPosX_OC = acc_sihit_startPosX(*closestCluster);
      const std::vector<float> &sihit_startPosY_OC = acc_sihit_startPosY(*closestCluster);
      const std::vector<float> &sihit_endPosX_OC   = acc_sihit_endPosX(*closestCluster);
      const std::vector<float> &sihit_endPosY_OC   = acc_sihit_endPosY(*closestCluster);
      
      otherParticles.clear();  
      for (unsigned int i=0; i<sihit_barcode_OC.size(); ++i) {
        if (sihit_barcode_OC[i]==truth->barcode()) {
          //std::cout << "Truth particle has same barcode as the one i am interested in" << std::endl;  
          //float dx = myhitX - (sihit_endPosX_OC[i] + sihit_startPosX_OC[i]) *0.5;
          //float dy = myhitY - (sihit_endPosY_OC[i] + sihit_startPosY_OC[i]) *0.5;
          //std::cout << "dx " << dx << " dy " <<  dy << std::endl; 
          continue;
        }
        //only consider primary particles
        if (sihit_barcode_OC[i]>0 && sihit_barcode_OC[i]<200000) {
          otherParticles.push_back(i);
        }
      }

      //Multiple paricles in the cluster which is the closeest one 
      if (otherParticles.size()>0) {
        float minDeltaR2(100000);
        int closest(0);
        for (auto i:otherParticles) {
          float dx = myhitX-(sihit_endPosX_OC[i]+sihit_startPosX_OC[i])*0.5;
          float dy = myhitY-(sihit_endPosY_OC[i]+sihit_startPosY_OC[i])*0.5;
          float deltaR2 = dx*dx+dy*dy;
          if (deltaR2<minDeltaR2) {
            closest =  (int)i;
            minDeltaR2 = deltaR2;
          }
        }
        //std::cout << "Found another truth particle in another cluster" << std::endl;  
        
        resultsR[layer] = minDeltaR2;
        resultsX[layer] = myhitX-(sihit_endPosX_OC[closest]+sihit_startPosX_OC[closest])*0.5;
        resultsY[layer] = myhitY-(sihit_endPosY_OC[closest]+sihit_startPosY_OC[closest])*0.5;
        continue;
      } 
      else {
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

/*
void TwoTrkEffStudy::filldetIdtoMsosMap(const xAOD::TrackParticleContainer* tc, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const {
	for (auto trk:(*tc)) {
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
*/

/*
void TwoTrkEffStudy::MinDistanceOtherTrack(const xAOD::TrackParticle* tp, const xAOD::TrackParticleContainer* tc, std::multimap<uint64_t, const xAOD::TrackStateValidation*>& detIdtoMsosMap) const {

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
*/

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

/*
}
*/

/*
void TwoTrkEffStudy::fillPixelClusterHists(const xAOD::TrackStateValidation* msos, xAOD::TrackMeasurementValidationContainer* pixClusters, std::vector<ClusterHistManager*> &clusHists) {

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

  for (auto histIt = clusHists.begin(); histIt != clusHists.end(); ++histIt){
    (*histIt)->FillHists( pixprd );
    // fill different components of histogram vector based on prd detType
  }
} // fillPixelClusterHists
*/

void TwoTrkEffStudy::countContributingParticles(xAOD::TrackMeasurementValidation* prd) {
  std::vector<int> bcodes;
  int nPrimaryParticles(0);
  const std::vector< std::vector<int> > sdo_barcodes = prd->auxdata<std::vector<std::vector<int>>>("sdo_depositsBarcode");
  for (unsigned int i=0; i<sdo_barcodes.size(); i++) { 
    for (unsigned int j=0; j<sdo_barcodes.at(i).size(); j++) { 
      int bc(sdo_barcodes.at(i).at(j));
      std::vector<int>::iterator it = find(bcodes.begin(), bcodes.end(), bc);
      if (it==bcodes.end()) { // if at end then new barcode
        bcodes.push_back(bc);
        if (bc<=10e3 || bc>=200e3) { continue; }
        nPrimaryParticles++;
      }
    }
  } // loop over sdo barcodes
  prd->auxdata<int>("nParticles") = bcodes.size();
  prd->auxdata<int>("nPrimaryParticles") = nPrimaryParticles;
  prd->auxdata<std::vector<int> >("SDOBarcodes") = bcodes;
}


void TwoTrkEffStudy::dressLayerHitContent(xAOD::TrackParticleContainer* trks, xAOD::TrackMeasurementValidationContainer* pixClusters, xAOD::TrackMeasurementValidationContainer* sctClusters) {

  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
//  static const char* measurementNames = "IDDET1_msosLink";   //Note the prefix could change
//  static const char* measurementNames = "IDTIDE1msosLink";   //Note the prefix could change
  static const char* measurementNames = "msosLink";   //Note the prefix could change

  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack(measurementNames);

  // for each track, fill in vector for how many particles are contributing to each hit 
  xAOD::TrackParticleContainer::iterator trk_itr = trks->begin();
  xAOD::TrackParticleContainer::iterator trk_end = trks->end();

  for (xAOD::TrackParticleContainer::iterator trk_itr=trks->begin(); trk_itr!=trks->end(); ++trk_itr) {
    std::vector<int> nPartPix    = { -1, -1, -1, -1 };
    std::vector<int> nPriPartPix = { -1, -1, -1, -1 };
    std::vector<int> nPartSCT    = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> nPriPartSCT = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::vector<int> nHitsOnPixelLayer = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector<int> nHitsOnSCTLayer   = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };	

    (*trk_itr)->auxdata<float>("IBLCharge")    = -1;
    (*trk_itr)->auxdata<int>("IBLHitSize")     = -1;
    (*trk_itr)->auxdata<int>("IBLHitSizeZ")    = -1;
    (*trk_itr)->auxdata<int>("IBLHitSizePhi")  = -1;
    (*trk_itr)->auxdata<int>("nIBLMergedHits") = 0;
    (*trk_itr)->auxdata<float>("IBLlocalPhi")    = -100.0;
    (*trk_itr)->auxdata<float>("IBLlocalTheta")  = -100.0;
    (*trk_itr)->auxdata<float>("IBLbiasedResidualX")   = -100.0;
    (*trk_itr)->auxdata<float>("IBLbiasedResidualY")   = -100.0;
    (*trk_itr)->auxdata<float>("IBLunbiasedResidualX") = -100.0;
    (*trk_itr)->auxdata<float>("IBLunbiasedResidualY") = -100.0;
    (*trk_itr)->auxdata<float>("IBLResidualX")   = -100.0;
    (*trk_itr)->auxdata<float>("IBLResidualY")   = -100.0;

    (*trk_itr)->auxdata<float>("BLayerCharge")    = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSize")     = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSizeZ")    = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSizePhi")  = -1;
    (*trk_itr)->auxdata<float>("BLayerlocalPhi")    = -100.0;
    (*trk_itr)->auxdata<float>("BLayerlocalTheta")  = -100.0;
    (*trk_itr)->auxdata<float>("BLayerbiasedResidualX")   = -100.0;
    (*trk_itr)->auxdata<float>("BLayerbiasedResidualY")   = -100.0;
    (*trk_itr)->auxdata<float>("BLayerunbiasedResidualX") = -100.0;
    (*trk_itr)->auxdata<float>("BLayerunbiasedResidualY") = -100.0;
    (*trk_itr)->auxdata<float>("BLayerResidualX")   = -100.0;
    (*trk_itr)->auxdata<float>("BLayerResidualY")   = -100.0;

    // Get the MSOS's
    if (!(*trk_itr)->isAvailable<MeasurementsOnTrack>(measurementNames)) { 
      (*trk_itr)->auxdata<std::vector<int>>("pixClusPartContent")    = nPartPix;
      (*trk_itr)->auxdata<std::vector<int>>("pixClusPriPartContent") = nPriPartPix;
      (*trk_itr)->auxdata<std::vector<int>>("sctClusPartContent")    = nPartSCT;
      (*trk_itr)->auxdata<std::vector<int>>("sctClusPriPartContent") = nPriPartSCT;
      (*trk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer") = nHitsOnPixelLayer;
      (*trk_itr)->auxdata<std::vector<int>>("nHitsOnSCTLayer")   = nHitsOnSCTLayer;
      continue; 
    }

    const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack(*(*trk_itr));
    // Loop over track TrackStateValidation's
    for (MeasurementsOnTrackIter msos_iter=measurementsOnTrack.begin(); msos_iter!=measurementsOnTrack.end(); ++msos_iter) {  

      //Check if the element link is valid
      if (!(*msos_iter).isValid()) { continue; }

      const xAOD::TrackStateValidation* msos = *(*msos_iter); 
      if (!msos->trackMeasurementValidationLink().isValid()) { continue; }
      if (!(*(msos->trackMeasurementValidationLink())))      { continue; }

      // The detector type enum is defined in TrkEventPrimitives/TrackStateDefs.h
      // The mesurement state on surface type is defined in the header file TrkTrack/TrackStateOnSurface.h
      // We will port the enums to xAOD world in the future
      if (msos->detType()==1) { // its a pixel 
        if (msos->type()==6) { continue; } // not a hole

        //Get pixel cluster then match to clusters in container
        const xAOD::TrackMeasurementValidation* msosClus =  *(msos->trackMeasurementValidationLink());        

        // find this cluster in the shallow copy container
        for (xAOD::TrackMeasurementValidationContainer::iterator clus_itr=(pixClusters)->begin(); clus_itr!=(pixClusters)->end(); ++clus_itr) {
          if ((*clus_itr)->identifier()!=(msosClus)->identifier()) { continue; }
          if ((*clus_itr)->auxdata<float>("charge")!=(msosClus)->auxdata<float>("charge")) { continue; }

          countContributingParticles((*clus_itr));

          // multiple hits on one layer? FIXME - this should fix it
          if (nHitsOnPixelLayer.at((*clus_itr)->auxdata<int>("layer"))<2) {
            nPartPix.at((*clus_itr)->auxdata<int>("layer"))=(*clus_itr)->auxdata<int>("nParticles");
            nPriPartPix.at((*clus_itr)->auxdata<int>("layer"))=(*clus_itr)->auxdata<int>("nPrimaryParticles");
          }

          int detailPixelLayer = (*clus_itr)->auxdata<int>("layer");
          if ((*clus_itr)->auxdataConst<int>("bec")!=0) { detailPixelLayer += 4; }  // endcap +4

          switch (detailPixelLayer) {
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
            case 4 :
              nHitsOnPixelLayer[4]++;
              break;	    
            case 5 :
              nHitsOnPixelLayer[5]++;
              break;	    
            case 6 :
              nHitsOnPixelLayer[6]++;
              break;	    
            default:
              std::cout << "Could not find pixel layer" << std::endl;
              break;
          }

          if ((*clus_itr)->auxdata<int>("layer")==0 && (*clus_itr)->auxdata<int>("bec")==0 && (*trk_itr)->auxdata<float>("IBLCharge")==-1) {
            (*trk_itr)->auxdata<float>("IBLCharge")   = (*clus_itr)->auxdata<float>("charge");
            (*trk_itr)->auxdata<int>("IBLHitSize")    = (*clus_itr)->auxdata<int>("size");
            (*trk_itr)->auxdata<int>("IBLHitSizeZ")   = (*clus_itr)->auxdata<int>("sizeZ");
            (*trk_itr)->auxdata<int>("IBLHitSizePhi") = (*clus_itr)->auxdata<int>("sizePhi");
            if ((*clus_itr)->auxdata<int>("nPrimaryParticles")>1) {
              (*trk_itr)->auxdata<int>("nIBLMergedHits") = (*trk_itr)->auxdata<int>("nIBLMergedHits")+1;
            }
            (*trk_itr)->auxdata<float>("IBLlocalPhi")    = msos->localPhi();
            (*trk_itr)->auxdata<float>("IBLlocalTheta")  = msos->localTheta();

            (*trk_itr)->auxdata<float>("IBLbiasedResidualX")   = msos->biasedResidualX();
            (*trk_itr)->auxdata<float>("IBLbiasedResidualY")   = msos->biasedResidualY();
            (*trk_itr)->auxdata<float>("IBLunbiasedResidualX") = msos->unbiasedResidualX();
            (*trk_itr)->auxdata<float>("IBLunbiasedResidualY") = msos->unbiasedResidualY();
//            (*trk_itr)->auxdata<float>("IBLResidualX")   = (*clus_itr)->localX()-(msos->localX()-msos->biasedResidualX()+msos->unbiasedResidualX());
//            (*trk_itr)->auxdata<float>("IBLResidualY")   = (*clus_itr)->localY()-(msos->localY()-msos->biasedResidualY()+msos->unbiasedResidualY());
            (*trk_itr)->auxdata<float>("IBLResidualX")   = (*clus_itr)->localX()-msos->localX();
            (*trk_itr)->auxdata<float>("IBLResidualY")   = (*clus_itr)->localY()-msos->localY();

          }

          if ((*clus_itr)->auxdata<int>("layer")==1 && (*clus_itr)->auxdata<int>("bec")==0 && (*trk_itr)->auxdata<float>("BLayerCharge")==-1) {
            (*trk_itr)->auxdata<float>("BLayerCharge")   = (*clus_itr)->auxdata<float>("charge");
            (*trk_itr)->auxdata<int>("BLayerHitSize")    = (*clus_itr)->auxdata<int>("size");
            (*trk_itr)->auxdata<int>("BLayerHitSizeZ")   = (*clus_itr)->auxdata<int>("sizeZ");
            (*trk_itr)->auxdata<int>("BLayerHitSizePhi") = (*clus_itr)->auxdata<int>("sizePhi");
            (*trk_itr)->auxdata<float>("BLayerlocalPhi")    = msos->localPhi();
            (*trk_itr)->auxdata<float>("BLayerlocalTheta")  = msos->localTheta();

            (*trk_itr)->auxdata<float>("BLayerbiasedResidualX")   = msos->biasedResidualX();
            (*trk_itr)->auxdata<float>("BLayerbiasedResidualY")   = msos->biasedResidualY();
            (*trk_itr)->auxdata<float>("BLayerunbiasedResidualX") = msos->unbiasedResidualX();
            (*trk_itr)->auxdata<float>("BLayerunbiasedResidualY") = msos->unbiasedResidualY();
//            (*trk_itr)->auxdata<float>("BLayerResidualX")   = (*clus_itr)->localX()-(msos->localX()-msos->biasedResidualX()+msos->unbiasedResidualX());
//            (*trk_itr)->auxdata<float>("BLayerResidualY")   = (*clus_itr)->localY()-(msos->localY()-msos->biasedResidualY()+msos->unbiasedResidualY());
            (*trk_itr)->auxdata<float>("BLayerResidualX")   = (*clus_itr)->localX()-msos->localX();
            (*trk_itr)->auxdata<float>("BLayerResidualY")   = (*clus_itr)->localY()-msos->localY();

          }
 
          break; // once is enough! Move to next MSOS...
        } // loop over cluster container looking for a match
      } // pixels
      else if(msos->detType()==2) { // SCT
        if (msos->type()==6) { continue; } // not a hole
      
        //Get sct cluster then match to clusters in container
        const xAOD::TrackMeasurementValidation* msosClus = *(msos->trackMeasurementValidationLink());
	
        // find this cluster in the shallow copy container
        for (xAOD::TrackMeasurementValidationContainer::iterator sct_itr=(sctClusters)->begin(); sct_itr!=(sctClusters)->end(); ++sct_itr) {
          if ((*sct_itr)->identifier()!=(msosClus)->identifier()) { continue; }
          //if( (*sct_itr)->auxdata< float >( "charge" ) != (msosClus)->auxdata< float >( "charge" ) ) { continue; }

          // what is layer number of 1st SCT layer ? 3 for Run I 4 for Run II...
          int sct_layer = (*sct_itr)->auxdata<int>("layer")/* - 3*/;

// ****  how is sct_layer indexed?

          countContributingParticles((*sct_itr));

          // consider only first measurement on layer
          if (nHitsOnSCTLayer.at(sct_layer)<2) {	  
            nPartSCT.at(sct_layer)    = (*sct_itr)->auxdata<int>("nParticles");
            nPriPartSCT.at(sct_layer) = (*sct_itr)->auxdata<int>("nPrimaryParticles");
          }
  
          // for the track under exam, count the number of hits on each layer - NB: 4 (double) SCT layers!
          switch((*sct_itr)->auxdata<int>("layer")) {
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
          break; // once is enough! Move to next MSOS...
        } // sct
      } // pixel or sct msos  
    } // track's MSOS

    (*trk_itr)->auxdata<std::vector<int>>("pixClusPartContent")    = nPartPix;
    (*trk_itr)->auxdata<std::vector<int>>("pixClusPriPartContent") = nPriPartPix;
    (*trk_itr)->auxdata<std::vector<int>>("sctClusPartContent")    = nPartSCT;
    (*trk_itr)->auxdata<std::vector<int>>("sctClusPriPartContent") = nPriPartSCT;

    (*trk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer") = nHitsOnPixelLayer;
    (*trk_itr)->auxdata<std::vector<int>>("nHitsOnSCTLayer")   = nHitsOnSCTLayer;

  } // loop over tracks
} // dressLayerHitContent

/*
const xAOD::TrackParticle* TwoTrkEffStudy::getHighestProbTrack( std::multimap<int, const xAOD::TrackParticle* >& bcodeMap, int bcode ) {

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

std::vector< const xAOD::TrackMeasurementValidation* > TwoTrkEffStudy :: getAllPixelClusters( const std::vector<ElementLink< xAOD::TrackStateValidationContainer > >& mot ) {
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
*/

//=============
bool TwoTrkEffStudy::checkQuality(const xAOD::TrackParticle* trk, int outerLayer) {
  // start with loose track quality cuts. (see https://indico.cern.ch/event/368819/contribution/4/material/slides/0.pdf )
  uint8_t getInt(0);
  int nSiHits = xAOD::TrackHelper::numberOfSiHits(trk);
  trk->summaryValue(getInt,xAOD::numberOfPixelSharedHits); 
  int nPixSharedHits = getInt;
  trk->summaryValue(getInt,xAOD::numberOfSCTSharedHits); 
  int nSCTSharedHits = getInt;
  trk->summaryValue(getInt,xAOD::numberOfPixelHoles); 
  int nPixHoles = getInt;
  trk->summaryValue(getInt,xAOD::numberOfSCTHoles); 
  int nSCTHoles = getInt;

  if (nSiHits>=7 && nPixSharedHits+nSCTSharedHits<=1 && nPixHoles+nSCTHoles<=2 && nPixHoles<=1) {
    // selection - require at least 1 hit in just outer pixel layer.
    if (trk->isAvailable< std::vector<int> >("nHitsOnPixelLayer")) {
      if (trk->auxdata<std::vector<int>>("nHitsOnPixelLayer")[outerLayer]>0) {
        return true;
      }
    }
  }
  return false;
}

void TwoTrkEffStudy::findPixelClusterOnLayer(const xAOD::TrackParticle* trk, 
                                             xAOD::TrackMeasurementValidationContainer* pixClusters, 
                                             int targetLayer, 
                                             int *nPenetrateLayer, 
                                             const xAOD::TrackStateValidation *msosPenetrate[4], 
                                             int *nClusterLayer, 
                                             const xAOD::TrackMeasurementValidation *pixPenetrateClus[4]) {
    
  int nMOSLayer = 0;
  int nLayer = 0;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > > MeasurementsOnTrack;
  typedef std::vector<ElementLink< xAOD::TrackStateValidationContainer > >::const_iterator MeasurementsOnTrackIter;
//  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack("IDDET1_msosLink");
//  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack("IDTIDE1msosLink");
  static SG::AuxElement::ConstAccessor< MeasurementsOnTrack >  acc_MeasurementsOnTrack("msosLink");

//  if (trk->isAvailable<MeasurementsOnTrack>("IDDET1_msosLink")) {  
//  if (trk->isAvailable<MeasurementsOnTrack>("IDTIDE1msosLink")) {  
  if (trk->isAvailable<MeasurementsOnTrack>("msosLink")) {  
    const MeasurementsOnTrack& measurementsOnTrack = acc_MeasurementsOnTrack(*trk);

    for (MeasurementsOnTrackIter msos_iter = measurementsOnTrack.begin(); msos_iter != measurementsOnTrack.end(); ++msos_iter) {
      if (!(*msos_iter).isValid()) { continue; }	 

      const  xAOD::TrackStateValidation *msos = *(*msos_iter);
      if (msos->detType()!=1) { continue; }
      if (!msos->trackMeasurementValidationLink().isValid()) { continue; }

      // check if this is target layer 
      if ((*(msos->trackMeasurementValidationLink()))->auxdata<int>("layer")!=targetLayer)    { continue; }
      if (targetLayer<4 && (*(msos->trackMeasurementValidationLink()))->auxdataConst<int>("bec")!=0) { continue; }
      msosPenetrate[nMOSLayer] = msos;
      nMOSLayer++;

      // search if cluster exists
      if (!(msos->type()==0)) { continue; }  // confirm this is a measurement type
      bool skip(false);
      for (MeasurementsOnTrackIter  msos_newiter = measurementsOnTrack.begin(); msos_newiter != msos_iter; ++msos_newiter) {
        const  xAOD::TrackStateValidation *this_msos = *(*msos_newiter);
        if ( xAOD::TrackHelper::areSiblings( msos, this_msos) ) { skip=true; break; }
      }
      if (skip) { continue; }

      xAOD::TrackMeasurementValidationContainer::iterator pixclus_itr = pixClusters->begin();
      xAOD::TrackMeasurementValidationContainer::iterator pixclus_end = pixClusters->end();
      xAOD::TrackMeasurementValidation* pixprd = 0;
      for (; pixclus_itr!=pixclus_end; ++pixclus_itr) {
        if ((*pixclus_itr)->identifier()!=(*(msos->trackMeasurementValidationLink()))->identifier()){ continue; }
        if ((*pixclus_itr)->auxdata<float>("charge")!=(*(msos->trackMeasurementValidationLink()))->auxdata<float>("charge")) { continue; }
        pixprd = (*pixclus_itr);
        break;
      }
      pixPenetrateClus[nLayer] = pixprd;
      nLayer++;
    }
  }
  *nPenetrateLayer = nMOSLayer;
  *nClusterLayer = nLayer;
}

//==============================================================================
// Cluster classification 
// ( https://indico.cern.ch/event/377016/contribution/3/material/slides/0.pdf )
//==============================================================================
void TwoTrkEffStudy::clusterPatternClassification(int nRefMSOSLayer, const xAOD::TrackStateValidation *msosRef[4], int nRefClusLayer, const xAOD::TrackMeasurementValidation *pixRefClus[4],
                                                  int nTarMSOSLayer, const xAOD::TrackStateValidation *msosTar[4], int nTarClusLayer, const xAOD::TrackMeasurementValidation *pixTarClus[4],
                                                  int *xclusPatternIndex, double xminMSOSLocalDistance[2], double xminClusLocalDistance[2], double xminClusGlobalDistance[3]) {

  int isSplit1 = 0;
  int isSplit2 = 0;
  int isSplit3 = 0;
  unsigned long id1,id2;
  int clusPatternIndex = -1;
  double minMSOSLocalDistance[2]  = {-1.0,-1.0};
  double minClusLocalDistance[2]  = {-1.0,-1.0};
  double minClusGlobalDistance[3] = {-1.0,-1.0,-1.0};

  if (nRefClusLayer+nTarClusLayer==0) {
    if (nRefMSOSLayer==1 && nTarMSOSLayer==1) {
      if (msosRef[0]->detElementId()==msosTar[0]->detElementId()) { clusPatternIndex = 1; }
      else                                                        { clusPatternIndex = 2; }
    }
    else if ((nRefMSOSLayer==2 && nTarMSOSLayer==1) || (nRefMSOSLayer==1 && nTarMSOSLayer==2)) {
      clusPatternIndex = 3;
    }
    else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) {
      clusPatternIndex = 4;
    }
  }
  else if (nRefClusLayer+nTarClusLayer==1) {
    if (nRefClusLayer==1) { isSplit1 = pixRefClus[0]->auxdata<int>("isSplit"); }
    else                  { isSplit1 = pixTarClus[0]->auxdata<int>("isSplit"); }

    if (isSplit1==0) {
      if (nRefMSOSLayer==1 && nTarMSOSLayer==1) {
        if (msosRef[0]->detElementId()==msosTar[0]->detElementId()) { clusPatternIndex = 5; }
        else                                                        { clusPatternIndex = 6; }
      }
      else if ((nRefMSOSLayer==2 && nTarMSOSLayer==1) || (nRefMSOSLayer==1 && nTarMSOSLayer==2)) {
        if ((nRefMSOSLayer==1 && nRefClusLayer==1) || (nTarMSOSLayer==1 && nTarClusLayer==1)) { clusPatternIndex = 7; }
        else                                                                                  { clusPatternIndex = 8; }
      }
      else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) {
        clusPatternIndex = 9;
      }
    }
    else {
      if (nRefMSOSLayer==1 && nTarMSOSLayer==1) {
        if (msosRef[0]->detElementId()==msosTar[0]->detElementId()) { clusPatternIndex = 10; }
        else                                                        { clusPatternIndex = 11; }
      }
      else if ((nRefMSOSLayer==2 && nTarMSOSLayer==1) || (nRefMSOSLayer==1 && nTarMSOSLayer==2)) {
        if ((nRefMSOSLayer==1 && nRefClusLayer==1) || (nTarMSOSLayer==1 && nTarClusLayer==1)) { clusPatternIndex = 12; }
        else                                                                                  { clusPatternIndex = 13; }
      }
      else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) {
        clusPatternIndex = 14;
      }
    }
  }
  else if (nRefClusLayer+nTarClusLayer==2) {
    if (nRefClusLayer==0 || nTarClusLayer==0) { 
      if (nTarClusLayer==0) { 
        isSplit1 = pixRefClus[0]->auxdata<int>("isSplit"); 
        isSplit2 = pixRefClus[1]->auxdata<int>("isSplit"); 
      }
      if (nRefClusLayer==0) { 
        isSplit1 = pixTarClus[0]->auxdata<int>("isSplit"); 
        isSplit2 = pixTarClus[1]->auxdata<int>("isSplit"); 
      }

      if (isSplit1==0 && isSplit2==0) {
        if      (nRefMSOSLayer==1 || nTarMSOSLayer==1) { clusPatternIndex = 17; }
        else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) { clusPatternIndex = 20; }
      }
      else if (isSplit1==1 && isSplit2==1) {
        if      (nRefMSOSLayer==1 || nTarMSOSLayer==1) { clusPatternIndex = 25; }
        else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) { clusPatternIndex = 28; }
      }
      else {
        if      (nRefMSOSLayer==1 || nTarMSOSLayer==1) { clusPatternIndex = 33; }
        else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) { clusPatternIndex = 36; }
      }
    }
    else {
      double unbiasLocalMSOS_dX = (msosTar[0]->localX()-msosTar[0]->biasedResidualX()+msosTar[0]->unbiasedResidualX())-(msosRef[0]->localX()-msosRef[0]->biasedResidualX()+msosRef[0]->unbiasedResidualX());
      double unbiasLocalMSOS_dY = (msosTar[0]->localY()-msosTar[0]->biasedResidualY()+msosTar[0]->unbiasedResidualY())-(msosRef[0]->localY()-msosRef[0]->biasedResidualY()+msosRef[0]->unbiasedResidualY());
      double clusLocal_dX = pixTarClus[0]->localX()-pixRefClus[0]->localX();
      double clusLocal_dY = pixTarClus[0]->localY()-pixRefClus[0]->localY();
      double clusGlobal_dX = pixTarClus[0]->globalX()-pixRefClus[0]->globalX();
      double clusGlobal_dY = pixTarClus[0]->globalY()-pixRefClus[0]->globalY();
      double clusGlobal_dZ = pixTarClus[0]->globalZ()-pixRefClus[0]->globalZ();

//      minMSOSLocalDistance  = TMath::Sqrt(unbiasLocalMSOS_dX*unbiasLocalMSOS_dX+unbiasLocalMSOS_dY*unbiasLocalMSOS_dY);
//      minClusLocalDistance  = TMath::Sqrt(clusLocal_dX*clusLocal_dX+clusLocal_dY*clusLocal_dY);
//      minClusGlobalDistance = TMath::Sqrt(clusGlobal_dX*clusGlobal_dX+clusGlobal_dY*clusGlobal_dY+clusGlobal_dZ*clusGlobal_dZ);

      minMSOSLocalDistance[0]  = unbiasLocalMSOS_dX;
      minClusLocalDistance[0]  = clusLocal_dX;
      minClusGlobalDistance[0] = clusGlobal_dX;

      minMSOSLocalDistance[1]  = unbiasLocalMSOS_dY;
      minClusLocalDistance[1]  = clusLocal_dY;
      minClusGlobalDistance[1] = clusGlobal_dY;

      minClusGlobalDistance[2] = clusGlobal_dZ;

      isSplit1 = pixRefClus[0]->auxdata<int>("isSplit");
      isSplit2 = pixTarClus[0]->auxdata<int>("isSplit");
      id1 = pixRefClus[0]->auxdata<unsigned long>("identifier");
      id2 = pixTarClus[0]->auxdata<unsigned long>("identifier");

      if (nRefMSOSLayer==1 && nTarMSOSLayer==1) {
        if (msosRef[0]->detElementId()==msosTar[0]->detElementId()) { 
          if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 15; }
          else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 39; }
          else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 23; }
          else                                             { clusPatternIndex = 31; }
        }
        else { 
          if      (isSplit1==0 && isSplit2==0) { clusPatternIndex = 16; }
          else if (isSplit1==1 && isSplit2==1) { clusPatternIndex = 24; }
          else                                 { clusPatternIndex = 32; }
        }
      }
      else if ((nRefMSOSLayer==2 && nTarMSOSLayer==1) || (nRefMSOSLayer==1 && nTarMSOSLayer==2)) {
        if (buildPixelModuleID(pixRefClus[0])==buildPixelModuleID(pixTarClus[0])) { 
          if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 18; }
          else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 40; }
          else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 26; }
          else                                             { clusPatternIndex = 34; }
        }
        else { 
          if      (isSplit1==0 && isSplit2==0) { clusPatternIndex = 19; }
          else if (isSplit1==1 && isSplit2==1) { clusPatternIndex = 27; }
          else                                 { clusPatternIndex = 35; }
        }
      }
      else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) {
        if (buildPixelModuleID(pixRefClus[0])==buildPixelModuleID(pixTarClus[0])) { 
          if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 21; }
          else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 41; }
          else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 29; }
          else                                             { clusPatternIndex = 37; }
        }
        else { 
          if      (isSplit1==0 && isSplit2==0) { clusPatternIndex = 22; }
          else if (isSplit1==1 && isSplit2==1) { clusPatternIndex = 30; }
          else                                 { clusPatternIndex = 38; }
        }
      }
    }
  }
  else if (nRefClusLayer+nTarClusLayer==3) {

    if (nRefClusLayer==0 || nTarClusLayer==0) {
      std::cout << "Strange... three clusters per layer in one track.." << std::endl;
    }

    int iref1 = 0; int iref2 = 1;
    int itar1 = 0; int itar2 = 1;
    if (buildPixelModuleID(pixRefClus[0])!=buildPixelModuleID(pixTarClus[0])) {
      if (nRefClusLayer==2) { iref1=1; itar1=0; iref2=0; itar2=1; }
      else                  { iref1=0; itar1=1; iref2=1; itar2=0; }
    }

    double unbiasLocalMSOS_dX = (msosTar[itar1]->localX()-msosTar[itar1]->biasedResidualX()+msosTar[itar1]->unbiasedResidualX())-(msosRef[iref1]->localX()-msosRef[iref1]->biasedResidualX()+msosRef[iref1]->unbiasedResidualX());
    double unbiasLocalMSOS_dY = (msosTar[itar1]->localY()-msosTar[itar1]->biasedResidualY()+msosTar[itar1]->unbiasedResidualY())-(msosRef[iref1]->localY()-msosRef[iref1]->biasedResidualY()+msosRef[iref1]->unbiasedResidualY());
    double clusLocal_dX = pixTarClus[itar1]->localX()-pixRefClus[iref1]->localX();
    double clusLocal_dY = pixTarClus[itar1]->localY()-pixRefClus[iref1]->localY();
    double clusGlobal_dX = pixTarClus[itar1]->globalX()-pixRefClus[iref1]->globalX();
    double clusGlobal_dY = pixTarClus[itar1]->globalY()-pixRefClus[iref1]->globalY();
    double clusGlobal_dZ = pixTarClus[itar1]->globalZ()-pixRefClus[iref1]->globalZ();

//    minMSOSLocalDistance  = TMath::Sqrt(unbiasLocalMSOS_dX*unbiasLocalMSOS_dX+unbiasLocalMSOS_dY*unbiasLocalMSOS_dY);
//    minClusLocalDistance  = TMath::Sqrt(clusLocal_dX*clusLocal_dX+clusLocal_dY*clusLocal_dY);
//    minClusGlobalDistance = TMath::Sqrt(clusGlobal_dX*clusGlobal_dX+clusGlobal_dY*clusGlobal_dY+clusGlobal_dZ*clusGlobal_dZ);

    minMSOSLocalDistance[0]  = unbiasLocalMSOS_dX;
    minClusLocalDistance[0]  = clusLocal_dX;
    minClusGlobalDistance[0] = clusGlobal_dX;

    minMSOSLocalDistance[1]  = unbiasLocalMSOS_dY;
    minClusLocalDistance[1]  = clusLocal_dY;
    minClusGlobalDistance[1] = clusGlobal_dY;

    minClusGlobalDistance[2] = clusGlobal_dZ;


    isSplit1 = pixRefClus[iref1]->auxdata<int>("isSplit");
    isSplit2 = pixTarClus[itar1]->auxdata<int>("isSplit");
    id1 = pixRefClus[iref1]->auxdata<unsigned long>("identifier");
    id2 = pixTarClus[itar1]->auxdata<unsigned long>("identifier");

    if (nRefClusLayer==2) { isSplit3 = pixRefClus[iref2]->auxdata<int>("isSplit"); }
    else                  { isSplit3 = pixTarClus[itar2]->auxdata<int>("isSplit"); }

    if ((nRefMSOSLayer==2 && nTarMSOSLayer==1) || (nRefMSOSLayer==1 && nTarMSOSLayer==2)) {
      if (isSplit3==0) {
        if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 42; }
        else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 54; }
        else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 51; }
        else                                             { clusPatternIndex = 46; }
      }
      else {
        if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 47; }
        else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 56; }
        else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 44; }
        else                                             { clusPatternIndex = 50; }
      }
    }
    else if (nRefMSOSLayer==2 && nTarMSOSLayer==2) {
      if (isSplit3==0) {
        if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 43; }
        else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 55; }
        else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 53; }
        else                                             { clusPatternIndex = 48; }
      }
      else {
        if      (isSplit1==0 && isSplit2==0)             { clusPatternIndex = 49; }
        else if (isSplit1==1 && isSplit2==1 && id1==id2) { clusPatternIndex = 57; }
        else if (isSplit1==1 && isSplit2==1 && id1!=id2) { clusPatternIndex = 45; }
        else                                             { clusPatternIndex = 52; }
      }
    }
  }
  else if (nRefClusLayer+nTarClusLayer==4) {

    if (nRefClusLayer<2 || nTarClusLayer<2) {
      std::cout << "Strange... three or four clusters per layer in one track.." << std::endl;
    }

    // make pair (1a,2a), (1b,2b)
    int iref1a = 0; int iref2a = 0;
    int iref1b = 1; int iref2b = 1;
    if (buildPixelModuleID(pixRefClus[0])!=buildPixelModuleID(pixTarClus[0])) { iref1a=0; iref2a=1; iref1b=1; iref2b=0; }

    double unbiasLocalMSOS_dX1 = (msosTar[iref2a]->localX()-msosTar[iref2a]->biasedResidualX()+msosTar[iref2a]->unbiasedResidualX())-(msosRef[iref1a]->localX()-msosRef[iref1a]->biasedResidualX()+msosRef[iref1a]->unbiasedResidualX());
    double unbiasLocalMSOS_dY1 = (msosTar[iref2a]->localY()-msosTar[iref2a]->biasedResidualY()+msosTar[iref2a]->unbiasedResidualY())-(msosRef[iref1a]->localY()-msosRef[iref1a]->biasedResidualY()+msosRef[iref1a]->unbiasedResidualY());
    double clusLocal_dX1 = pixTarClus[iref2a]->localX()-pixRefClus[iref1a]->localX();
    double clusLocal_dY1 = pixTarClus[iref2a]->localY()-pixRefClus[iref1a]->localY();
    double clusGlobal_dX1 = pixTarClus[iref2a]->globalX()-pixRefClus[iref1a]->globalX();
    double clusGlobal_dY1 = pixTarClus[iref2a]->globalY()-pixRefClus[iref1a]->globalY();
    double clusGlobal_dZ1 = pixTarClus[iref2a]->globalZ()-pixRefClus[iref1a]->globalZ();

    double unbiasLocalMSOS_dX2 = (msosTar[iref2b]->localX()-msosTar[iref2b]->biasedResidualX()+msosTar[iref2b]->unbiasedResidualX())-(msosRef[iref1b]->localX()-msosRef[iref1b]->biasedResidualX()+msosRef[iref1b]->unbiasedResidualX());
    double unbiasLocalMSOS_dY2 = (msosTar[iref2b]->localY()-msosTar[iref2b]->biasedResidualY()+msosTar[iref2b]->unbiasedResidualY())-(msosRef[iref1b]->localY()-msosRef[iref1b]->biasedResidualY()+msosRef[iref1b]->unbiasedResidualY());
    double clusLocal_dX2 = pixTarClus[iref2b]->localX()-pixRefClus[iref1b]->localX();
    double clusLocal_dY2 = pixTarClus[iref2b]->localY()-pixRefClus[iref1b]->localY();
    double clusGlobal_dX2 = pixTarClus[iref2b]->globalX()-pixRefClus[iref1b]->globalX();
    double clusGlobal_dY2 = pixTarClus[iref2b]->globalY()-pixRefClus[iref1b]->globalY();
    double clusGlobal_dZ2 = pixTarClus[iref2b]->globalZ()-pixRefClus[iref1b]->globalZ();

//    minMSOSLocalDistance  = TMath::Max(TMath::Sqrt(unbiasLocalMSOS_dX1*unbiasLocalMSOS_dX1+unbiasLocalMSOS_dY1*unbiasLocalMSOS_dY1), 
//                                       TMath::Sqrt(unbiasLocalMSOS_dX2*unbiasLocalMSOS_dX2+unbiasLocalMSOS_dY2*unbiasLocalMSOS_dY2));
//    minClusLocalDistance  = TMath::Max(TMath::Sqrt(clusLocal_dX1*clusLocal_dX1+clusLocal_dY1*clusLocal_dY1), 
//                                       TMath::Sqrt(clusLocal_dX2*clusLocal_dX2+clusLocal_dY2*clusLocal_dY2));
//    minClusGlobalDistance = TMath::Max(TMath::Sqrt(clusGlobal_dX1*clusGlobal_dX1+clusGlobal_dY1*clusGlobal_dY1+clusGlobal_dZ1*clusGlobal_dZ1), 
//                                       TMath::Sqrt(clusGlobal_dX2*clusGlobal_dX2+clusGlobal_dY2*clusGlobal_dY2+clusGlobal_dZ2*clusGlobal_dZ2));

    if (TMath::Sqrt(unbiasLocalMSOS_dX1*unbiasLocalMSOS_dX1+unbiasLocalMSOS_dY1*unbiasLocalMSOS_dY1)>TMath::Sqrt(unbiasLocalMSOS_dX2*unbiasLocalMSOS_dX2+unbiasLocalMSOS_dY2*unbiasLocalMSOS_dY2)) {
      minMSOSLocalDistance[0] = unbiasLocalMSOS_dX1;
      minMSOSLocalDistance[1] = unbiasLocalMSOS_dY1;
    }
    else {
      minMSOSLocalDistance[0] = unbiasLocalMSOS_dX2;
      minMSOSLocalDistance[1] = unbiasLocalMSOS_dY2;
    }

    if (TMath::Sqrt(clusLocal_dX1*clusLocal_dX1+clusLocal_dY1*clusLocal_dY1)>TMath::Sqrt(clusLocal_dX2*clusLocal_dX2+clusLocal_dY2*clusLocal_dY2)) {
      minClusLocalDistance[0] = clusLocal_dX1;
      minClusLocalDistance[1] = clusLocal_dY1;
    }
    else {
      minClusLocalDistance[0] = clusLocal_dX2;
      minClusLocalDistance[1] = clusLocal_dY2;
    }

    if (TMath::Sqrt(clusGlobal_dX1*clusGlobal_dX1+clusGlobal_dY1*clusGlobal_dY1+clusGlobal_dZ1*clusGlobal_dZ1)>TMath::Sqrt(clusGlobal_dX2*clusGlobal_dX2+clusGlobal_dY2*clusGlobal_dY2+clusGlobal_dZ2*clusGlobal_dZ2)) {
      minClusGlobalDistance[0] = clusGlobal_dX1;
      minClusGlobalDistance[1] = clusGlobal_dY1;
      minClusGlobalDistance[2] = clusGlobal_dZ1;
    }
    else {
      minClusGlobalDistance[0] = clusGlobal_dX2;
      minClusGlobalDistance[1] = clusGlobal_dY2;
      minClusGlobalDistance[2] = clusGlobal_dZ2;
    }

    int isSplit1a = pixRefClus[iref1a]->auxdata<int>("isSplit");
    int isSplit2a = pixTarClus[iref2a]->auxdata<int>("isSplit");
    unsigned long id1a = pixRefClus[iref1a]->auxdata<unsigned long>("identifier");
    unsigned long id2a = pixTarClus[iref2a]->auxdata<unsigned long>("identifier");

    int isSplit1b = pixRefClus[iref1b]->auxdata<int>("isSplit");
    int isSplit2b = pixTarClus[iref2b]->auxdata<int>("isSplit");
    unsigned long id1b = pixRefClus[iref1b]->auxdata<unsigned long>("identifier");
    unsigned long id2b = pixTarClus[iref2b]->auxdata<unsigned long>("identifier");

    if (isSplit1a==0 && isSplit2a==0) {
      if      (isSplit1b==0 && isSplit2b==0)               { clusPatternIndex = 58; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b==id2b) { clusPatternIndex = 64; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b!=id2b) { clusPatternIndex = 62; }
      else                                                 { clusPatternIndex = 63; }
    }
    else if (isSplit1a==1 && isSplit2a==1 && id1a==id2a) {
      if      (isSplit1b==0 && isSplit2b==0)               { clusPatternIndex = 64; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b==id2b) { clusPatternIndex = 67; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b!=id2b) { clusPatternIndex = 65; }
      else                                                 { clusPatternIndex = 66; }
    }
    else if (isSplit1a==1 && isSplit2a==1 && id1a!=id2a) {
      if      (isSplit1b==0 && isSplit2b==0)               { clusPatternIndex = 62; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b==id2b) { clusPatternIndex = 65; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b!=id2b) { clusPatternIndex = 59; }
      else                                                 { clusPatternIndex = 60; }
    }
    else {
      if      (isSplit1b==0 && isSplit2b==0)               { clusPatternIndex = 63; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b==id2b) { clusPatternIndex = 66; }
      else if (isSplit1b==1 && isSplit2b==1 && id1b!=id2b) { clusPatternIndex = 60; }
      else                                                 { clusPatternIndex = 61; }
    }
  }


  if (clusPatternIndex==-1) {
    std::cout << "clusters are not assigned : (nRefMSOS, nRefClus, nTarMSOS, nTarClus) = ( " << nRefMSOSLayer << ", " << nRefClusLayer << ", " << nTarMSOSLayer << ", " << nTarClusLayer << " )" << std::endl;
  }


  (*xclusPatternIndex)      = clusPatternIndex;
//  (*xminMSOSLocalDistance)  = minMSOSLocalDistance;
//  (*xminClusLocalDistance)  = minClusLocalDistance;
//  (*xminClusGlobalDistance) = minClusGlobalDistance;

  xminMSOSLocalDistance[0]  = minMSOSLocalDistance[0];
  xminClusLocalDistance[0]  = minClusLocalDistance[0];
  xminClusGlobalDistance[0] = minClusGlobalDistance[0];

  xminMSOSLocalDistance[1]  = minMSOSLocalDistance[1];
  xminClusLocalDistance[1]  = minClusLocalDistance[1];
  xminClusGlobalDistance[1] = minClusGlobalDistance[1];

  xminClusGlobalDistance[2] = minClusGlobalDistance[2];

}

