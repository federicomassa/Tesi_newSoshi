#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <ParticleAnalysis/ClusterAnalysis.h>
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
ClassImp(ClusterAnalysis)

// switches - made static because need to be set (in job steering macro) before calling constructor
char* ClusterAnalysis::pileupWeightFile = "";

ClusterAnalysis::ClusterAnalysis() {

  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

}


EL::StatusCode ClusterAnalysis::setupJob(EL::Job& job) {
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

EL::StatusCode ClusterAnalysis::histInitialize() {
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
  trkHist_reco_all            = new TrackHists("reco_all");
  trkHist_reco_selected    = new TrackHists("reco_selected");

  clusHist_reco_selected    = new ClusterHists("reco_selected");

  // event histograms
  eventHist_all               = new EventHists("all");

  trkHist_reco_all            -> Init( wk() );
  trkHist_reco_selected    -> Init( wk() );

  clusHist_reco_selected    -> Init( wk() );

  eventHist_all               ->Init( wk() );

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ClusterAnalysis::fileExecute() {
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ClusterAnalysis::changeInput (bool /* firstFile */) {
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  TFile* infile = wk()->inputFile ();
  std::cout << infile->GetName() << std::endl;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ClusterAnalysis::initialize() {
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
  const char* fullGRLFilePath = gSystem->ExpandPathName("data15_13TeV.periodAllYear_DetStatus-v64-pro19_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml");
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

EL::StatusCode ClusterAnalysis::execute() {
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
//      std::cout << "Look in here to change -> ClusterAnalysis :: execute()" << std::endl;
//    }
//    return EL::StatusCode::SUCCESS; 
//  }

  EventFeatures eventFeatures;

  //---------------------------
  // Reset barcode map in histogram class
  //---------------------------
  trkHist_reco_all            -> resetBarcodes();
  trkHist_reco_selected    -> resetBarcodes();

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
//  bool flag_jet850 = false;
//  for (xAOD::JetContainer::const_iterator jet_itr=jetAntiKt4EMTopo->begin(); jet_itr!=jetAntiKt4EMTopo->end(); jet_itr++) {
//    if ((*jet_itr)->pt()/1000.0>850.0) { flag_jet850 = true; }
//  }
//  if (flag_jet850) { return EL::StatusCode::SUCCESS; }

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
  uint32_t lumibk = 0.0;

  if (isMC) {
    pileupWeight = eventInfo->auxdata<double>("PileupWeight");
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
    lumibk = eventInfo->lumiBlock();
  }

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
  }


  // Trigger selection
//  if (eventFeatures.trigCounter==0) { return EL::StatusCode::SUCCESS; }


  // weight for histogram
  float histWeight = pileupWeight*xsecWeight*mcEvWeight*triggerWeight;

  //==========================
  // Count number of vertices
  //==========================
  int numVtx = 0;
  for (xAOD::VertexContainer::const_iterator vtx_itr=primaryVertex->begin(); vtx_itr!=primaryVertex->end(); vtx_itr++) {
    if ((*vtx_itr)->vertexType()!=0) { numVtx++; }
  }

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
  // Build reconstructed tracks class
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

  int cut_ntrackIsolation = 100;             // 0 track except target track.
  double cut_TrkPt = 0.0;                    // [GeV]

  int numRecoTracks = 0;
  int numRecoTracks1GeV = 0;
  int nIsoTrk = 0;

  xAOD::TrackParticleContainer* selectedRecoTracks = new xAOD::TrackParticleContainer();
  xAOD::TrackParticleAuxContainer* selectedRecoTracksAux = new xAOD::TrackParticleAuxContainer();
  selectedRecoTracks->setStore(selectedRecoTracksAux); //< Connect the two

  for (xAOD::TrackParticleContainer::const_iterator recoTrk_itr=recoTracks->begin(); recoTrk_itr!=recoTracks->end(); recoTrk_itr++) {

    if (!checkQuality(*recoTrk_itr))           { continue; }    // Quality cut
    numRecoTracks++;

    if ((*recoTrk_itr)->pt()/1000.0<cut_TrkPt) { continue; } // pt cut
    numRecoTracks1GeV++;

    int nonIsoTrk = 0;
    for (xAOD::TrackParticleContainer::const_iterator nonIsoTrk_itr=recoTracks->begin(); nonIsoTrk_itr!=recoTracks->end(); nonIsoTrk_itr++) {
      if ((*nonIsoTrk_itr)->pt()/1000.0<cut_TrkPt) { continue; }   // pt cut
      if (!checkQuality(*nonIsoTrk_itr))           { continue; } // Quality cut
      double dR = deltaR(*recoTrk_itr,*nonIsoTrk_itr);
      if (dR<0.1) { nonIsoTrk++; }
    }

    if (nonIsoTrk<=1+cut_ntrackIsolation) {
      nIsoTrk++;

      // Find cluster on surface
      int nRefMSOSLayer = 0;
      int nRefClusLayer = 0;
      const xAOD::TrackStateValidation *msosRef[4];
      const xAOD::TrackMeasurementValidation *pixRefClus[4];
      findPixelClusterOnLayer(*recoTrk_itr, pixClusters, m_targetLayer, &nRefMSOSLayer, &msosRef[0], &nRefClusLayer, &pixRefClus[0]);

      // Copy this track to the new container:
      xAOD::TrackParticle* tp = new xAOD::TrackParticle();
      tp->makePrivateStore( **recoTrk_itr);

      static SG::AuxElement::Decorator<int> nRefMSOSAtLayerOfInterest("nRefMSOSAtLayerOfInterest");
      static SG::AuxElement::Decorator<int> nRefHitAtLayerOfInterest("nRefHitAtLayerOfInterest");
      static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> RefMSOS1("RefMSOS1");
      static SG::AuxElement::Decorator<const xAOD::TrackStateValidation*> RefMSOS2("RefMSOS2");
      static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelRefClus1("PixelRefClus1");
      static SG::AuxElement::Decorator<const xAOD::TrackMeasurementValidation*> PixelRefClus2("PixelRefClus2");
      nRefMSOSAtLayerOfInterest(*tp) = nRefMSOSLayer;
      nRefHitAtLayerOfInterest(*tp) = nRefClusLayer;
      RefMSOS1(*tp) = msosRef[0];
      RefMSOS2(*tp) = msosRef[1];
      PixelRefClus1(*tp) = pixRefClus[0];
      PixelRefClus2(*tp) = pixRefClus[1];

      static SG::AuxElement::Decorator<int>   TargetLayerOfInterest("TargetLayerOfInterest");
      TargetLayerOfInterest(*tp)        = m_targetLayer;

      static SG::AuxElement::Decorator<float>   z0Corr1("z0Corr1");
      static SG::AuxElement::Decorator<float>   z0Corr2("z0Corr2");
      z0Corr1(*tp) = (*recoTrk_itr)->z0()+(*recoTrk_itr)->vz()-vtxZ;

      static SG::AuxElement::Decorator<float>   muTrk("muTrk");
      static SG::AuxElement::Decorator<int>     lbTrk("lbTrk");
      static SG::AuxElement::Decorator<int>     pixClusTrk("pixClusTrk");
      muTrk(*tp) = muValue;
      lbTrk(*tp) = lumibk;
      pixClusTrk(*tp) = pixClustersOrig->size();

      selectedRecoTracks->push_back( tp );
    }
  }
  //=========================================================================================================================

  // Fill all pixel clusters

//  std::cout << "ok1" << std::endl;

  int numPixelClustersOnTarget = 0;
  for (xAOD::TrackMeasurementValidationContainer::iterator clus_itr=(pixClusters)->begin(); clus_itr!=(pixClusters)->end(); ++clus_itr) {
    int detailPixelLayer = (*clus_itr)->auxdata<int>("layer");
    if ((*clus_itr)->auxdataConst<int>("bec")!=0) { detailPixelLayer += 4; }  // endcap +4
    if (detailPixelLayer==m_targetLayer) {
      numPixelClustersOnTarget++;
    }
  }

//  std::cout << "ok2" << std::endl;

  //---------------------
  // FILL PLOTS!!!
  //---------------------

  //=========================================================================================================================
  // Fill histograms for track-pair for reconstructed tracks
  //=========================================================================================================================
  for (xAOD::TrackParticleContainer::const_iterator recoTrk_itr=selectedRecoTracks->begin(); recoTrk_itr!=selectedRecoTracks->end(); recoTrk_itr++) {

    int TargetLayer  = (*recoTrk_itr)->auxdata<int>("TargetLayerOfInterest");
    int nRefMSOSLayer = (*recoTrk_itr)->auxdata<int>("nRefMSOSAtLayerOfInterest");
    int nRefClusLayer = (*recoTrk_itr)->auxdata<int>("nRefHitAtLayerOfInterest");
    const xAOD::TrackStateValidation *msosRef1 = (*recoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("RefMSOS1");
    const xAOD::TrackStateValidation *msosRef2 = (*recoTrk_itr)->auxdata<const xAOD::TrackStateValidation *>("RefMSOS2");
    const xAOD::TrackMeasurementValidation *pixRefClus1 = (*recoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus1");
    const xAOD::TrackMeasurementValidation *pixRefClus2 = (*recoTrk_itr)->auxdata<const xAOD::TrackMeasurementValidation *>("PixelRefClus2");

    // Now fill track histogram; avoid overlap region
    trkHist_reco_all -> FillHists((*recoTrk_itr),histWeight);
//    if (nRefMSOSLayer==1) { 
    if (nRefMSOSLayer>0) { 
      if ((*recoTrk_itr)->isAvailable< std::vector<int> >("nHitsOnPixelLayer")) {
/*
        if ((*recoTrk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer")[m_targetLayer-1]>0 && (*recoTrk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer")[m_targetLayer+1]>0) {
          trkHist_reco_all -> FillHists((*recoTrk_itr),histWeight);
          if (nRefClusLayer>0) { 
            trkHist_reco_selected  -> FillHists((*recoTrk_itr),histWeight);
            clusHist_reco_selected -> FillHists(pixRefClus1,histWeight); 
          }
        }
        if ((*recoTrk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer")[m_targetLayer-1]>0) {
          trkHist_reco_all -> FillHists((*recoTrk_itr),histWeight);
          if (nRefClusLayer>0) { 
            trkHist_reco_selected  -> FillHists((*recoTrk_itr),histWeight);
            clusHist_reco_selected -> FillHists(pixRefClus1,histWeight); 
          }
        }
        if ((*recoTrk_itr)->auxdata<std::vector<int>>("nHitsOnPixelLayer")[m_targetLayer+1]>0) {
          trkHist_reco_all -> FillHists((*recoTrk_itr),histWeight);
          if (nRefClusLayer>0) { 
            trkHist_reco_selected  -> FillHists((*recoTrk_itr),histWeight);
            clusHist_reco_selected -> FillHists(pixRefClus1,histWeight); 
          }
        }
//        trkHist_reco_all -> FillHists((*recoTrk_itr),histWeight);
*/
        if (nRefClusLayer>0) { 
          trkHist_reco_selected  -> FillHists((*recoTrk_itr),histWeight);
          clusHist_reco_selected -> FillHists(pixRefClus1,histWeight); 
        }

      }
    }

  }

  eventFeatures.muValue     = muValue;
  eventFeatures.muValueCorr = muValueCorr;
  eventFeatures.nRecoTracks = numRecoTracks;
  eventFeatures.nRecoTracks1GeV = numRecoTracks1GeV;
  eventFeatures.nPixelClusters = pixClustersOrig->size();
  eventFeatures.nSCTClusters = sctClustersOrig->size();
  eventFeatures.nSiClusters = pixClustersOrig->size() + sctClustersOrig->size();
  eventFeatures.nIsoTrack = nIsoTrk;
  eventFeatures.nVertex = numVtx;
  eventFeatures.nPixelClustersAtIBL = numPixelClustersOnTarget;

  eventHist_all->FillHists(eventFeatures,histWeight); 

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


EL::StatusCode ClusterAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ClusterAnalysis :: finalize ()
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



EL::StatusCode ClusterAnalysis :: histFinalize ()
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

int ClusterAnalysis::buildPixelModuleID(const xAOD::TrackMeasurementValidation* prd) const {
  int bec        = prd->auxdataConst< int >("bec");   // -4,4
  int layer      = prd->auxdataConst< int >("layer");      // 0,4
  int eta_module = prd->auxdataConst< int >("eta_module"); // -20,20
  int phi_module = prd->auxdataConst< int >("phi_module"); // 0-50
  return bec * 1000000 + layer * 100000 + phi_module * 100 + eta_module + 50; 
}

void ClusterAnalysis::countContributingParticles(xAOD::TrackMeasurementValidation* prd) {
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


void ClusterAnalysis::dressLayerHitContent(xAOD::TrackParticleContainer* trks, xAOD::TrackMeasurementValidationContainer* pixClusters, xAOD::TrackMeasurementValidationContainer* sctClusters) {

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

    (*trk_itr)->auxdata<float>("BLayerCharge")    = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSize")     = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSizeZ")    = -1;
    (*trk_itr)->auxdata<int>("BLayerHitSizePhi")  = -1;
    (*trk_itr)->auxdata<int>("nBLayerMergedHits") =  0;

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

          if ((*clus_itr)->auxdata<int>("layer")==0 && (*trk_itr)->auxdata<float>("BLayerCharge")==-1) {
            (*trk_itr)->auxdata<float>("BLayerCharge")   = (*clus_itr)->auxdata<float>("charge");
            (*trk_itr)->auxdata<int>("BLayerHitSize")    = (*clus_itr)->auxdata<int>("size");
            (*trk_itr)->auxdata<int>("BLayerHitSizeZ")   = (*clus_itr)->auxdata<int>("sizeZ");
            (*trk_itr)->auxdata<int>("BLayerHitSizePhi") = (*clus_itr)->auxdata<int>("sizePhi");
            if ((*clus_itr)->auxdata<int>("nPrimaryParticles")>1) {
              (*trk_itr)->auxdata<int>("nBLayerMergedHits") = (*trk_itr)->auxdata<int>("nBLayerMergedHits")+1;
            }
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

//=============
bool ClusterAnalysis::checkQuality(const xAOD::TrackParticle* trk) {
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
    return true;
  }
  return false;
}

void ClusterAnalysis::findPixelClusterOnLayer(const xAOD::TrackParticle* trk, 
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


