/******************************************
 * Plot of tracks and tracking quantities
 *
 * G. Facini
 * Sep 10 10:03:25 CEST 2014
 *
 ******************************************/

#include <HistManager/TrackHistManager.h>

TrackHistManager::TrackHistManager(TString name, bool doPrimary, bool doSecondary, bool doFake, bool doDups, bool splitOnHits) {
   m_name = name;
   
   m_doPrimary   = doPrimary;
   m_doSecondary = doSecondary;
   m_doFake      = doFake;
   m_doDups      = doDups;
   m_splitOnHits = splitOnHits;
   m_all               = new TrackHists(name+"_all");
   if(m_doPrimary) {
      m_primary        = new TrackHists(name+"_primary");
   }
   if(m_doSecondary) { 
      m_secondary      = new TrackHists(name+"_secondary");
   }
   if(m_doFake) {
      m_fake           = new TrackHists(name+"_fake");
      m_true           = new TrackHists(name+"_true");
   }
   if( (m_doPrimary && m_doSecondary) && m_doDups) { 
     m_primary_dup     = new TrackHists(name+"_primary_dup");
     m_secondary_dup   = new TrackHists(name+"_secondary_dup");
   }
   if( (m_doPrimary) && m_splitOnHits) {
     m_primary_bhit         = new TrackHists(name+"_primary_bhit");
     m_primary_bhitShare    = new TrackHists(name+"_primary_bhitShare");
     m_primary_bhitSplit    = new TrackHists(name+"_primary_bhitSplit");
     m_primary_bhitMerge    = new TrackHists(name+"_primary_bhitMerge");
   }
}

TrackHistManager::~TrackHistManager() {}

void TrackHistManager::Init(EL::Worker* wk) {

  std::cout << "Init TrackHistManager " << m_name << std::endl;

  //Now initialize the worker and add TrackHists objects to vector, needed for RunHistManager

  m_all           -> Init( wk );
  m_histList       . push_back(m_all);
  if(m_doPrimary) {
    m_primary       -> Init( wk );
    m_histList       . push_back(m_primary);
  }
  if(m_doSecondary) { 
    m_secondary     -> Init( wk );
    m_histList       . push_back(m_secondary);
  }
  if(m_doFake) {   
    m_fake          -> Init( wk );
    m_true          -> Init( wk );

    m_histList       . push_back(m_fake);
    m_histList       . push_back(m_true);
  }
  if( (m_doPrimary && m_doSecondary) && m_doDups) {
    m_primary_dup   -> Init( wk );
    m_secondary_dup -> Init( wk );

    m_histList       . push_back(m_primary_dup);
    m_histList       . push_back(m_secondary_dup);
  }
  if( (m_doPrimary) && m_splitOnHits) {
    m_primary_bhit         -> Init( wk );
    m_primary_bhitShare    -> Init( wk );
    m_primary_bhitSplit    -> Init( wk );
    m_primary_bhitMerge    -> Init( wk );

    m_histList              . push_back(m_primary_bhit);
    m_histList              . push_back(m_primary_bhitShare);
    m_histList              . push_back(m_primary_bhitSplit);
    m_histList              . push_back(m_primary_bhitMerge);
  }

} // BookHists

void TrackHistManager::resetBarcodes() {

  m_all             -> resetBarcodes();
  if(m_doPrimary) {
    m_primary       -> resetBarcodes();
  }
  if(m_doSecondary) {   
    m_secondary     -> resetBarcodes();
  }
  if(m_doFake) {     
    m_fake          -> resetBarcodes();
    m_true          -> resetBarcodes();
  }
  if( (m_doPrimary && m_doSecondary) && m_doDups) {
    m_primary_dup   -> resetBarcodes();
    m_secondary_dup -> resetBarcodes();
  }
  if( (m_doPrimary) && m_splitOnHits) {
    m_primary_bhit         -> resetBarcodes();
    m_primary_bhitShare    -> resetBarcodes();
    m_primary_bhitSplit    -> resetBarcodes();
    m_primary_bhitMerge    -> resetBarcodes();
  }

} // BookHists

void TrackHistManager::FillHists(const xAOD::TrackParticle* trk, float weight) const {

  m_all->FillHists( trk,weight);
  if( m_doFake && xAOD::TrackHelper::isFake( trk ) ){
    m_fake->FillHists( trk,weight);
  }
  else if (m_doFake && !xAOD::TrackHelper::isFake( trk ) ) {
    m_true->FillHists( trk,weight);
  }
  else if( m_doPrimary && xAOD::TrackHelper::isPrimary( trk ) ){
    if( !m_primary->hasBeenUsed( trk ) ){
      m_primary->FillHists( trk,weight);
      if(m_splitOnHits) {
        uint8_t getInt(0);   // for accessing summary information
        //has a b-layer hit?
        trk->summaryValue( getInt, xAOD::numberOfBLayerHits );
        if( getInt > 0 ) { 
          m_primary_bhit     ->FillHists( trk,weight);
          // has a shared b-layer hit
          trk->summaryValue( getInt, xAOD::numberOfBLayerSharedHits);
          if( getInt > 0 ) { m_primary_bhitShare->FillHists( trk,weight); }
          // has a split b-layer hit          
          trk->summaryValue( getInt, xAOD::numberOfBLayerSplitHits);
          if( getInt > 0 ) { m_primary_bhitSplit->FillHists( trk,weight); }
        }
        // is b-layer hit merged? Need to go to cluster and count!!
        // already done above for all hits on track!
	if (trk->isAvailable<std::vector<int>>("pixClustPriPartContent"))
	  if( (trk->auxdata< std::vector<int> >( "pixClusPriPartContent" )).at(0) > 1 ) {
          m_primary_bhitMerge->FillHists( trk,weight);
        }
      } // splitonhits
    } else if(m_doDups) {
      m_primary_dup->FillHists( trk,weight);
    } 
  } else if( m_doSecondary && xAOD::TrackHelper::isSecondary( trk ) ){
    if( !m_secondary->hasBeenUsed( trk ) ){
      m_secondary->FillHists( trk,weight);
    } else if(m_doDups) {
      m_secondary_dup->FillHists( trk,weight);
    }
  }

} // FillHists

