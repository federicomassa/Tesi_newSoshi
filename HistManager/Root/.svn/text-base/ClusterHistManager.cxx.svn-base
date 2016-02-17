/******************************************
 * Plot of tracks and tracking quantities
 *
 * G. Facini
 * Sep 10 10:03:25 CEST 2014
 *
 * Migrate into rel20. -- Mar  2 11:44:17 CET 2015 S. Tsuno
 *
 ******************************************/

#include <HistManager/ClusterHistManager.h>

ClusterHistManager::ClusterHistManager(TString name, int layer, bool more) {
   m_name = name;
   m_layer = layer;
   m_more = more;
   m_all          = new ClusterHists(name);
   m_single       = new ClusterHists(name+"_single");
   m_single_split = new ClusterHists(name+"_single_split");
   m_merged       = new ClusterHists(name+"_merged");
   m_merged_split = new ClusterHists(name+"_merged_split");
   if ( m_more ) {
     m_merged2       = new ClusterHists(name+"_merged2");
     m_merged2_split = new ClusterHists(name+"_merged2_split");
     m_merged3       = new ClusterHists(name+"_merged3");
     m_merged3_split = new ClusterHists(name+"_merged3_split");
     m_merged4       = new ClusterHists(name+"_merged4");
     m_merged4_split = new ClusterHists(name+"_merged4_split");
   }
}

ClusterHistManager::~ClusterHistManager() {}

void ClusterHistManager::Init(EL::Worker* wk) {

  std::cout << "Init ClusterHistManager " << m_name << " for layer " << m_layer << std::endl;
  m_all          -> Init( wk );
  m_single       -> Init( wk );
  m_single_split -> Init( wk );
  m_merged       -> Init( wk );
  m_merged_split -> Init( wk );
  if( m_more ) {
    m_merged2       -> Init( wk );
    m_merged2_split -> Init( wk );
    m_merged3       -> Init( wk );
    m_merged3_split -> Init( wk );
    m_merged4       -> Init( wk );
    m_merged4_split -> Init( wk );
  }

} // BookHists

void ClusterHistManager::FillHists(const xAOD::TrackMeasurementValidation* prd, float weight) const {

  // check if have a layer cut and if do apply the requirement
  if( m_layer >= 0 ) {
    if( prd->auxdata< int >( "layer" ) != m_layer ) { return; }
  }

  // use nParticles here or nPrimaryParticles ??
  
  m_all -> FillHists( prd,weight);
  if( prd->auxdata< int >( "nParticles" ) == 1 ) { 
    m_single->FillHists( prd,weight);
    if( prd->auxdata< int >( "isSplit" ) ) { 
      m_single_split->FillHists( prd,weight);
    }
  } else {
    m_merged->FillHists( prd,weight);
    if( prd->auxdata< int >( "isSplit" ) ) { 
      m_merged_split->FillHists( prd,weight);
    }
    if( m_more ) {
      switch (prd->auxdata< int >( "nParticles" ))
      {
        case 2 :
          m_merged2->FillHists( prd,weight);
          if( prd->auxdata< int >( "isSplit" ) ) { 
            m_merged2_split->FillHists( prd,weight);
          }
          break;
          
        case 3 :
          m_merged3->FillHists( prd,weight);
          if( prd->auxdata< int >( "isSplit" ) ) { 
            m_merged3_split->FillHists( prd,weight);
          }
          break;

        default : // 4+
          m_merged4->FillHists( prd,weight);
          if( prd->auxdata< int >( "isSplit" ) ) { 
            m_merged4_split->FillHists( prd,weight);
          }
          break;
      } // switch
    }
  }

} // FillHists

