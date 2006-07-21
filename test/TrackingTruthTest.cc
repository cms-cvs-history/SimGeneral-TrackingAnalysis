#include "SimGeneral/TrackingAnalysis/test/TrackingTruthTest.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackVertexMap.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

typedef edm::RefVector< std::vector<TrackingParticle> > TrackingParticleContainer;
typedef std::vector<TrackingParticle> TrackingParticleCollection;

TrackingTruthTest::TrackingTruthTest(const edm::ParameterSet& conf){
  conf_ = conf;
}

void TrackingTruthTest::analyze(const edm::Event& event, const edm::EventSetup& c){
  using namespace std;

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  edm::Handle<TrackingVertexCollection>    TruthVertexContainer;
  edm::Handle<VertexTrackAssociationCollection>    vertexToTrackHandle;
  event.getByType(TruthTrackContainer );
  event.getByType(TruthVertexContainer);
  event.getByType(vertexToTrackHandle);

  const TrackingParticleCollection *tPC   = TruthTrackContainer.product();
  const TrackingVertexCollection   *tVC   = TruthVertexContainer.product();
  const VertexTrackAssociationCollection   *vTMap = vertexToTrackHandle.product();

  cout << "Found " << tPC->size() << " tracks and " << tVC->size() << " vertices."<<endl;
  cout << "Found " << vTMap -> size() << " vertex to track associations " <<
      endl;  
// Loop over TrackingParticle's  
  int count = 0; 
  
  cout << "Dumping sample track info" << endl;
  for (TrackingParticleCollection::const_iterator t = tPC -> begin(); 
      t != tPC -> end(); ++t, ++count) {
    for (TrackingParticle::genp_iterator hepT = t -> genParticle_begin();
         hepT !=  t -> genParticle_end(); ++hepT) {
      cout << "  Gen Track PDG ID   " <<  (*hepT)->pdg_id() << endl;    
      cout << "  Gen Track Momentum " << (*hepT)->momentum() << endl;    
    }
    for (TrackingParticle::g4t_iterator g4T = t -> g4Track_begin();
         g4T !=  t -> g4Track_end(); ++g4T) {
      cout << "  G4  Track Momentum " << (*g4T)->momentum() << endl;    
    }
  }  
 
// Loop over TrackingVertex's  
  
  cout << "Dumping some vertex info" << endl;
  for (TrackingVertexCollection::const_iterator v = tVC -> begin(); v != tVC ->
      end(); ++v) {
    cout << " Vertex Position " << v-> position() << endl;
    for (TrackingParticleContainer::iterator t =  v -> tracks_begin(); 
                                             t != v -> tracks_end(); ++t) {
      cout << "  Track PDG ID " << (*t)->pdgId() << endl;
      // Get info from SimTrack
      for (TrackingParticle::g4t_iterator g4T = (*t) -> g4Track_begin();
           g4T != (*t) -> g4Track_end(); ++g4T) {
        cout << "   G4 Track Momentum " << (*g4T)->momentum() << endl;    
      }
    }  
  }  
  
// Loop over vertex to track map
  cout << "Getting info from map" << endl;
  for (VertexTrackAssociationCollection::const_iterator v2t = vTMap->begin(); v2t
      != vTMap->end(); ++v2t) {
    TrackingVertexRef         vertex = v2t -> key;
    
    TrackingParticleRefVector tracks = v2t -> val;
    cout << " Vertex Position " << vertex->position() << endl;
    for (TrackingParticleContainer::const_iterator track = tracks.begin();
         track != tracks.end(); ++track) {
      cout << "  Track PDG ID " << (*track)->pdgId() << endl;
      for (TrackingParticle::g4t_iterator g4T = (*track) -> g4Track_begin();
           g4T != (*track) -> g4Track_end(); ++g4T) {
        cout << "   G4 Track Momentum " << (*g4T)->momentum() << endl;    
      }
    }     
  }  
      
  cout << "Done with event" << endl;
  
}
