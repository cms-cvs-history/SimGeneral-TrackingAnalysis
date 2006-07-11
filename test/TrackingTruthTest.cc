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
  event.getByType(TruthTrackContainer );
  event.getByType(TruthVertexContainer);

  const TrackingParticleCollection *tPC = TruthTrackContainer.product();
  const TrackingVertexCollection   *tVC = TruthVertexContainer.product();
  cout << "Found " << tPC->size() << " tracks and " << tVC->size() << " vertices."<<endl;
  
  for (TrackingVertexCollection::const_iterator v = tVC -> begin(); v != tVC ->
      end(); ++v) {
    cout << " Vertex Position " << v-> position() << endl;
    for (TrackingParticleContainer::iterator t =  v -> tracks_begin(); 
                                             t != v -> tracks_end(); ++t) {
      cout << "  Track PDG ID " << (*t)->pdgId() << endl;
      // Get info from SimTrack
    }  
  }  
  
  
//  edm::Handle<edm::EmbdSimTrackContainer> G4TrkContainer;
//  e.getByType(G4TrkContainer);
  
//  delete tPC;
  
}
