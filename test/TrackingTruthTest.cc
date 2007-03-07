#include "SimGeneral/TrackingAnalysis/test/TrackingTruthTest.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

typedef edm::RefVector< std::vector<TrackingParticle> > TrackingParticleContainer;
typedef std::vector<TrackingParticle>                   TrackingParticleCollection;

typedef TrackingParticleRefVector::iterator               tp_iterator;
typedef TrackingParticle::g4t_iterator                   g4t_iterator;
typedef TrackingParticle::genp_iterator                 genp_iterator;
typedef TrackingVertex::genv_iterator                   genv_iterator;
typedef TrackingVertex::g4v_iterator                     g4v_iterator;

TrackingTruthTest::TrackingTruthTest(const edm::ParameterSet& conf){
  conf_ = conf;
}

void TrackingTruthTest::analyze(const edm::Event& event, const edm::EventSetup& c){
  using namespace std;

  edm::Handle<TrackingParticleCollection>  TruthTrackContainer ;
  edm::Handle<TrackingVertexCollection>    TruthVertexContainer;
  event.getByType(TruthTrackContainer );
  event.getByType(TruthVertexContainer);

  const TrackingParticleCollection *tPC   = TruthTrackContainer.product();
  const TrackingVertexCollection   *tVC   = TruthVertexContainer.product();

// Get and print HepMC event for comparison 
  edm::Handle<edm::HepMCProduct> hepMC;
  event.getByLabel("VtxSmeared",hepMC);
  const edm::HepMCProduct *mcp = hepMC.product();
  const HepMC::GenEvent *genEvent = mcp -> GetEvent();
  genEvent -> print();

  edm::Handle<CrossingFrame> cf;
  event.getByType(cf);      
  std::auto_ptr<MixCollection<SimTrack> >   trackCollection (new MixCollection<SimTrack>(cf.product()));
  std::auto_ptr<MixCollection<SimVertex> > vertexCollection (new MixCollection<SimVertex>(cf.product()));

// Dump GEANT tracks and vertices  
  
/*  cout << "Dumping GEANT tracks" << endl;
  cout << "PDG Momentum          Vtx genPart" << endl;
  for (MixCollection<SimTrack>::MixItr itP = trackCollection->begin(); itP !=  trackCollection->end(); ++itP){
    cout  << itP -> eventId().bunchCrossing() << " " << itP ->  eventId().event() << " " << (*itP) << endl;
  }  
  

  cout << "Dumping GEANT vertices" << endl;
  cout << "Position         track #" << endl;
  for (MixCollection<SimVertex>::MixItr itV = vertexCollection->begin(); itV != vertexCollection->end(); ++itV) {
    cout << itV -> eventId().bunchCrossing() << " " << itV ->  eventId().event() << " " << (*itV) << endl;
  }  
*/
             
  cout << "Found " << tPC -> size() << " tracks and " << tVC -> size() << " vertices." <<endl;
// Loop over TrackingParticle's

  cout << "Dumping sample track info" << endl;
  for (TrackingParticleCollection::const_iterator t = tPC -> begin(); t != tPC -> end(); ++t) {
    
    // Compare momenta from sources
    cout << "T.P.   Track Momentum, q , ID, & Event # " 
          << t -> p4()    << " " << t -> charge() << " "
          << t -> pdgId() << " " 
          << t -> eventId().bunchCrossing() << "." << t -> eventId().event() << endl;
    cout << " Hits for this track: " << t -> trackPSimHit().size() << endl;

    for (TrackingParticle::genp_iterator hepT = t -> genParticle_begin();
         hepT !=  t -> genParticle_end(); ++hepT) {
      cout << " HepMC Track Momentum " << (*hepT)->momentum() << endl;    
    }
    for (TrackingParticle::g4t_iterator g4T = t -> g4Track_begin();
         g4T !=  t -> g4Track_end(); ++g4T) {
      cout << " Geant Track Momentum  " << g4T->momentum() << endl;   
      cout << " Geant Track ID & type " << g4T->trackId() << " " << g4T->type() << endl;   
      if (g4T->type() !=  t -> pdgId()) {
        cout << " Mismatch b/t TrackingParticle and Geant types" << endl;
      }
    }

    // Compare starting and ending points
    TrackingVertexRef parentV = t -> parentVertex();
//    TrackingVertexRef decayV  = t -> decayVertex();
    
    cout << " Track start position " << t -> vertex() << endl;
    if (parentV.isNull()) {
      cout << " No parent vertex" << endl;
    } else {  
      cout << " Parent  vtx position " << parentV -> position() << endl;
    }  
// Reinstate in 1.4.0
//   for (TrackingParticle::tv_iterator decayV = t -> decayVertices_begin();
//        decayV !=  t -> decayVertices_end(); ++decayV) {
//     cout << " Decay   vtx position " << (*decayV)  -> position() << endl;
//   }
  }  // End loop over TrackingParticle
 
// Loop over TrackingVertex's
  
  cout << "Dumping sample vertex info" << endl;
  for (TrackingVertexCollection::const_iterator v = tVC -> begin(); v != tVC -> end(); ++v) {
    cout << " Vertex Position & Event #" << v -> position() << " " << v -> eventId().bunchCrossing() << "." << v -> eventId().event() << endl; 
    cout << "  Associated with " << v -> daughterTracks().size() << " tracks" << endl;
    // Get Geant and HepMC positions
    for (genv_iterator genV = v -> genVertices_begin(); genV != v -> genVertices_end(); ++genV) {
      cout << "  HepMC vertex position " << (*(*genV)).position() << endl; 
    }  
    for (g4v_iterator g4V = v -> g4Vertices_begin(); g4V != v -> g4Vertices_end(); ++g4V) {
      cout << "  Geant vertex position " << (*g4V).position() << endl; 
      // Probably empty all the time, currently
    }  
    
    // Loop over daughter track(s)
    for (tp_iterator iTP = v -> daughterTracks_begin(); iTP != v -> daughterTracks_end(); ++iTP) {
      cout << "  Daughter starts:      " << (*(*iTP)).vertex();
      for (g4t_iterator g4T  = (*(*iTP)).g4Track_begin(); g4T != (*(*iTP)).g4Track_end(); ++g4T) {
        cout << " p " << g4T->momentum();    
      }
      cout << endl;
    }   
    
    // Loop over source track(s) (can be multiple since vertices are collapsed)
    for (tp_iterator iTP = v -> sourceTracks_begin(); iTP != v -> sourceTracks_end(); ++iTP) {
      cout << "  Source   starts: " << (*(*iTP)).vertex();
      for (g4t_iterator g4T  = (*iTP)->g4Track_begin(); g4T != (*iTP)->g4Track_end(); ++g4T) {
        cout << ", p " <<  g4T ->momentum();    
      }
      cout << endl;
    }   
  }  // End loop over TrackingVertex
  
}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TrackingTruthTest);
