#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimGeneral/TrackingAnalysis/interface/EncodedTruthId.h"
#include "SimGeneral/TrackingAnalysis/interface/TrackingTruthProducer.h"

#include <map>

using namespace edm;
using namespace std;
using CLHEP::HepLorentzVector;

typedef edm::Ref<edm::HepMCProduct, HepMC::GenParticle > GenParticleRef;
typedef edm::Ref<edm::HepMCProduct, HepMC::GenVertex >   GenVertexRef;

string MessageCategory = "TrackingTruthProducer";

TrackingTruthProducer::TrackingTruthProducer(const edm::ParameterSet &conf) {
  produces<TrackingVertexCollection>();
  produces<TrackingParticleCollection>();

  conf_ = conf;
  distanceCut_           = conf_.getParameter<double>("vertexDistanceCut");
  dataLabels_            = conf_.getParameter<vector<string> >("HepMCDataLabels");
  simHitLabel_           = conf_.getParameter<string>("simHitLabel");
  hitLabelsVector_       = conf_.getParameter<vector<string> >("TrackerHitLabels");
  volumeRadius_          = conf_.getParameter<double>("volumeRadius");
  volumeZ_               = conf_.getParameter<double>("volumeZ");
  discardOutVolume_      = conf_.getParameter<bool>("discardOutVolume");
  discardHitsFromDeltas_ = conf_.getParameter<bool>("DiscardHitsFromDeltas");

  edm::LogInfo (MessageCategory) << "Setting up TrackingTruthProducer";
  edm::LogInfo (MessageCategory) << "Vertex distance cut set to " << distanceCut_  << " mm";
  edm::LogInfo (MessageCategory) << "Volume radius set to "       << volumeRadius_ << " mm";
  edm::LogInfo (MessageCategory) << "Volume Z      set to "       << volumeZ_      << " mm";
  edm::LogInfo (MessageCategory) << "Discard out of volume? "     << discardOutVolume_;
  edm::LogInfo (MessageCategory) << "Discard Hits from Deltas? "  << discardHitsFromDeltas_;

}

void TrackingTruthProducer::produce(Event &event, const EventSetup &) {
//  TimerStack timers;  // Don't need the timers now, left for example
//  timers.push("TrackingTruth:Producer");
//  timers.push("TrackingTruth:Setup");
  // Get information out of event record
  edm::Handle<edm::HepMCProduct> hepMC;
  for (vector<string>::const_iterator source = dataLabels_.begin();
       source != dataLabels_.end(); ++source) {
    try {
      event.getByLabel(*source,hepMC);
      edm::LogInfo (MessageCategory) << "Using HepMC source " << *source;
      break;
    } catch (std::exception &e) {
      // No error since we have other sources to try
    }
  }

  const edm::HepMCProduct *mcp = hepMC.product();
  if (mcp == 0) {
    edm::LogWarning (MessageCategory) << "No HepMC source found";
    return;
  }
  const HepMC::GenEvent *genEvent = mcp -> GetEvent();

  edm::Handle<CrossingFrame> cf;
  try {
    event.getByType(cf);
  } catch (std::exception &e) {
    edm::LogWarning (MessageCategory) << "Crossing frame not found.";
    return;
  }
  std::auto_ptr<MixCollection<SimTrack> >   trackCollection(new MixCollection<SimTrack>(cf.product()));
  std::auto_ptr<MixCollection<SimVertex> > vertexCollection(new MixCollection<SimVertex>(cf.product()));
  std::auto_ptr<MixCollection<PSimHit> >      hitCollection(new MixCollection<PSimHit>(cf.product(),hitLabelsVector_));

// Create collections of things we will put in event,
  auto_ptr<TrackingParticleCollection> tPC(new TrackingParticleCollection);
  auto_ptr<TrackingVertexCollection>   tVC(new TrackingVertexCollection  );

// Get references before put so we can cross reference
  TrackingParticleRefProd refTPC = event.getRefBeforePut<TrackingParticleCollection>();
  TrackingVertexRefProd   refTVC = event.getRefBeforePut<TrackingVertexCollection>();

  map<EncodedTruthId,EncodedTruthId> simTrack_sourceV; // Encoded SimTrack to encoded source vertex
  map<EncodedTruthId,int>            simTrack_tP;      // Encoded SimTrack to TrackingParticle index
  multimap<EncodedTruthId,PSimHit>   simTrack_hit;
  typedef multimap<EncodedTruthId,PSimHit>::const_iterator hitItr;
//  timers.pop();
  for (MixCollection<PSimHit>::MixItr hit = hitCollection -> begin();
       hit != hitCollection -> end(); ++hit) {
    EncodedTruthId simTrackId = EncodedTruthId(hit->eventId(),hit->trackId());
    simTrack_hit.insert(make_pair(simTrackId,*hit));
  }

  for (MixCollection<SimTrack>::MixItr itP = trackCollection->begin();
       itP !=  trackCollection->end(); ++itP){
    int                       q = (int)(itP -> charge()); // Check this
    CLHEP::HepLorentzVector   p = itP -> momentum();
    unsigned int     simtrackId = itP -> trackId();
    int                 genPart = itP -> genpartIndex(); // The HepMC particle number
    int                 genVert = itP -> vertIndex();    // The SimVertex #
    int                   pdgId = itP -> type();
    EncodedEventId trackEventId = itP -> eventId();
    EncodedTruthId      trackId = EncodedTruthId(trackEventId,simtrackId);

    bool signalEvent = (trackEventId.event() == 0 && trackEventId.bunchCrossing() == 0);
    const TrackingParticle::LorentzVector theMomentum(p.x(), p.y(), p.z(), p.t());
    double  time = 0;

    const HepMC::GenParticle *gp = 0;

    if (genPart >= 0 && signalEvent) {
      gp = genEvent -> barcode_to_particle(genPart);  // Pointer to the generating particle.
      pdgId = gp -> pdg_id();
    }

    math::XYZPoint theVertex;

    if (genVert >= 0){ // Add to useful maps
      EncodedTruthId vertexId = EncodedTruthId(trackEventId,genVert);
      simTrack_sourceV.insert(make_pair(trackId,vertexId));
    }

    TrackingParticle tp(q, theMomentum, theVertex, time, pdgId, trackEventId);

// Counting the TP hits using the layers (as in ORCA).
// Does seem to find less hits. maybe b/c layer is a number now, not a pointer
    int totsimhit = 0;
    int oldlay = 0;
    int newlay = 0;
    int olddet = 0;
    int newdet = 0;

// Using simTrack_hit map makes this very fast
    for (hitItr iHit  = simTrack_hit.lower_bound(trackId);
                iHit != simTrack_hit.upper_bound(trackId); ++iHit) {
      PSimHit hit = iHit->second;
      float pratio = hit.pabs()/(itP->momentum().v().mag());

// Discard hits from delta rays if requested

      if (!discardHitsFromDeltas_ || ( discardHitsFromDeltas_ &&  0.5 < pratio && pratio < 2) ) {
        tp.addPSimHit(hit);
        unsigned int detid = hit.detUnitId();
        DetId detId = DetId(detid);
        oldlay = newlay;
        olddet = newdet;
        newlay = LayerFromDetid(detid);
        newdet = detId.subdetId();

// Count hits using layers for glued detectors

        if (oldlay != newlay || (oldlay==newlay && olddet!=newdet) ) {
          totsimhit++;
        }
      }
    }
    tp.setMatchedHit(totsimhit);

    tp.addG4Track(*itP);
    if (genPart >= 0 && signalEvent) {
      tp.addGenParticle(GenParticleRef(hepMC,genPart));
    }

// Add indices to map and add to collection
    simTrack_tP.insert(make_pair(trackId,tPC->size()));
    tPC -> push_back(tp);
  } // Loop on MixCollection<SimTrack>

// Find and loop over EmbdSimVertex vertices

  int vertexIndex = 0;        // Needed for
  int oldTrigger = -1;        // renumbering
  int oldBX      = -999999;   // of vertices
  for (MixCollection<SimVertex>::MixItr itV = vertexCollection->begin();
       itV != vertexCollection->end(); ++itV) {

    CLHEP::HepLorentzVector position = itV -> position();  // Get position of ESV
    bool inVolume = (position.perp() < volumeRadius_ && abs(position.z()) < volumeZ_); // In or out of Tracker
    if (!inVolume && discardOutVolume_) { continue; }        // Skip if desired

    EncodedEventId vertEvtId = itV -> eventId();

// Begin renumbering vertices if we move from signal to pileup or change bunch crossings
    if (oldTrigger !=  itV.getTrigger() || oldBX !=  vertEvtId.bunchCrossing()) {
      vertexIndex = 0;
      oldTrigger =  itV.getTrigger();
      oldBX =  vertEvtId.bunchCrossing();
    }
    EncodedTruthId vertexId  = EncodedTruthId(vertEvtId,vertexIndex);

// Figure out the barcode of the HepMC Vertex if there is one by
// getting incoming SimTtrack (if any), finding corresponding HepMC track and
// then decay (HepMC) vertex of that track.  HepMC data only exists for signal sub-event
    int vertexBarcode = 0;
    unsigned int vtxParent = itV -> parentIndex();
    if (vtxParent >= 0 && itV.getTrigger() ) {
      for (MixCollection<SimTrack>::MixItr itP = trackCollection->begin(); itP != trackCollection->end(); ++itP){
        if (vtxParent==itP->trackId() && itP->eventId() == vertEvtId){
          int parentBC = itP->genpartIndex();
          HepMC::GenParticle *parentParticle = genEvent -> barcode_to_particle(parentBC);
          if (parentParticle != 0) {
            HepMC::GenVertex *hmpv = parentParticle -> end_vertex();
            if (hmpv != 0) {
              vertexBarcode = hmpv  -> barcode();
            }
          }
          break;
        }
      }
    }

// Find closest vertex to this one in same sub-event, save in nearestVertex
    int indexTV = 0;
    double closest = 9e99;
    TrackingVertexCollection::iterator nearestVertex;

    int tmpTV = 0;
    for (TrackingVertexCollection::iterator iTrkVtx = tVC -> begin(); iTrkVtx != tVC ->end(); ++iTrkVtx, ++tmpTV) {
      double distance = HepLorentzVector(iTrkVtx -> position() - position).v().mag();
      if (distance <= closest && vertEvtId == iTrkVtx -> eventId()) { // flag which one so we can associate them
        closest = distance;
        nearestVertex = iTrkVtx;
        indexTV = tmpTV;
      }
    }

// If outside cutoff, create another TrackingVertex, set nearestVertex to it

    if (closest > distanceCut_) {
      indexTV = tVC -> size();
      tVC -> push_back(TrackingVertex(position,inVolume,vertEvtId));
      nearestVertex = --(tVC -> end());  // Last entry of vector
    }

// Add data to closest vertex

    (*nearestVertex).addG4Vertex(*itV); // Add G4 vertex
    if (vertexBarcode != 0) {
      (*nearestVertex).addGenVertex(GenVertexRef(hepMC,vertexBarcode)); // Add HepMC vertex
    }

// Identify and add child tracks
    for (std::map<EncodedTruthId,EncodedTruthId>::iterator mapIndex = simTrack_sourceV.begin();
         mapIndex != simTrack_sourceV.end(); ++mapIndex) {
      EncodedTruthId mapTrackId  = mapIndex -> first;
      EncodedTruthId mapVertexId = mapIndex -> second;
      if (mapVertexId == vertexId) {
        if (simTrack_tP.count(mapTrackId)) {
          int indexTP = simTrack_tP[mapTrackId];
          (*nearestVertex).addDaughterTrack(TrackingParticleRef(refTPC,indexTP));
          (tPC->at(indexTP)).setParentVertex(TrackingVertexRef(refTVC,indexTV));
          const CLHEP::HepLorentzVector &v = (*nearestVertex).position();
          math::XYZPoint xyz = math::XYZPoint(v.x(), v.y(), v.z());
          double t = v.t();
          (tPC->at(indexTP)).setVertex(xyz,t);
        }
      }
    }

// Identify and add parent tracks
    if (vtxParent > 0) {
      EncodedTruthId trackId =  EncodedTruthId(vertEvtId,vtxParent);
      if (simTrack_tP.count(trackId) > 0) {
        int indexTP = simTrack_tP[trackId];
        (tPC->at(indexTP)).addDecayVertex(TrackingVertexRef(refTVC,indexTV));
        (*nearestVertex).addParentTrack(TrackingParticleRef(refTPC,indexTP));
      }
    }
    ++vertexIndex;
  } // Loop on MixCollection<SimVertex>

  edm::LogInfo(MessageCategory) << "TrackingTruth found "  << tVC -> size()
                                << " unique vertices and " << tPC -> size() << " tracks.";

// Put TrackingParticles and TrackingVertices in event
  event.put(tPC);
  event.put(tVC);
//  timers.pop();
//  timers.pop();
}

int TrackingTruthProducer::LayerFromDetid(const unsigned int& detid ) {
  DetId detId = DetId(detid);
  int layerNumber=0;
  unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());
  if ( subdetId == StripSubdetector::TIB)
    {
      TIBDetId tibid(detId.rawId());
      layerNumber = tibid.layer();
    }
  else if ( subdetId ==  StripSubdetector::TOB )
    {
      TOBDetId tobid(detId.rawId());
      layerNumber = tobid.layer();
    }
  else if ( subdetId ==  StripSubdetector::TID)
    {
      TIDDetId tidid(detId.rawId());
      layerNumber = tidid.wheel();
    }
  else if ( subdetId ==  StripSubdetector::TEC )
    {
      TECDetId tecid(detId.rawId());
      layerNumber = tecid.wheel();
    }
  else if ( subdetId ==  PixelSubdetector::PixelBarrel )
    {
      PXBDetId pxbid(detId.rawId());
      layerNumber = pxbid.layer();
    }
  else if ( subdetId ==  PixelSubdetector::PixelEndcap )
    {
      PXFDetId pxfid(detId.rawId());
      layerNumber = pxfid.disk();
    }
  else
    edm::LogVerbatim("TrackingTruthProducer") << "Unknown subdetid: " <<  subdetId;

  return layerNumber;
}

DEFINE_FWK_MODULE(TrackingTruthProducer);
