#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixModFactory.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"

#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "SimGeneral/TrackingAnalysis/interface/TrackingTruthProducer.h"

typedef edm::Ref<edm::HepMCProduct, HepMC::GenParticle > GenParticleRef;
typedef edm::Ref<edm::HepMCProduct, HepMC::GenVertex >   GenVertexRef;
typedef math::XYZTLorentzVectorD    LorentzVector;
typedef math::XYZPoint Vector;

TrackingTruthProducer::TrackingTruthProducer(edm::ParameterSet const& config, edm::EDProducer& mixMod) :
    // Initialize global parameters
    minBunch_(config.getParameter<int>("minBunch")),
    maxBunch_(config.getParameter<int>("maxBunch")),
    dataLabels_(config.getParameter<std::vector<std::string> >("HepMCDataLabels")),
    useMultipleHepMCLabels_(config.getParameter<bool>("useMultipleHepMCLabels")),
    distanceCut_(config.getParameter<double>("vertexDistanceCut")),
    volumeRadius_(config.getParameter<double>("volumeRadius")),
    volumeZ_(config.getParameter<double>("volumeZ")),
    mergedBremsstrahlung_(config.getParameter<bool>("mergedBremsstrahlung")),
    removeDeadModules_(config.getParameter<bool>("removeDeadModules")),
    simHitLabel_(config.getParameter<std::string>("simHitLabel")),
    MessageCategory_("TrackingTruthProducer"),
    pSimHitSelector_(config),
    // Initialize selection for building TrackingParticles
    selectorFlag_(config.exists("select")),
    selector_(selectorFlag_ ?
        TrackingParticleSelector(
		config.getParameterSet("select").getParameter<double>("ptMinTP"),
		config.getParameterSet("select").getParameter<double>("minRapidityTP"),
		config.getParameterSet("select").getParameter<double>("maxRapidityTP"),
		config.getParameterSet("select").getParameter<double>("tipTP"),
		config.getParameterSet("select").getParameter<double>("lipTP"),
		config.getParameterSet("select").getParameter<int>("minHitTP"),
		config.getParameterSet("select").getParameter<bool>("signalOnlyTP"),
		config.getParameterSet("select").getParameter<bool>("chargedOnlyTP"),
		config.getParameterSet("select").getParameter<bool>("stableOnlyTP"),
		config.getParameterSet("select").getParameter<std::vector<int> >("pdgIdTP")

       ) : TrackingParticleSelector()) {

    edm::LogInfo (MessageCategory_) << "Setting up TrackingTruthProducer";
    edm::LogInfo (MessageCategory_) << "Vertex distance cut set to " << distanceCut_  << " mm";
    edm::LogInfo (MessageCategory_) << "Volume radius set to "       << volumeRadius_ << " mm";
    edm::LogInfo (MessageCategory_) << "Volume Z      set to "       << volumeZ_      << " mm";

    if(useMultipleHepMCLabels_) edm::LogInfo (MessageCategory_) << "Collecting generator information from pileup.";
    if(mergedBremsstrahlung_) edm::LogInfo (MessageCategory_) << "Merging electrom bremsstralung";
    if(removeDeadModules_) edm::LogInfo (MessageCategory_) << "Removing psimhit from dead modules";
 

    for(int i = minBunch_; i <= maxBunch_; ++i) {
      std::ostringstream ostr;
      if(i != 0) ostr << i;
      std::string const& str = ostr.str();

      mixMod.produces<TrackingVertexCollection>(str);
      mixMod.produces<TrackingParticleCollection>(str);

      if(mergedBremsstrahlung_) {
        mixMod.produces<TrackingVertexCollection>("MergedTrackTruth" + str);
        mixMod.produces<TrackingParticleCollection>("MergedTrackTruth" + str);
      }
    }
}


void TrackingTruthProducer::initializeEvent(edm::Event const& event, edm::EventSetup const&) {

    // Clean the list of hepmc products
    hepMCProducts_.clear();

    // Collect all the HepMCProducts
    edm::Handle<edm::HepMCProduct> hepMCHandle;

    for(std::vector<std::string>::const_iterator source = dataLabels_.begin(); source != dataLabels_.end(); ++source) {
        if(event.getByLabel(*source, hepMCHandle)) {
            hepMCProducts_.push_back(hepMCHandle);
            edm::LogInfo (MessageCategory_) << "Using HepMC source " << *source;
            if(!useMultipleHepMCLabels_) break;
        }
    }

    if(hepMCProducts_.empty()) {
        edm::LogWarning (MessageCategory_) << "No HepMC source found";
    } else if(hepMCProducts_.size() > 1 || useMultipleHepMCLabels_) {
        edm::LogInfo (MessageCategory_) << "You are using more than one HepMC source.";
        edm::LogInfo (MessageCategory_) << "If the labels are not in the same order as the events in the crossing frame (i.e. signal, pileup(s) ) ";
        edm::LogInfo (MessageCategory_) << "or there are fewer labels than events in the crossing frame";
        edm::LogInfo (MessageCategory_) << MessageCategory_ << " may try to access data in the wrong HepMCProduct and crash.";
    }
}

void TrackingTruthProducer::initializeBunchCrossing(edm::Event const& event, edm::EventSetup const& setup, int bunchCrossing) {
  if(bunchCrossing < minBunch_ || bunchCrossing > maxBunch_) {
    return;
  }
   // If this is bunch crossing zero, process the signal event.
  if(bunchCrossing == 0) {
    accumulateEvent(event, setup);
  }
}

void TrackingTruthProducer::finalizeBunchCrossing(edm::Event& event, edm::EventSetup const&, int bunchCrossing) {
  if(bunchCrossing < minBunch_ || bunchCrossing > maxBunch_) {
    return;
  }
  std::ostringstream ostr;
  if(bunchCrossing != 0) ostr << bunchCrossing;
  std::string const& str = ostr.str();
  std::string const mstr = "MergedTrackTruth" + str;

  // Create collections of things we will put in event,
  std::auto_ptr<TrackingParticleCollection> trackingParticles(new TrackingParticleCollection);
  std::auto_ptr<TrackingVertexCollection> trackingVertexes(new TrackingVertexCollection);

  EncodedTruthIdToIndexes trackIdToHits;
  EncodedTruthIdToIndex   trackIdToIndex;
  EncodedTruthIdToIndex   vertexIdToIndex;

  // Create a multimap between trackId and hit indices
  associator(pSimHits_, trackIdToHits);

  // Create a map between trackId and track index
  associator(simTracks_, trackIdToIndex);

  // Create a map between vertexId and vertex index
  associator(simVertexes_, vertexIdToIndex);

  createTrackingTruth(event, trackIdToHits, trackIdToIndex, vertexIdToIndex, *trackingParticles, *trackingVertexes);

  if(mergedBremsstrahlung_) {
    // Create collections of things we will put in event,
    std::auto_ptr<TrackingParticleCollection> mergedTrackingParticles(new TrackingParticleCollection);
    std::auto_ptr<TrackingVertexCollection> mergedTrackingVertexes(new TrackingVertexCollection);

    // Merged brem electrons
    mergeBremsstrahlung(event, *trackingParticles, *trackingVertexes, *mergedTrackingParticles, *mergedTrackingVertexes);

    // Put merged TrackingParticles and TrackingVertices in event
    event.put(mergedTrackingParticles, mstr);
    event.put(mergedTrackingVertexes, mstr);
  }

  // Put TrackingParticles and TrackingVertices in event
  event.put(trackingParticles, str);
  event.put(trackingVertexes, str);

  // clear and shrink the vectors
  simTracks_ = std::move(std::vector<SimTrack>());
  simVertexes_ = std::move(std::vector<SimVertex>());
  pSimHits_ = std::move(std::vector<PSimHit>());
}

void TrackingTruthProducer::accumulate(edm::Event const& event, edm::EventSetup const&) {
   // We do nothing here, because the signal event is processed when the zero bunch crossing pileups are processed.
}

void TrackingTruthProducer::accumulate(PileUpEventPrincipal const& event, edm::EventSetup const& setup) {
  int bunchCrossing = event.bunchCrossing();
  if(bunchCrossing < minBunch_ || bunchCrossing > maxBunch_) {
    return;
  }
  edm::Event const ev(const_cast<edm::EventPrincipal&>(event.principal()), edm::ModuleDescription());
  accumulateEvent(ev, setup);
}

void TrackingTruthProducer::accumulateEvent(edm::Event const& event, edm::EventSetup const& setup) {

   edm::InputTag tag(simHitLabel_);

   // Collect all the simTracks
   edm::Handle<std::vector<SimTrack> > hSimTracks;
   event.getByLabel(tag, hSimTracks);
   std::vector<SimTrack> const& simTracks = *hSimTracks.product();
   simTracks_.reserve(simTracks_.capacity() + simTracks.size());
   simTracks_.insert(simTracks_.end(), simTracks.begin(), simTracks.end());

   // Collect all the simVertexes
   edm::Handle<std::vector<SimVertex> > hSimVertexes;
   event.getByLabel(tag, hSimVertexes);
   std::vector<SimVertex> const& simVertexes = *hSimVertexes.product();
   simVertexes_.reserve(simVertexes_.capacity() + simVertexes.size());
   simVertexes_.insert(simVertexes_.end(), simVertexes.begin(), simVertexes.end());

    // Collect all the psimhits
    if(removeDeadModules_) {
        pSimHitSelector_.selectPixel(pSimHits_, event, setup);
        pSimHitSelector_.selectTracker(pSimHits_, event, setup);
        pSimHitSelector_.selectMuon(pSimHits_, event, setup);
    } else {
        pSimHitSelector_.select(pSimHits_, event, setup);
    }
}

void TrackingTruthProducer::finalizeEvent(edm::Event& event, edm::EventSetup const& setup) {
  // Nothing to do, because everything was done at the end of each bunch crossing.
}


void TrackingTruthProducer::associator(
    std::vector<PSimHit> const& pSimHits,
    EncodedTruthIdToIndexes& association)
{
    // Clear the association map
    association.clear();

    // Create a association from simtracks to overall index in the mix collection
    for(std::size_t i = 0; i < pSimHits.size(); ++i) {
        EncodedTruthIdToIndexes::key_type objectId = EncodedTruthIdToIndexes::key_type(pSimHits[i].eventId(), pSimHits[i].trackId());
        association.insert(std::make_pair(objectId, i));
    }
}


void TrackingTruthProducer::associator(
    SimTracks const& simTracks,
    EncodedTruthIdToIndex& association) {

    // Clear the association map
    association.clear();
    // Create a association from simtracks to overall index in the mix collection
    int index = 0;
    for(auto const& simTrack : simTracks) {
      EncodedTruthId objectId = EncodedTruthId(simTrack.eventId(), simTrack.trackId());
      association.insert(std::make_pair(objectId, index));
      ++index;
    }
}


void TrackingTruthProducer::associator(
    SimVertexes const& simVertexes,
    EncodedTruthIdToIndex& association) {

    // Solution to the problem of not having vertexId
    bool useVertexId = true;
    EncodedEventIdToIndex vertexId;
    EncodedEventId oldEventId;
    unsigned int oldVertexId = 0;

    bool first = true;
    // Loop for finding repeated vertexId (vertexId problem hack)
    for(auto const& simVertex : simVertexes) {
        if(first || simVertex.eventId() != oldEventId) {
            oldEventId = simVertex.eventId();
            oldVertexId = simVertex.vertexId();
            first = false;
            continue;
        }

        if(simVertex.vertexId() == oldVertexId) {
            edm::LogWarning(MessageCategory_) << "Multiple vertexId found, no using vertexId.";
            useVertexId = false;
            break;
        }
    }

    int index = 0;

    // Clear the association map
    association.clear();

    // Create a association from simVertexes to overall index
    for(auto const& simVertex : simVertexes) {
      EncodedTruthId objectId;
      if(useVertexId) {
        objectId = EncodedTruthId(simVertex.eventId(), simVertex.vertexId());
      } else {
        objectId = EncodedTruthId(simVertex.eventId(), vertexId[simVertex.eventId()]++);
      }
      association.insert(std::make_pair(objectId, index));
    }
}

void TrackingTruthProducer::mergeBremsstrahlung(edm::Event& event,
                                                TrackingParticleCollection& trackingParticles,
                                                TrackingVertexCollection& trackingVertexes,
                                                TrackingParticleCollection& mergedTrackingParticles,
                                                TrackingVertexCollection& mergedTrackingVertexes) {

    // Get references before put so we can cross reference
    TrackingParticleRefProd refMergedTrackingParticles = event.getRefBeforePut<TrackingParticleCollection>("MergedTrackTruth");
    TrackingVertexRefProd refMergedTrackingVertexes = event.getRefBeforePut<TrackingVertexCollection>("MergedTrackTruth");

    unsigned int index = 0;

    std::set<unsigned int> excludedTV, excludedTP;

    // Merge Bremsstrahlung vertexes
    for(TrackingVertexCollection::iterator iVC = trackingVertexes.begin(); iVC != trackingVertexes.end(); ++iVC, ++index) {
        // Check Bremsstrahlung vertex
        if(isBremsstrahlungVertex(*iVC, trackingParticles)) {
            // Get a pointer to the source track (A Ref<> cannot be use with a product!)
            TrackingParticle* track = &trackingParticles.at(iVC->sourceTracks_begin()->key());
            // Get a Ref<> to the source track
            TrackingParticleRef trackRef = *iVC->sourceTracks_begin();
            // Pointer to electron daughter
            TrackingParticle* daughter = 0;
            // Ref<> to electron daughter
            TrackingParticleRef daughterRef;

            // Select the electron daughter and redirect the photon
            for(TrackingVertex::tp_iterator idaughter = iVC->daughterTracks_begin(); idaughter != iVC->daughterTracks_end(); ++idaughter) {
                TrackingParticle* pointer = &trackingParticles.at(idaughter->key());
                if(std::abs(pointer->pdgId()) == 11) {
                    // Set pointer to the electron daughter
                    daughter = pointer;
                    // Set Ref<> to the electron daughter
                    daughterRef = *idaughter;
                } else if(pointer->pdgId() == 22) {
                    // Delete the photon original parent vertex
                    pointer->clearParentVertex();
                    // Set the new parent vertex to the vertex of the source track
                    pointer->setParentVertex(track->parentVertex());
                    // Get a non-const pointer to the parent vertex
                    TrackingVertex* vertex = &trackingVertexes.at(track->parentVertex().key());
                    // Add the photon to the daughter list of the parent vertex
                    vertex->addDaughterTrack(*idaughter);
                }
            }

            // Add the electron segments from the electron daughter
            // track must not be the same particle as daughter
            // if(track != daughter)
            for(TrackingParticle::g4t_iterator isegment = daughter->g4Track_begin(); isegment != daughter->g4Track_end(); ++isegment)
                track->addG4Track(*isegment);

            // Copy all the simhits to the new track
            for(std::vector<PSimHit>::const_iterator ihit = daughter->pSimHit_begin(); ihit != daughter->pSimHit_end(); ++ihit)
                track->addPSimHit(*ihit);

            // Make a copy of the decay vertexes of the track
            TrackingVertexRefVector decayVertices(track->decayVertices());

            // Clear the decay vertex list
            track->clearDecayVertices();

            // Add the remaining vertexes
            for(TrackingVertexRefVector::const_iterator idecay = decayVertices.begin(); idecay != decayVertices.end(); ++idecay)
                if((*idecay).key() != index) track->addDecayVertex(*idecay);

            // Redirect all the decay source vertexes to those in the electron daughter
            for(TrackingParticle::tv_iterator idecay = daughter->decayVertices_begin(); idecay != daughter->decayVertices_end(); ++idecay) {
                // Add the vertexes to the decay list of the source particles
                track->addDecayVertex(*idecay);
                // Get a reference to decay vertex
                TrackingVertex* vertex = &trackingVertexes.at(idecay->key());
                // Copy all the source tracks from of the decay vertex
                TrackingParticleRefVector sources(vertex->sourceTracks());
                // Clear the source track references
                vertex->clearParentTracks();
                // Add the new source tracks by excluding the one with the segment merged
                for(TrackingVertex::tp_iterator isource = sources.begin(); isource != sources.end(); ++isource)
                    if(*isource != daughterRef)
                        vertex->addParentTrack(*isource);
                // Add the track reference to the list of sources
                vertex->addParentTrack(trackRef);
            }

            // Adding the vertex to the exlusion list
            excludedTV.insert(index);

            // Adding the electron segment tp into the exlusion list
            excludedTP.insert(daughterRef.key());
        }
    }

    edm::LogInfo(MessageCategory_) << "Generating the merged collection." << std::endl;

    // Reserved the same amount of memory for the merged collections
    mergedTrackingParticles.reserve(trackingParticles.size());
    mergedTrackingVertexes.reserve(trackingVertexes.size());

    index = 0;
    std::map<unsigned int, unsigned int> vertexMap;

    // Copy non-excluded vertices discarding parent & child tracks
    for(TrackingVertexCollection::const_iterator iVC = trackingVertexes.begin(); iVC != trackingVertexes.end(); ++iVC, ++index)
    {
        if(excludedTV.find(index) != excludedTV.end()) continue;
        // Save the new location of the non excluded vertexes (take in consideration those were removed)
        vertexMap.insert(std::make_pair(index, mergedTrackingVertexes.size()));
        // Copy those vertexes are not excluded
        TrackingVertex newVertex = (*iVC);
        newVertex.clearDaughterTracks();
        newVertex.clearParentTracks();
        mergedTrackingVertexes.push_back(newVertex);
    }

    index = 0;

    // Copy and cross reference the non-excluded tp to the merged collection
    for(TrackingParticleCollection::const_iterator iTP = trackingParticles.begin(); iTP != trackingParticles.end(); ++iTP, ++index)
    {
        if(excludedTP.find(index) != excludedTP.end()) continue;

        TrackingVertexRef       sourceV = iTP->parentVertex();
        TrackingVertexRefVector decayVs = iTP->decayVertices();
        TrackingParticle newTrack = *iTP;

        newTrack.clearParentVertex();
        newTrack.clearDecayVertices();

        // Set vertex indices for new vertex product and track references in those vertices

        // Index of parent vertex in vertex container
        unsigned int parentIndex = vertexMap[sourceV.key()];
        // Index of this track in track container
        unsigned int tIndex = mergedTrackingParticles.size();

        // Add vertex to track
        newTrack.setParentVertex(TrackingVertexRef(refMergedTrackingVertexes, parentIndex));
        // Add track to vertex
        (mergedTrackingVertexes.at(parentIndex)).addDaughterTrack(TrackingParticleRef(refMergedTrackingParticles, tIndex));

        for(TrackingVertexRefVector::const_iterator iDecayV = decayVs.begin(); iDecayV != decayVs.end(); ++iDecayV)
        {
            // Index of decay vertex in vertex container
            unsigned int daughterIndex = vertexMap[iDecayV->key()];
            // Add vertex to track
            newTrack.addDecayVertex(TrackingVertexRef(refMergedTrackingVertexes, daughterIndex));
            // Add track to vertex
            (mergedTrackingVertexes.at(daughterIndex)).addParentTrack(TrackingParticleRef(refMergedTrackingParticles, tIndex));
        }

        mergedTrackingParticles.push_back(newTrack);
    }
}


bool TrackingTruthProducer::isBremsstrahlungVertex(
    TrackingVertex const& vertex,
    TrackingParticleCollection const& tPC) {

    TrackingParticleRefVector const parents(vertex.sourceTracks());

    // Check for the basic parent conditions
    if(parents.size() != 1)
        return false;

    // Check for the parent particle is a |electron| (electron or positron)
    if(std::abs(tPC.at(parents.begin()->key()).pdgId()) != 11)
        return false;

    unsigned int nElectrons = 0;
    unsigned int nOthers = 0;

    // Loop over the daughter particles and counts the number of |electrons|, others (non photons)
    for(TrackingVertex::tp_iterator it = vertex.daughterTracks_begin(); it != vertex.daughterTracks_end(); ++it) {
        // Stronger rejection for looping particles
        if(parents[0] == *it) {
            return false;
        }

        if(std::abs(tPC.at(it->key()).pdgId()) == 11) {
            ++nElectrons;
        } else if(tPC.at(it->key()).pdgId() != 22) {
            ++nOthers;
        }
    }

    // Condition to be a Bremsstrahlung Vertex
    if(nElectrons == 1 && nOthers == 0)
        return true;

    return false;
}


void TrackingTruthProducer::createTrackingTruth(edm::Event& event,
                                                EncodedTruthIdToIndexes const& trackIdToHits,
                                                EncodedTruthIdToIndex& trackIdToIndex,
                                                EncodedTruthIdToIndex& vertexIdToIndex,
                                                TrackingParticleCollection& trackingParticles,
                                                TrackingVertexCollection& trackingVertexes) {
    // Get references before put so we can cross reference
    TrackingParticleRefProd refTrackingParticles = event.getRefBeforePut<TrackingParticleCollection>();
    TrackingVertexRefProd refTrackingVertexes = event.getRefBeforePut<TrackingVertexCollection>();

    EncodedEventIdToIndex eventIdCounter;

    // Define a container of vetoed traks
    std::map<int,std::size_t> vetoedTracks;

    // Define map between parent simtrack and tv indexes
    std::map<int,std::size_t> vetoedSimVertexes;

    for(size_t simTrackIndex = 0; simTrackIndex != simTracks_.size(); ++simTrackIndex) {
        // Check if the simTrack is excluded (includes non traceable and recovered by history)
        if(vetoedTracks.find(simTrackIndex) != vetoedTracks.end()) continue;

        SimTrack const& simTrack = simTracks_.at(simTrackIndex);

        TrackingParticle trackingParticle;

        // Set a bare tp (only with the psimhit) with a given simtrack
        // the function return true if it is tracable
        if(setTrackingParticle(simTrack, trackIdToHits, trackingParticle)) {
            // Follows the path upward recovering the history of the particle
            SimTrack const* currentSimTrack = &simTrack;

            // Initial condition for the tp and tv indexes
            int trackingParticleIndex = -1;
            int trackingVertexIndex = -1;

            do {
                // Set a new tracking particle for the current simtrack
                // and add it to the list of parent tracks of previous vertex
                if(trackingParticleIndex >= 0) {
                    setTrackingParticle(*currentSimTrack, trackIdToHits, trackingParticle);

                    // Set the tp index to its new value
                    trackingParticleIndex = trackingParticles.size();
                    // Push the tp in to the collection
                    trackingParticles.push_back(trackingParticle);

                    // Add the previous track to the list of decay vertexes of the new tp
                    trackingParticles.at(trackingParticleIndex).addDecayVertex(
                        TrackingVertexRef(refTrackingVertexes, trackingVertexIndex)
                    );

                    // Add the new tp to the list of parent tracks of the previous tv
                    trackingVertexes.at(trackingVertexIndex).addParentTrack(
                        TrackingParticleRef(refTrackingParticles, trackingParticleIndex)
                    );
                } else {
                    // Set the tp index to its new value
                    trackingParticleIndex = trackingParticles.size();
                    // Push the tp in to the collection
                    trackingParticles.push_back(trackingParticle);
                    // Vetoed the simTrack
                    vetoedTracks.insert(std::make_pair(simTrackIndex, trackingParticleIndex));
                }

                // Verify if the parent simVertex has a simTrack or if the source is a vetoSimVertex
                if(currentSimTrack->noVertex()) break;

                // Get the simTrack parent index (it is implicit should be in the same event as current)
                unsigned int parentSimVertexIndex = vertexIdToIndex[
                                                        EncodedTruthId(
                                                            currentSimTrack->eventId(),
                                                            currentSimTrack->vertIndex()
                                                        )
                                                    ];
                // Create a new tv
                TrackingVertex trackingVertex;
                // Get the parent simVertex associated to the current simTrack
                SimVertex const* parentSimVertex = &simVertexes_.at(parentSimVertexIndex);

                bool vetoSimVertex = vetoedSimVertexes.find(parentSimVertexIndex) != vetoedSimVertexes.end();

                // Check for a already visited parent simTrack
                if(!vetoSimVertex) {
                    // Set the tv by using simvertex
                    trackingVertexIndex = setTrackingVertex(*parentSimVertex, eventIdCounter, trackingVertexes, trackingVertex);

                    // Check if a new vertex needs to be created
                    if(trackingVertexIndex < 0) {
                        // Set the tv index ot its new value
                        trackingVertexIndex = trackingVertexes.size();
                        // Push the new tv in to the collection
                        trackingVertexes.push_back(trackingVertex);
                    } else {
                        // Get the postion and time of the vertex
                        LorentzVector const& position = trackingVertexes.at(trackingVertexIndex).position();
                        Vector xyz = Vector(position.x(), position.y(), position.z());
                        double t = position.t();
                        // Set the vertex postion of the tp to the closest vertex
                        trackingParticles.at(trackingParticleIndex).setVertex(xyz, t);
                    }

                    vetoedSimVertexes.insert(std::make_pair(parentSimVertexIndex, trackingVertexIndex));
                } else {
                    trackingVertexIndex = vetoedSimVertexes[parentSimVertexIndex];
                }

                // Set the newly created tv as parent vertex
                trackingParticles.at(trackingParticleIndex).setParentVertex(
                    TrackingVertexRef(refTrackingVertexes, trackingVertexIndex)
                );

                // Add the newly created tp to the tv daughter list
                trackingVertexes.at(trackingVertexIndex).addDaughterTrack(
                    TrackingParticleRef(refTrackingParticles, trackingParticleIndex)
                );

                // Verify if the parent simVertex has a simTrack or if the source is a vetoSimVertex
                if(parentSimVertex->noParent() || vetoSimVertex) break;

                // Get the next simTrack index (it is implicit should be in the same event as current).
                unsigned int nextSimTrackIndex = trackIdToIndex[
                                                     EncodedTruthId(
                                                         currentSimTrack->eventId(),
                                                         parentSimVertex->parentIndex()
                                                     )
                                                 ];

                // Check if the next track exist
                if(vetoedTracks.find(nextSimTrackIndex) != vetoedTracks.end()) {
                    // Add to the newly created tv the existent next simtrack in to parent list.
                    trackingVertexes.at(trackingVertexIndex).addParentTrack(
                        TrackingParticleRef(refTrackingParticles, vetoedTracks[nextSimTrackIndex])
                    );
                    // Add the vertex to list of decay vertexes of the new tp
                    trackingParticles.at(vetoedTracks[nextSimTrackIndex]).addDecayVertex(
                        TrackingVertexRef(refTrackingVertexes, trackingVertexIndex)
                    );
                    break;
                }

                // Vetoed the next simTrack
                vetoedTracks.insert(std::make_pair(nextSimTrackIndex, trackingParticleIndex));

                // Set the current simTrack as the next simTrack
                currentSimTrack = &simTracks_.at(nextSimTrackIndex);
            }
            while (!currentSimTrack->noVertex());
        }
    }
}


bool TrackingTruthProducer::setTrackingParticle(
    SimTrack const& simTrack,
    EncodedTruthIdToIndexes const& trackIdToHits,
    TrackingParticle& trackingParticle) {

    // Get the eventid associated to the track
    EncodedEventId trackEventId = simTrack.eventId();
    // Get the simtrack id
    EncodedTruthId simTrackId = EncodedTruthId(trackEventId, simTrack.trackId());

    // Location of the parent vertex
    LorentzVector position;
    // If not parent then location is (0,0,0,0)
    if(simTrack.noVertex()) {
        position = LorentzVector(0, 0, 0, 0);
    } else {
        position = simVertexes_.at(simTrack.vertIndex()).position();
    }

    // Define the default status and pdgid
    int status = -99;
    int pdgId = simTrack.type();

    int genParticleIndex = simTrack.genpartIndex();
    bool signalEvent = (trackEventId.event() == 0 && trackEventId.bunchCrossing() == 0);

    // In the case of a existing generated particle and track
    // event is signal redefine status a pdgId

    edm::Handle<edm::HepMCProduct> hepmc;

    if(genParticleIndex >= 0 && (signalEvent || useMultipleHepMCLabels_)) {
        // Get the generated particle
        hepmc = (useMultipleHepMCLabels_) ? hepMCProducts_.at(trackEventId.rawId()) : hepmc = hepMCProducts_.at(0);

        HepMC::GenParticle const* genParticle = hepmc->GetEvent()->barcode_to_particle(genParticleIndex);

        if(genParticle) {
            status = genParticle->status();
            pdgId  = genParticle->pdg_id();
        }
    }

    // Create a tp from the simtrack
    trackingParticle = TrackingParticle(
                           (char) simTrack.charge(),
                           simTrack.momentum(),
                           Vector(position.x(), position.y(), position.z()),
                           position.t(),
                           pdgId,
                           status,
                           trackEventId
                       );

    bool init = true;

    int processType = 0;
    int particleType = 0;

    // Counting the TP hits using the layers (as in ORCA).
    // Does seem to find less hits. maybe b/c layer is a number now, not a pointer
    int totalSimHits = 0;
    int oldLayer = 0;
    int newLayer = 0;
    int oldDetector = 0;
    int newDetector = 0;

    // Loop over the associated hits per track
    for(EncodedTruthIdToIndexes::const_iterator iEntry = trackIdToHits.lower_bound(simTrackId);
          iEntry != trackIdToHits.upper_bound(simTrackId);
          ++iEntry) {
        // Get a constant reference to the simhit
        PSimHit const& pSimHit = pSimHits_.at(iEntry->second);

        // Initial condition for consistent simhit selection
        if(init) {
            processType = pSimHit.processType();
            particleType = pSimHit.particleType();
            init = false;
        }

        // Check for delta and interaction products discards
        if(processType == pSimHit.processType() && particleType == pSimHit.particleType() && pdgId == pSimHit.particleType()) {
            trackingParticle.addPSimHit(pSimHit);

            unsigned int detectorIdIndex = pSimHit.detUnitId();
            DetId detectorId = DetId(detectorIdIndex);
            oldLayer = newLayer;
            oldDetector = newDetector;
            newLayer = LayerFromDetid(detectorIdIndex);
            newDetector = detectorId.subdetId();

            // Count hits using layers for glued detectors
            // newlayer !=0 excludes Muon layers set to 0 by LayerFromDetid
            if((oldLayer != newLayer || (oldLayer==newLayer && oldDetector!=newDetector)) && newLayer != 0) ++totalSimHits;
        }
    }

    // Set the number of matched simhits
    trackingParticle.setMatchedHit(totalSimHits);

    // Add the simtrack associated to the tp
    trackingParticle.addG4Track(simTrack);

    // Add the generator information
    if(genParticleIndex >= 0 && (signalEvent || useMultipleHepMCLabels_)) {
        trackingParticle.addGenParticle(GenParticleRef(hepmc, genParticleIndex));
    }

    if(selectorFlag_) {
      return selector_(trackingParticle);
    }

    return true;
}


int TrackingTruthProducer::setTrackingVertex(
    SimVertex const& simVertex,
    EncodedEventIdToIndex& eventIdCounter,
    TrackingVertexCollection& trackingVertexes,
    TrackingVertex& trackingVertex) {
    LorentzVector const& position = simVertex.position();

    // Look for close by vertexes
    for(std::size_t trackingVertexIndex = 0; trackingVertexIndex < trackingVertexes.size(); ++trackingVertexIndex) {
        // Calculate the distance
        double distance = (position - trackingVertexes.at(trackingVertexIndex).position()).P();
        // If the distance is under a given cut return the trackingVertex index (vertex merging)
        if(distance <= distanceCut_) {
            // Add simvertex to the pre existent tv
            trackingVertexes.at(trackingVertexIndex).addG4Vertex(simVertex);
            // return tv index
            return trackingVertexIndex;
        }
    }

    // Get the event if from the simvertex
    EncodedEventId simVertexEventId = simVertex.eventId();

    // Initialize the event counter
    if(eventIdCounter.find(simVertexEventId) == eventIdCounter.end()) {
        eventIdCounter[simVertexEventId] = 0;
    }

    // Get the simVertex id
    EncodedTruthId simVertexId = EncodedTruthId(simVertexEventId, eventIdCounter[simVertexEventId]);

    // Calculate if the vertex is in the tracker volume (it needs to be review for other detectors)
    bool inVolume = (position.Pt() < volumeRadius_ && std::abs(position.z()) < volumeZ_); // In or out of Tracker

    // Initialize the new vertex
    trackingVertex = TrackingVertex(position, inVolume, simVertexId);

    // Find the the closest GenVertexes to the created tv
    addCloseGenVertexes(trackingVertex);

    // Add the g4 vertex to the tv
    trackingVertex.addG4Vertex(simVertex);

    // Initialize the event counter
    ++eventIdCounter[simVertexEventId];

    return -1;
}


void TrackingTruthProducer::addCloseGenVertexes(TrackingVertex& trackingVertex) {
    // Get the generated particle
    edm::Handle<edm::HepMCProduct> hepmc = (useMultipleHepMCLabels_) ? hepMCProducts_.at(trackingVertex.eventId().rawId()) : hepMCProducts_.at(0);
    HepMC::GenEvent const* genEvent = hepmc->GetEvent();

    // Get the postion of the tv
    Vector tvPosition(trackingVertex.position().x(), trackingVertex.position().y(), trackingVertex.position().z());

    // Find HepMC vertices, put them in a close TrackingVertex (this could conceivably add the same GenVertex to multiple TrackingVertices)
    for(HepMC::GenEvent::vertex_const_iterator iGenVertex = genEvent->vertices_begin();
        iGenVertex != genEvent->vertices_end();
        ++iGenVertex) {
        // Get the position of the genVertex
        HepMC::ThreeVector rawPosition = (*iGenVertex)->position();

        // Convert to cm
        Vector genPosition(rawPosition.x()/10.0, rawPosition.y()/10.0, rawPosition.z()/10.0);

        // Calculate the dis
        double distance = sqrt((tvPosition - genPosition).mag2());

        if(distance <= distanceCut_) {
            trackingVertex.addGenVertex(GenVertexRef(hepmc, (*iGenVertex)->barcode()));
        }
    }
}


int TrackingTruthProducer::LayerFromDetid(unsigned int const& detid) {
    DetId detId = DetId(detid);

    if(detId.det() != DetId::Tracker) return 0;

    int layerNumber=0;
    unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());

    if(subdetId == StripSubdetector::TIB) {
        TIBDetId tibid(detId.rawId());
        layerNumber = tibid.layer();
    } else if(subdetId == StripSubdetector::TOB) {
        TOBDetId tobid(detId.rawId());
        layerNumber = tobid.layer();
    } else if(subdetId == StripSubdetector::TID) {
        TIDDetId tidid(detId.rawId());
        layerNumber = tidid.wheel();
    } else if(subdetId ==  StripSubdetector::TEC) {
        TECDetId tecid(detId.rawId());
        layerNumber = tecid.wheel();
    } else if(subdetId ==  PixelSubdetector::PixelBarrel) {
        PXBDetId pxbid(detId.rawId());
        layerNumber = pxbid.layer();
    } else if(subdetId ==  PixelSubdetector::PixelEndcap) {
        PXFDetId pxfid(detId.rawId());
        layerNumber = pxfid.disk();
    } else {
        edm::LogVerbatim("TrackingTruthProducer") << "Unknown subdetid: " <<  subdetId;
    }

    return layerNumber;
}

DEFINE_DIGI_ACCUMULATOR(TrackingTruthProducer);
