#ifndef TrackingAnalysis_TrackingTruthProducer_h
#define TrackingAnalysis_TrackingTruthProducer_h

#include <map>
#include <vector>

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixMod.h"
#include "SimGeneral/TrackingAnalysis/interface/EncodedTruthId.h"
#include "SimGeneral/TrackingAnalysis/interface/PSimHitSelector.h"

#include "Utilities/Timing/interface/TimingReport.h"
#include "Utilities/Timing/interface/TimerStack.h"

class TrackingTruthProducer : public DigiAccumulatorMixMod {

public:

    TrackingTruthProducer(edm::ParameterSet const&, edm::EDProducer& mixMod);

    virtual void initializeEvent(edm::Event const&, edm::EventSetup const&);
    virtual void accumulate(edm::Event const&, edm::EventSetup const&);
    virtual void accumulate(PileUpEventPrincipal const&, edm::EventSetup const&);
    virtual void finalizeEvent(edm::Event&, edm::EventSetup const&);
    virtual void initializeBunchCrossing(edm::Event const&, edm::EventSetup const&, int bunchCrossing);
    virtual void finalizeBunchCrossing(edm::Event&, edm::EventSetup const&, int bunchCrossing);

private:

    typedef std::vector<SimTrack>  SimTracks;
    typedef std::vector<SimVertex> SimVertexes;

    const int         minBunch_;
    const int         maxBunch_;
    const std::vector<std::string> dataLabels_;
    const bool        useMultipleHepMCLabels_;
    const double      distanceCut_;
    const double      volumeRadius_;
    const double      volumeZ_;
    const bool        mergedBremsstrahlung_;
    const bool        removeDeadModules_;
    const std::string simHitLabel_;

    const std::string MessageCategory_;

    const PSimHitSelector        pSimHitSelector_;

    const bool selectorFlag_;
    const TrackingParticleSelector selector_;

    // Related to production

    std::vector<edm::Handle<edm::HepMCProduct> > hepMCProducts_;

    PSimHitSelector::PSimHitCollection pSimHits_;

    SimTracks   simTracks_;
    SimVertexes simVertexes_;

    typedef std::map<EncodedEventId, unsigned int> EncodedEventIdToIndex;
    typedef std::map<EncodedTruthId, unsigned int> EncodedTruthIdToIndex;
    typedef std::multimap<EncodedTruthId, unsigned int> EncodedTruthIdToIndexes;

    int LayerFromDetid(unsigned int const&);

    void accumulateEvent(edm::Event const&, edm::EventSetup const&);

    void associator(
        std::vector<PSimHit> const&,
        EncodedTruthIdToIndexes&
    );

    void associator(
        SimTracks const&,
        EncodedTruthIdToIndex&
    );

    void associator(
        SimVertexes const&,
        EncodedTruthIdToIndex&
    );

    void mergeBremsstrahlung(edm::Event& event,
                             TrackingParticleCollection& trackingParticles,
                             TrackingVertexCollection& trackingVertexes,
                             TrackingParticleCollection& mergedTrackingParticles,
                             TrackingVertexCollection& mergedTrackingVertexes);

    bool isBremsstrahlungVertex(
        TrackingVertex const& vertex,
        TrackingParticleCollection const& tPC
    );

    void createTrackingTruth(edm::Event& event,
                             EncodedTruthIdToIndexes const& trackIdToHits,
                             EncodedTruthIdToIndex& trackIdToIndex,
                             EncodedTruthIdToIndex& vertexIdToIndex,
                             TrackingParticleCollection& trackingParticles,
                             TrackingVertexCollection& trackingVertexes);

    bool setTrackingParticle(
        SimTrack const&,
        EncodedTruthIdToIndexes const& trackIdToHits,
        TrackingParticle&
    );

    int setTrackingVertex(
        SimVertex const&,
        EncodedEventIdToIndex& eventIdCounter,
        TrackingVertexCollection& trackingVertexes,
        TrackingVertex&
    );

    void addCloseGenVertexes(TrackingVertex&);
};

#endif
