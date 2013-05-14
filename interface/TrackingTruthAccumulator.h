#ifndef TrackingAnalysis_TrackingTruthAccumulator_h
#define TrackingAnalysis_TrackingTruthAccumulator_h

#include "SimGeneral/MixingModule/interface/DigiAccumulatorMixMod.h"
#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"
#include <memory> // required for std::auto_ptr
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


// Forward declarations
namespace edm
{
	class ParameterSet;
	class EDProducer;
	class Event;
	class EventSetup;
}
class PileUpEventPrincipal;
class PSimHit;



/** @brief Replacement for TrackingTruthProducer in the new pileup mixing setup.
 * @author Mark Grimes (mark.grimes@bristol.ac.uk)
 * @date 11/Oct/2012
 */
class TrackingTruthAccumulator : public DigiAccumulatorMixMod
{
public:
	explicit TrackingTruthAccumulator( const edm::ParameterSet& config, edm::EDProducer& mixMod );
private:
	virtual void initializeEvent( const edm::Event& event, const edm::EventSetup& setup );
	virtual void accumulate( const edm::Event& event, const edm::EventSetup& setup );
	virtual void accumulate( const PileUpEventPrincipal& event, const edm::EventSetup& setup );
	virtual void finalizeEvent( edm::Event& event, const edm::EventSetup& setup );

	/** @brief Both forms of accumulate() delegate to this templated method. */
	template<class T> void accumulateEvent( const T& event, const edm::EventSetup& setup );

	/** @brief Fills the supplied vector with pointers to the SimHits, checking for bad modules if required */
	template<class T> void fillSimHits( std::vector<const PSimHit*>& returnValue, const T& event, const edm::EventSetup& setup );

	const std::string messageCategory_; ///< The message category used to send messages to MessageLogger

	const double volumeRadius_;
	const double volumeZ_;
	const bool ignoreTracksOutsideVolume_;

	/** The maximum bunch crossing BEFORE the signal crossing to create TrackinParticles for. Use positive values. If set to zero no
	 * previous bunches are added and only in-time, signal and after bunches (defined by maximumSubsequentBunchCrossing_) are used.*/
	const unsigned int maximumPreviousBunchCrossing_;
	/** The maximum bunch crossing AFTER the signal crossing to create TrackinParticles for. E.g. if set to zero only
	 * uses the signal and in time pileup (and previous bunches defined by the maximumPreviousBunchCrossing_ parameter). */
	const unsigned int maximumSubsequentBunchCrossing_;
	/// If bremsstrahlung merging, whether to also add the unmerged collection to the event or not.
	const bool createUnmergedCollection_;
	const bool createMergedCollection_;
	/// Whether or not to add the full parentage of any TrackingParticle that is inserted in the collection.
	const bool addAncestors_;

	/// As of 11/Feb/2013 this option hasn't been implemented yet.
	const bool removeDeadModules_;
	const std::string simHitLabel_;
	edm::ParameterSet simHitCollectionConfig_;
	edm::InputTag genParticleLabel_;

	bool selectorFlag_;
	TrackingParticleSelector selector_;
	/// Uses the same config as selector_, but can be used to drop out early since selector_ requires the TrackingParticle to be created first.
	bool chargedOnly_;
	/// Uses the same config as selector_, but can be used to drop out early since selector_ requires the TrackingParticle to be created first.
	bool signalOnly_;

public:
	// These always go hand in hand, and I need to pass them around in the internal
	// functions, so I might as well package them up in a struct.
	struct OutputCollections
	{
		std::auto_ptr<TrackingParticleCollection> pTrackingParticles;
		std::auto_ptr<TrackingVertexCollection> pTrackingVertices;
		TrackingParticleRefProd refTrackingParticles;
		TrackingVertexRefProd refTrackingVertexes;
	};
private:
	OutputCollections unmergedOutput_;
	OutputCollections mergedOutput_;
};

#endif // end of "#ifndef TrackingAnalysis_TrackingTruthAccumulator_h"
