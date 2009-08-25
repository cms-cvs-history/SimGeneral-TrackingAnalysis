
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "SimGeneral/TrackingAnalysis/interface/PSimHitSelector.h"

PSimHitSelector::PSimHitSelector(edm::ParameterSet const & config)
{
    // Initilize psimhit collection discriminated by sub systems
    edm::ParameterSet pSimHitCollections = config.getParameter<edm::ParameterSet>("simHitCollections");

    std::vector<std::string> subdetectors( pSimHitCollections.getParameterNames() );

    for (size_t i = 0; i < subdetectors.size(); ++i)
    {
        pSimHitCollectionMap_.insert(
            std::pair<std::string, std::vector<std::string> >(
                subdetectors[i],
                pSimHitCollections.getParameter<std::vector<std::string> >(subdetectors[i])
            )
        );
    }
}


void PSimHitSelector::newEvent(edm::Event const & event, edm::EventSetup const & setup)
{
    // Clear all the psimhit in the list
    pSimHits_.clear();

    // Look for all psimhit collections
    PSimHitCollectionMap::const_iterator pSimHitCollections = pSimHitCollectionMap_.begin();

    std::vector<const CrossingFrame<PSimHit> *> cfPSimHitProductPointers;

    for (; pSimHitCollections != pSimHitCollectionMap_.end(); ++pSimHitCollections)
    {
        // Grab all the PSimHit from the different sencitive volumes
        edm::Handle<CrossingFrame<PSimHit> > cfPSimHits;

        // Collect the product pointers to the different psimhit collection
        for (std::size_t i = 0; i < pSimHitCollections->second.size(); ++i)
        {
            event.getByLabel("mix", pSimHitCollections->second[i], cfPSimHits);
            cfPSimHitProductPointers.push_back(cfPSimHits.product());
        }
    }

    // Create a mix collection from the different psimhit collections
    std::auto_ptr<MixCollection<PSimHit> > pSimHits(new MixCollection<PSimHit>(cfPSimHitProductPointers));

    // Select all psimhits
    for (MixCollection<PSimHit>::MixItr pSimHit = pSimHits->begin(); pSimHit != pSimHits->end(); ++pSimHit)
        pSimHits_.push_back(*pSimHit);
}
