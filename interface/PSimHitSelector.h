#ifndef SimGeneral_TrackingAnalysis_PSimHitSelector_h
#define SimGeneral_TrackingAnalysis_PSimHitSelector_h

#include <map>
#include <string>
#include <vector>

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
}

//! PSimHitSelector class
class PSimHitSelector {

public:

    typedef std::vector<PSimHit> PSimHitCollection;

    PSimHitSelector(edm::ParameterSet const &);

    ~PSimHitSelector() {}

    //! Select the psimhits and add them to a PSimHitCollection
    void select(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const;
    void selectMuon(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const;
    void selectPixel(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const;
    void selectTracker(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const;

protected:

    typedef std::map<std::string, std::vector<std::string> > PSimHitCollectionMap;

    std::string const simHitLabel_;
    PSimHitCollectionMap pSimHitCollectionMap_;
};

#endif
