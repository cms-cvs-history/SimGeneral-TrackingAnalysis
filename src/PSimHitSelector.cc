#include "SimGeneral/TrackingAnalysis/interface/PSimHitSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CondFormats/CSCObjects/interface/CSCBadChambers.h"
#include "CondFormats/DataRecord/interface/CSCBadChambersRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelQuality.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"

#include <string>

PSimHitSelector::PSimHitSelector(edm::ParameterSet const & config) :
  simHitLabel_(config.getParameter<std::string>("simHitLabel")) {
    // Initilize psimhit collection discriminated by sub systems
    edm::ParameterSet pSimHitCollections = config.getParameter<edm::ParameterSet>("simHitCollections");

    std::vector<std::string> subdetectors(pSimHitCollections.getParameterNames());

    for(auto const& subdetector : subdetectors) {
        pSimHitCollectionMap_.insert(
            std::pair<std::string, std::vector<std::string> >(
                subdetector,
                pSimHitCollections.getParameter<std::vector<std::string> >(subdetector)
            )
        );
    }
}

void PSimHitSelector::select(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const {
    // Look for all psimhit collections
    for(auto const& simHitCollections : pSimHitCollectionMap_) {
      // Grab all the PSimHits from the different sensitive volumes
      // Collect the product pointers to the different psimhit collections
      for(auto const& simHitCollection : simHitCollections.second) {
        edm::Handle<std::vector<PSimHit> > evPSimHits;
        event.getByLabel(edm::InputTag(simHitLabel_, simHitCollection), evPSimHits);
        // Select all psimhits
        for(auto const& simHit : *evPSimHits.product()) {
          simHits.push_back(simHit);
        }
      }
    }
}

void PSimHitSelector::selectMuon(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const {
    // Look for psimhit collection associated to the muon system
    PSimHitCollectionMap::const_iterator pSimHitCollections = pSimHitCollectionMap_.find("muon");

    // Check that there are psimhit collections defined for the tracker
    if(pSimHitCollections == pSimHitCollectionMap_.end()) return;

    // Get CSC Bad Chambers (ME4/2)
    edm::ESHandle<CSCBadChambers> cscBadChambers;
    setup.get<CSCBadChambersRcd>().get(cscBadChambers);

    // Grab all the PSimHits from the different sensitive volumes
    for(auto const& simHitCollection : pSimHitCollections->second) {
      edm::Handle<std::vector<PSimHit> > evPSimHits;
      event.getByLabel(edm::InputTag(simHitLabel_, simHitCollection), evPSimHits);
      // Select only psimhits from alive modules
      for(auto const& simHit : *evPSimHits.product()) {
        DetId dId = DetId(simHit.detUnitId());
        if(dId.det() == DetId::Muon && dId.subdetId() == MuonSubdetId::CSC) {
          if(!cscBadChambers->isInBadChamber(CSCDetId(dId))) {
            simHits.push_back(simHit);
          }
        } else {
            simHits.push_back(simHit);
        }
      }
    }
}

void PSimHitSelector::selectPixel(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const {
    // Look for psimhit collection associated o the tracker
    PSimHitCollectionMap::const_iterator pSimHitCollections = pSimHitCollectionMap_.find("pixel");

    // Check that there are psimhit collections defined for the tracker
    if(pSimHitCollections == pSimHitCollectionMap_.end()) return;

    // Accessing dead pixel modules from DB:
    edm::ESHandle<SiPixelQuality> siPixelBadModule;
    setup.get<SiPixelQualityRcd>().get(siPixelBadModule);

    // Reading the DB information
    std::vector<SiPixelQuality::disabledModuleType> badModules(siPixelBadModule->getBadComponentList());
    SiPixelQuality pixelQuality(badModules);

    // Grab all the PSimHits from the different sensitive volumes
    for(auto const& simHitCollection : pSimHitCollections->second) {
      edm::Handle<std::vector<PSimHit> > evPSimHits;
      event.getByLabel(edm::InputTag(simHitLabel_, simHitCollection), evPSimHits);
      // Select only psimhits from alive modules
      for(auto const& simHit : *evPSimHits.product()) {
        if(!pixelQuality.IsModuleBad(simHit.detUnitId())) {
            simHits.push_back(simHit);
        }
      }
    }
}

void PSimHitSelector::selectTracker(PSimHitCollection& simHits, edm::Event const& event, edm::EventSetup const& setup) const {
    // Look for psimhit collection associated with the tracker
    PSimHitCollectionMap::const_iterator pSimHitCollections = pSimHitCollectionMap_.find("tracker");

    // Check that there are psimhit collections defined for the tracker
    if (pSimHitCollections == pSimHitCollectionMap_.end()) return;

    // Setup the cabling mapping
    std::map<uint32_t, std::vector<int> > theDetIdList;
    edm::ESHandle<SiStripDetCabling> detCabling;
    setup.get<SiStripDetCablingRcd>().get( detCabling );
    detCabling->addConnected(theDetIdList);

    // Grab all the PSimHits from the different sensitive volumes
    for(auto const& simHitCollection : pSimHitCollections->second) {
      edm::Handle<std::vector<PSimHit> > evPSimHits;
      event.getByLabel(edm::InputTag(simHitLabel_, simHitCollection), evPSimHits);
      // Select only psimhits from alive modules
      for(auto const& simHit : *evPSimHits.product()) {
        if(theDetIdList.empty()) {
          simHits.push_back(simHit);
        } else {
          uint32_t tkid = simHit.detUnitId();
          if(theDetIdList.find(tkid) != theDetIdList.end()){
            simHits.push_back(simHit);
          }
        }
      }
    }
}
