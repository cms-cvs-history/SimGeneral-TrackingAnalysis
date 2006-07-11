#ifndef TrackingTruthTest_h
#define TrackingTruthTest_h
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TrackingTruthTest  : public edm::EDAnalyzer {
 public:

  explicit TrackingTruthTest(const edm::ParameterSet& conf);

  virtual ~TrackingTruthTest(){}

  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

 private:
  edm::ParameterSet conf_;

};

#endif
