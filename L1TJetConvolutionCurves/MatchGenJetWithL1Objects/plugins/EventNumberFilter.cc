// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/EventNumberFilter
// Class:      EventNumberFilter
// 
/**\class EventNumberFilter EventNumberFilter.cc L1TJetConvolutionCurves/EventNumberFilter/plugins/EventNumberFilter.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 19 Jul 2017 15:38:49 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


//
// class declaration
//

class EventNumberFilter : public edm::stream::EDFilter<> {
  public:
    explicit EventNumberFilter(const edm::ParameterSet&);
    ~EventNumberFilter();
  private:
    bool filter(edm::Event&, const edm::EventSetup&) override;
    // ----------member data ---------------------------

    unsigned int _eventNumber;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
EventNumberFilter::EventNumberFilter(const edm::ParameterSet& iConfig):
_eventNumber(0)
{
}


EventNumberFilter::~EventNumberFilter()
{
}

// ------------ method called to produce the data  ------------
bool
EventNumberFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (this -> _eventNumber == 159) {
    this -> _eventNumber++;
    return true;
  }
  this -> _eventNumber++;
  return false; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventNumberFilter);
