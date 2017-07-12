// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/PtPlotter
// Class:      PtPlotter
// 
/**\class PtPlotter PtPlotter.cc L1TJetConvolutionCurves/PtPlotter/plugins/PtPlotter.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 08:08:38 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PtPlotter : public edm::one::EDAnalyzer</*edm::one::SharedResources*/>  {
  public:
    explicit PtPlotter(const edm::ParameterSet&);
    ~PtPlotter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    edm::EDGetTokenT< std::vector< reco::GenJet > >* _jetCollectionTag;
    int _iEv;
    // ----------member data ---------------------------
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
PtPlotter::PtPlotter(const edm::ParameterSet& iConfig):
_iEv(0)
{
  try {
    _jetCollectionTag =  new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("GenxJetTag")));
  } catch (std::exception const & ex) {
    std::cout << "EXCEPTION!!!" << std::endl;
  }
  //now do what ever initialization is needed
  //usesResource("TFileService");
}


PtPlotter::~PtPlotter()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PtPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle < std::vector<reco::GenJet> > genJetHandle;
  iEvent.getByToken(*(this -> _jetCollectionTag), genJetHandle);
  for (const reco::GenJet & lGenJet : (*genJetHandle)){
    std::cout << this -> _iEv << " " << lGenJet.pt() << std::endl;
  }
  this -> _iEv ++;
}


// ------------ method called once each job just before starting event loop  ------------
void 
PtPlotter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PtPlotter::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PtPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PtPlotter);
