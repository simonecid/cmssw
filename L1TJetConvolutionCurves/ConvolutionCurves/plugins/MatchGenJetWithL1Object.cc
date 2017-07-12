// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchGenJetWithL1Objects
// Class:      MatchGenJetWithL1Objects
// 
/**\class MatchGenJetWithL1Objects MatchGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchGenJetWithL1Objects.cc

 Description: Matches a gen jet with a L1T object and creates a tree

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 14:21:56 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MatchGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MatchGenJetWithL1Objects(const edm::ParameterSet&);
      ~MatchGenJetWithL1Objects();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  
      edm::EDGetTokenT< std::vector< reco::GenJet > > _genJetCollectionTag;
      edm::EDGetTokenT< std::vector< l1t::Muon > > _l1tMuonCollectionTag;
      edm::EDGetTokenT< std::vector< l1t::Jet > > _l1tJetCollectionTag;
      //!TODO: It should splitted in two different objects in the Stage-2
      edm::EDGetTokenT< std::vector< l1t::EGamma > > _l1tEGammaCollectionTag;
      //!TODO: Like this
      //edm::EDGetTokenT< std::vector< l1t::Electron > > _l1tElectronCollectionTag;
      //edm::EDGetTokenT< std::vector< l1t::Gamma > > _l1tGammaCollectionTag;
      edm::EDGetTokenT< std::vector< l1t::Tau > > _l1tTauCollectionTag;
      
      // ----------member data ---------------------------
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig)

{
   //usesResource("TFileService");
}


MatchGenJetWithL1Objects::~MatchGenJetWithL1Objects()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
MatchGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MatchGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MatchGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchGenJetWithL1Objects);
