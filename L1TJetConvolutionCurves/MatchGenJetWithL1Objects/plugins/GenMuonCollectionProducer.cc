// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/GenMuonCollectionProducer
// Class:      GenMuonCollectionProducer
// 
/**\class GenMuonCollectionProducer GenMuonCollectionProducer.cc L1TJetConvolutionCurves/GenMuonCollectionProducer/plugins/GenMuonCollectionProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Simone Bologna
//         Created:  Sun, 20 May 2018 16:06:34 GMT
//
//


// system include files
#include <memory>
#include <csignal>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

class GenMuonCollectionProducer : public edm::stream::EDProducer<> {
   public:
      explicit GenMuonCollectionProducer(const edm::ParameterSet&);
      ~GenMuonCollectionProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genParticleCollectionTag;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
GenMuonCollectionProducer::GenMuonCollectionProducer(const edm::ParameterSet& iConfig)
{
  this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  produces<std::vector<reco::GenParticle> >( "GenMuons" ).setBranchAlias("GenMuons");
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


GenMuonCollectionProducer::~GenMuonCollectionProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
  delete (this -> _genParticleCollectionTag);
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenMuonCollectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
  iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);

  std::unique_ptr< std::vector<reco::GenParticle> > genMuonCollection(new std::vector<reco::GenParticle>);
  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    //Looking for a gen muon
    if (abs(genParticleIterator->pdgId()) == 13){
      genMuonCollection -> push_back(*genParticleIterator);
    }
  }

  iEvent.put(std::move(genMuonCollection), "GenMuons");
  
  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenMuonCollectionProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenMuonCollectionProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenMuonCollectionProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenMuonCollectionProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenMuonCollectionProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenMuonCollectionProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenMuonCollectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMuonCollectionProducer);
