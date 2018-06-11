// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/PropagateGenMuonAndGenJetsTo2ndMuonStation
// Class:      PropagateGenMuonAndGenJetsTo2ndMuonStation
// 
/**\class PropagateGenMuonAndGenJetsTo2ndMuonStation PropagateGenMuonAndGenJetsTo2ndMuonStation.cc L1TJetConvolutionCurves/PropagateGenMuonAndGenJetsTo2ndMuonStation/plugins/PropagateGenMuonAndGenJetsTo2ndMuonStation.cc

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
#include "DataFormats/JetReco/interface/GenJet.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//
// class declaration
//

class PropagateGenMuonAndGenJetsTo2ndMuonStation : public edm::stream::EDProducer<> {
   public:
      explicit PropagateGenMuonAndGenJetsTo2ndMuonStation(const edm::ParameterSet&);
      ~PropagateGenMuonAndGenJetsTo2ndMuonStation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genParticleCollectionTag;
      edm::EDGetTokenT<std::vector<reco::GenJet>> *_genJetCollectionTag;

      PropagateToMuon _propagator2nd;

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
PropagateGenMuonAndGenJetsTo2ndMuonStation::PropagateGenMuonAndGenJetsTo2ndMuonStation(const edm::ParameterSet& iConfig) :
_propagator2nd(iConfig.getParameter<edm::ParameterSet>("prop2nd"))
{
  this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  produces<std::vector<reco::GenParticle> >( "PropagatedGenMuons" ).setBranchAlias("PropagatedGenMuons");
  produces<std::vector<reco::GenJet> >( "PropagatedGenJets" ).setBranchAlias("PropagatedGenJets");
  
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


PropagateGenMuonAndGenJetsTo2ndMuonStation::~PropagateGenMuonAndGenJetsTo2ndMuonStation()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)
  delete (this -> _genParticleCollectionTag);
  delete (this -> _genJetCollectionTag);
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  this -> _propagator2nd.init(iSetup);

  edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
  iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);

  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  std::unique_ptr< std::vector<reco::GenParticle> > genMuonCollection(new std::vector<reco::GenParticle>);
  std::unique_ptr< std::vector<reco::GenJet> > genJetCollection(new std::vector<reco::GenJet>);

  for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
  {
    //Looking for a gen muon
    if (abs(genParticleIterator->pdgId()) == 13){
      //Propagating it to the muon chamber
      reco::GenParticle* propagatedMuon = genParticleIterator -> clone();
      TrajectoryStateOnSurface stateAt2ndMuonStation = this -> _propagator2nd.extrapolate(*propagatedMuon);
      if (stateAt2ndMuonStation.isValid()){
        math::PtEtaPhiMLorentzVector momentum(propagatedMuon -> polarP4());
        momentum.SetPhi(stateAt2ndMuonStation.globalPosition().phi());
        momentum.SetEta(stateAt2ndMuonStation.globalPosition().eta());
        propagatedMuon->setP4(momentum);
      }
      genMuonCollection -> push_back(*propagatedMuon);
    }
  }

  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    reco::GenJet* propagatedJet = genJetIterator -> clone();
    TrajectoryStateOnSurface stateAt2ndMuonStation = this -> _propagator2nd.extrapolate(*propagatedJet);
    if (stateAt2ndMuonStation.isValid()){
      math::PtEtaPhiMLorentzVector momentum(propagatedJet -> polarP4());      
      momentum.SetPhi(stateAt2ndMuonStation.globalPosition().phi());
      momentum.SetEta(stateAt2ndMuonStation.globalPosition().eta());
      propagatedJet->setP4(momentum);      
    }
    genJetCollection -> push_back(*propagatedJet);
  }

  iEvent.put(std::move(genMuonCollection), "PropagatedGenMuons");
  iEvent.put(std::move(genJetCollection), "PropagatedGenJets");
  
  return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PropagateGenMuonAndGenJetsTo2ndMuonStation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PropagateGenMuonAndGenJetsTo2ndMuonStation);
