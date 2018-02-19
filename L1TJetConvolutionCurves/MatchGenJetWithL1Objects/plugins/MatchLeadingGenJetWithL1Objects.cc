// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchLeadingGenJetWithL1Objects
// Class:      MatchLeadingGenJetWithL1Objects
// 
/**\class MatchLeadingGenJetWithL1Objects MatchLeadingGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchLeadingGenJetWithL1Objects/plugins/MatchLeadingGenJetWithL1Objects.cc

 Description: Matches a gen jet with a L1T object and creates a tree

 Implementation:
    Receives tags of the various L1T objects
*/
//
// Original Author:  Simone Bologna
//         Created:  Wed, 12 Jul 2017 14:21:56 GMT
//
//


// system include files
#include <memory>
#include <signal.h>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"

#include <utility>

struct Particle {
  unsigned int id;
  float pt, eta, phi;
};

struct TriggerObject {
  unsigned int id;
  float pt, eta, phi;
  int hwQual;
};

class MatchLeadingGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchLeadingGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchLeadingGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector< std::tuple < TriggerObject, Particle, float> >
    _matchGenJetWithL1Object(
      const reco::GenJet &,
      const edm::Handle<BXVector<T>> &,
      float = 0.25
    );

    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector < std::tuple < TriggerObject, Particle, float > > &
    );

    edm::EDGetTokenT<std::vector<reco::GenParticle>> *_genParticleCollectionTag;
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Jet > > *_l1tJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Tau > > *_l1tTauCollectionTag;
    //!TODO: It should splitted in two different objects in the Stage-2
    edm::EDGetTokenT< BXVector< l1t::EGamma > > *_l1tEGammaCollectionTag;
    //!TODO: Like this
    //edm::EDGetTokenT< BXVector< l1t::Electron > > _l1tElectronCollectionTag;
    //edm::EDGetTokenT< BXVector< l1t::Gamma > > _l1tGammaCollectionTag;

    TTree * _l1tMuonLeadingGenJetTree;
    TTree * _l1tJetLeadingGenJetTree;
    TTree * _l1tEGammaLeadingGenJetTree;
    TTree * _l1tTauLeadingGenJetTree;
    TTree * _leadingGenJetTree;

    Particle _genJetParticle;
    TriggerObject _l1tObjectParticle;
    int _hwQual;
    float _deltaR2;
};

MatchLeadingGenJetWithL1Objects::MatchLeadingGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TMuonTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tJetLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tEGammaLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TEGammaTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tTauLeadingGenJetTree = fs -> make<TTree>("matchedLeadingGenJetL1TTauTree", "TTree with generator-level jet / L1T Tau information");
  this -> _leadingGenJetTree = fs -> make<TTree>("leadingGenJetTree", "TTree with generator-level jet information");

  this -> _l1tMuonLeadingGenJetTree -> Branch("leadingGenJet_id", &(this -> _genJetParticle.id), "leadingGenJet_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("leadingGenJet_pt", &(this -> _genJetParticle.pt), "leadingGenJet_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("leadingGenJet_eta", &(this -> _genJetParticle.eta), "leadingGenJet_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("leadingGenJet_phi", &(this -> _genJetParticle.phi), "leadingGenJet_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonLeadingGenJetTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  this -> _l1tMuonLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tJetLeadingGenJetTree -> Branch("leadingGenJet_id", &(this -> _genJetParticle.id), "leadingGenJet_id/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("leadingGenJet_pt", &(this -> _genJetParticle.pt), "leadingGenJet_pt/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("leadingGenJet_eta", &(this -> _genJetParticle.eta), "leadingGenJet_eta/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("leadingGenJet_phi", &(this -> _genJetParticle.phi), "leadingGenJet_phi/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetLeadingGenJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");
  this -> _l1tJetLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tEGammaLeadingGenJetTree -> Branch("leadingGenJet_id", &(this -> _genJetParticle.id), "leadingGenJet_id/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("leadingGenJet_pt", &(this -> _genJetParticle.pt), "leadingGenJet_pt/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("leadingGenJet_eta", &(this -> _genJetParticle.eta), "leadingGenJet_eta/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("leadingGenJet_phi", &(this -> _genJetParticle.phi), "leadingGenJet_phi/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_id", &(this -> _l1tObjectParticle.id), "l1tEGamma_id/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectParticle.pt), "l1tEGamma_pt/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectParticle.eta), "l1tEGamma_eta/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectParticle.phi), "l1tEGamma_phi/F");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("l1tEGamma_qual", &(this -> _l1tObjectParticle.hwQual), "l1tEGamma_qual/i");
  this -> _l1tEGammaLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tTauLeadingGenJetTree -> Branch("leadingGenJet_id", &(this -> _genJetParticle.id), "leadingGenJet_id/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("leadingGenJet_pt", &(this -> _genJetParticle.pt), "leadingGenJet_pt/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("leadingGenJet_eta", &(this -> _genJetParticle.eta), "leadingGenJet_eta/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("leadingGenJet_phi", &(this -> _genJetParticle.phi), "leadingGenJet_phi/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_id", &(this -> _l1tObjectParticle.id), "l1tTau_id/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_pt", &(this -> _l1tObjectParticle.pt), "l1tTau_pt/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_eta", &(this -> _l1tObjectParticle.eta), "l1tTau_eta/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_phi", &(this -> _l1tObjectParticle.phi), "l1tTau_phi/F");
  this -> _l1tTauLeadingGenJetTree -> Branch("l1tTau_qual", &(this -> _l1tObjectParticle.hwQual), "l1tTau_qual/i");
  this -> _l1tTauLeadingGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _leadingGenJetTree -> Branch("leadingGenJet_id", &(this -> _genJetParticle.id), "leadingGenJet_id/i");
  this -> _leadingGenJetTree -> Branch("leadingGenJet_pt", &(this -> _genJetParticle.pt), "leadingGenJet_pt/F");
  this -> _leadingGenJetTree -> Branch("leadingGenJet_eta", &(this -> _genJetParticle.eta), "leadingGenJet_eta/F");
  this -> _leadingGenJetTree -> Branch("leadingGenJet_phi", &(this -> _genJetParticle.phi), "leadingGenJet_phi/F");

}

void MatchLeadingGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
  try
  {
    this -> _genParticleCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenParticle > >(consumes< std::vector< reco::GenParticle > > (iConfig.getParameter< edm::InputTag >("genParticleCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> genParticleCollectionTag not found." << std::endl;
    this -> _genParticleCollectionTag = NULL;
  }

  try
  {
    this -> _l1tMuonCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Muon > >(consumes< BXVector< l1t::Muon > > (iConfig.getParameter< edm::InputTag >("l1tMuonCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tMuonCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tMuonCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tJetCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Jet > >(consumes< BXVector< l1t::Jet > > (iConfig.getParameter< edm::InputTag >("l1tJetCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tJetCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tJetCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tEGammaCollectionTag = new edm::EDGetTokenT< BXVector < l1t::EGamma > >(consumes< BXVector< l1t::EGamma > > (iConfig.getParameter< edm::InputTag >("l1tEGammaCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tEGammaCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tEGammaCollectionTag = NULL;
  }
  
  try
  {
    this -> _l1tTauCollectionTag = new edm::EDGetTokenT< BXVector < l1t::Tau > >(consumes< BXVector< l1t::Tau > > (iConfig.getParameter< edm::InputTag >("l1tTauCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> l1tTauCollectionTag not found, proceeding without computing the corresponding convolution." << std::endl;
    this -> _l1tTauCollectionTag = NULL;
  }
  
  return;
}

void MatchLeadingGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _genParticleCollectionTag) delete this -> _genParticleCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  if (this -> _l1tJetCollectionTag) delete this -> _l1tJetCollectionTag;
  if (this -> _l1tEGammaCollectionTag) delete this -> _l1tEGammaCollectionTag;
  if (this -> _l1tTauCollectionTag) delete this -> _l1tTauCollectionTag;
}

MatchLeadingGenJetWithL1Objects::~MatchLeadingGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchLeadingGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  //Looking for the leading gen jet

  float maxPt = 0;
  bool save = false;
  const reco::GenJet* leadingGenJet = NULL;

  for (auto genJetIterator = genJetCollectionHandle->begin(); genJetIterator != genJetCollectionHandle->end(); genJetIterator++)
  {
    // Only barrel leading l1t jets
    if ((genJetIterator->pt() > maxPt) && (genJetIterator->eta() < 1.44) && (genJetIterator->eta() > -1.44))
    {
      save = true;
      maxPt = genJetIterator->pt();
      leadingGenJet = &(*genJetIterator);
      this->_genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
      this->_genJetParticle.pt = genJetIterator->pt();
      this->_genJetParticle.eta = genJetIterator->eta();
      this->_genJetParticle.phi = genJetIterator->phi();
    }
  }
  if (save)
    this->_leadingGenJetTree->Fill();

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes

  if (leadingGenJet != NULL){
    if (this -> _l1tMuonCollectionTag)
    {
      // Muon veto. For muon matching
      edm::Handle < std::vector< reco::GenParticle > > genParticleCollectionHandle;
      iEvent.getByToken(*(this -> _genParticleCollectionTag), genParticleCollectionHandle);
      bool hasMuons = false;
      for (auto genParticleIterator = genParticleCollectionHandle -> begin(); genParticleIterator != genParticleCollectionHandle -> end(); genParticleIterator++ )
      {
        if (abs(genParticleIterator->pdgId()) == 13)
          hasMuons = true;
      }

      if (!hasMuons)
      {
        edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
        iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
        auto l1tMuonGenJetPairs = this -> _matchGenJetWithL1Object<>(*leadingGenJet, l1tMuonCollectionHandle);
        this -> _fillTreeWithMatchedPairs(*(this -> _l1tMuonLeadingGenJetTree), l1tMuonGenJetPairs);     
      }
    }
    
    if (this -> _l1tJetCollectionTag)
    {
      edm::Handle < BXVector< l1t::Jet > > l1tJetCollectionHandle;
      iEvent.getByToken(*(this -> _l1tJetCollectionTag), l1tJetCollectionHandle);
      auto l1tJetGenJetPairs = this -> _matchGenJetWithL1Object<>(*leadingGenJet, l1tJetCollectionHandle);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tJetLeadingGenJetTree), l1tJetGenJetPairs);
    }

    if (this -> _l1tEGammaCollectionTag)
    {
      edm::Handle < BXVector< l1t::EGamma > > l1tEGammaCollectionHandle;
      iEvent.getByToken(*(this -> _l1tEGammaCollectionTag), l1tEGammaCollectionHandle);
      auto l1tEGammaGenJetPairs = this -> _matchGenJetWithL1Object<>(*leadingGenJet, l1tEGammaCollectionHandle);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tEGammaLeadingGenJetTree), l1tEGammaGenJetPairs);
    }
    
    if (this -> _l1tTauCollectionTag)
    {
      edm::Handle < BXVector< l1t::Tau > > l1tTauCollectionHandle;
      iEvent.getByToken(*(this -> _l1tTauCollectionTag), l1tTauCollectionHandle);
      auto l1tTauGenJetPairs = this -> _matchGenJetWithL1Object<>(*leadingGenJet, l1tTauCollectionHandle);
      this -> _fillTreeWithMatchedPairs(*(this -> _l1tTauLeadingGenJetTree), l1tTauGenJetPairs);
    }
  }

}

void
MatchLeadingGenJetWithL1Objects::_fillTreeWithMatchedPairs
(
  TTree & aTree,
  const std::vector < std::tuple <TriggerObject, Particle, float > > & matchedL1TObjectJetPairs
)
{
  for (const auto & matchTuple : matchedL1TObjectJetPairs)
  {
    const TriggerObject & l1tObject = std::get<0>(matchTuple);
    const Particle & genJet = std::get<1>(matchTuple);
    float deltaR2 = std::get<2>(matchTuple);
    this -> _genJetParticle.id = genJet.id;
    this -> _genJetParticle.pt = genJet.pt;
    this -> _genJetParticle.eta = genJet.eta;
    this -> _genJetParticle.phi = genJet.phi;
    this -> _l1tObjectParticle.id = l1tObject.id;
    this -> _l1tObjectParticle.pt = l1tObject.pt;
    this -> _l1tObjectParticle.eta = l1tObject.eta;
    this -> _l1tObjectParticle.phi = l1tObject.phi;
    this -> _l1tObjectParticle.hwQual = l1tObject.hwQual;
    this -> _deltaR2 = deltaR2;
    aTree.Fill();
  }
}

// Considers a genJet and looks for the closest l1tObject within that radius (1 genjet can have only one l1tobject)

template <class T> // <3
const std::vector <std::tuple < TriggerObject, Particle, float > > 
MatchLeadingGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const reco::GenJet & genJet,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle,
  float dr2Min
)
{

  std::vector< std::tuple < TriggerObject, Particle , float > > l1tObjectGenJetPairs;
  // for each object in the genJet collection we look for the closest l1tobject in a wide range
  
  bool foundMatch = false;
  std::tuple<TriggerObject, Particle, float > l1tObjectGenJetPair;
  
  Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
  matchedGenJet.id = 0;
  matchedGenJet.pt = genJet.pt();
  matchedGenJet.phi = genJet.phi();
  matchedGenJet.eta = genJet.eta();

  for (
        typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle->begin(0);
        bx0Iterator != l1tObjectCollectionHandle->end(0);
        bx0Iterator++
  )
  {
    float dr2 = reco::deltaR2(*bx0Iterator, genJet);

    if ((dr2 < dr2Min))
    {
      TriggerObject &matchedL1TObject = std::get<0>(l1tObjectGenJetPair);
      matchedL1TObject.id = (bx0Iterator - l1tObjectCollectionHandle->begin(0));
      matchedL1TObject.pt = bx0Iterator->pt();
      matchedL1TObject.phi = bx0Iterator->phi();
      matchedL1TObject.eta = bx0Iterator->eta();
      matchedL1TObject.hwQual = bx0Iterator->hwQual();

      std::get<2>(l1tObjectGenJetPair) = dr2;
      dr2Min = dr2;
      foundMatch = true;
    }
  
  } 
  //if we have found a compatible l1t object we push the object-jet pair in a vector which will be our result
  if (foundMatch)
    l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchLeadingGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchLeadingGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchLeadingGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchLeadingGenJetWithL1Objects);
