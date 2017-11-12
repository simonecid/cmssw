// -*- C++ -*-
//
// Package:    L1TJetConvolutionCurves/MatchGenJetWithL1Objects
// Class:      MatchGenJetWithL1Objects
// 
/**\class MatchGenJetWithL1Objects MatchGenJetWithL1Objects.cc L1TJetConvolutionCurves/MatchGenJetWithL1Objects/plugins/MatchGenJetWithL1Objects.cc

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

class MatchGenJetWithL1Objects : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit MatchGenJetWithL1Objects(const edm::ParameterSet&);
    ~MatchGenJetWithL1Objects();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void _getTokens(const edm::ParameterSet&);
    void _freeTokens();

    // Beautiful template function to perform matching between a gen jet and whatever L1T object
    template <class T> // <3
    const std::vector < std::tuple < TriggerObject, Particle, float > > 
    _matchGenJetWithL1Object
    (
      const edm::Handle< std::vector< reco::GenJet > > &,
      const edm::Handle < BXVector < T > > & 
    );
    
    void
    _fillTreeWithMatchedPairs
    (
      TTree &,
      const std::vector < std::tuple < TriggerObject, Particle, float > > &
    );
  
    edm::EDGetTokenT< std::vector< reco::GenJet > > *_genJetCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Muon > > *_l1tMuonCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Jet > > *_l1tJetCollectionTag;
    //!TODO: It should splitted in two different objects in the Stage-2
    edm::EDGetTokenT< BXVector< l1t::EGamma > > *_l1tEGammaCollectionTag;
    //!TODO: Like this
    //edm::EDGetTokenT< BXVector< l1t::Electron > > _l1tElectronCollectionTag;
    //edm::EDGetTokenT< BXVector< l1t::Gamma > > _l1tGammaCollectionTag;
    edm::EDGetTokenT< BXVector< l1t::Tau > > *_l1tTauCollectionTag;

    TTree * _l1tMuonGenJetTree;
    TTree * _l1tMuonTree;
    TTree * _l1tLeadingMuonTree;
    TTree * _l1tJetGenJetTree;
    TTree * _l1tJetTree;
    TTree * _l1tLeadingJetTree;
    TTree * _l1tEGammaGenJetTree;
    TTree * _l1tEGammaTree;
    TTree * _l1tLeadingEGammaTree;
    TTree * _l1tTauGenJetTree;
    TTree * _l1tTauTree;
    TTree * _l1tLeadingTauTree;
    TTree * _genJetTree;

    Particle _genJetParticle;
    TriggerObject _l1tObjectParticle;
    int _hwQual;
    float _deltaR2;
};

MatchGenJetWithL1Objects::MatchGenJetWithL1Objects(const edm::ParameterSet& iConfig)
{  
  this -> _getTokens(iConfig);
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  this -> _l1tMuonGenJetTree = fs -> make<TTree>("matchedL1TMuonGenJetTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tJetGenJetTree = fs -> make<TTree>("matchedL1TJetGenJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tEGammaGenJetTree = fs -> make<TTree>("matchedL1TEGammaGenJetTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tTauGenJetTree = fs -> make<TTree>("matchedL1TTauGenJetTree", "TTree with generator-level jet / L1T Tau information");
  this -> _l1tMuonTree = fs -> make<TTree>("l1tMuonTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tLeadingMuonTree = fs -> make<TTree>("l1tLeadingMuonTree", "TTree with generator-level jet / L1T Muon information");
  this -> _l1tJetTree = fs -> make<TTree>("l1tJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tLeadingJetTree = fs -> make<TTree>("l1tLeadingJetTree", "TTree with generator-level jet / L1T Jet information");
  this -> _l1tEGammaTree = fs -> make<TTree>("l1tEGammaTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tLeadingEGammaTree = fs -> make<TTree>("l1tLeadingEGammaTree", "TTree with generator-level jet / L1T EGamma information");
  this -> _l1tTauTree = fs -> make<TTree>("l1tTauTree", "TTree with generator-level jet / L1T Tau information");
  this -> _l1tLeadingTauTree = fs -> make<TTree>("l1tLeadingTauTree", "TTree with generator-level jet / L1T Tau information");
  this -> _genJetTree = fs -> make<TTree>("genJetTree", "TTree with generator-level jet information");

  this -> _l1tMuonGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tMuonGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tMuonGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tMuonGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonGenJetTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  this -> _l1tMuonGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tJetGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tJetGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetGenJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");
  this -> _l1tJetGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tEGammaGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tEGammaGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_id", &(this -> _l1tObjectParticle.id), "l1tEGamma_id/i");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectParticle.pt), "l1tEGamma_pt/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectParticle.eta), "l1tEGamma_eta/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectParticle.phi), "l1tEGamma_phi/F");
  this -> _l1tEGammaGenJetTree -> Branch("l1tEGamma_qual", &(this -> _l1tObjectParticle.hwQual), "l1tEGamma_qual/i");
  this -> _l1tEGammaGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  this -> _l1tTauGenJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _l1tTauGenJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _l1tTauGenJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _l1tTauGenJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_id", &(this -> _l1tObjectParticle.id), "l1tTau_id/i");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_pt", &(this -> _l1tObjectParticle.pt), "l1tTau_pt/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_eta", &(this -> _l1tObjectParticle.eta), "l1tTau_eta/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_phi", &(this -> _l1tObjectParticle.phi), "l1tTau_phi/F");
  this -> _l1tTauGenJetTree -> Branch("l1tTau_qual", &(this -> _l1tObjectParticle.hwQual), "l1tTau_qual/i");
  this -> _l1tTauGenJetTree -> Branch("deltaR2", &(this -> _deltaR2), "deltaR2/F");

  //Used to detemine the prob that a jet will be misidentified binned in pt
  this -> _genJetTree -> Branch("genJet_id", &(this -> _genJetParticle.id), "genJet_id/i");
  this -> _genJetTree -> Branch("genJet_pt", &(this -> _genJetParticle.pt), "genJet_pt/F");
  this -> _genJetTree -> Branch("genJet_eta", &(this -> _genJetParticle.eta), "genJet_eta/F");
  this -> _genJetTree -> Branch("genJet_phi", &(this -> _genJetParticle.phi), "genJet_phi/F");

  this -> _l1tEGammaTree -> Branch("l1tEGamma_id", &(this -> _l1tObjectParticle.id), "l1tEGamma_id/i");
  this -> _l1tEGammaTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectParticle.pt), "l1tEGamma_pt/F");
  this -> _l1tEGammaTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectParticle.eta), "l1tEGamma_eta/F");
  this -> _l1tEGammaTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectParticle.phi), "l1tEGamma_phi/F");
  this -> _l1tEGammaTree -> Branch("l1tEGamma_qual", &(this -> _l1tObjectParticle.hwQual), "l1tEGamma_qual/i");

  this -> _l1tTauTree -> Branch("l1tTau_id", &(this -> _l1tObjectParticle.id), "l1tTau_id/i");
  this -> _l1tTauTree -> Branch("l1tTau_pt", &(this -> _l1tObjectParticle.pt), "l1tTau_pt/F");
  this -> _l1tTauTree -> Branch("l1tTau_eta", &(this -> _l1tObjectParticle.eta), "l1tTau_eta/F");
  this -> _l1tTauTree -> Branch("l1tTau_phi", &(this -> _l1tObjectParticle.phi), "l1tTau_phi/F");
  this -> _l1tTauTree -> Branch("l1tTau_qual", &(this -> _l1tObjectParticle.hwQual), "l1tTau_qual/i");

  this -> _l1tMuonTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tMuonTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tMuonTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tMuonTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tMuonTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");
  

  this -> _l1tJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");

  this -> _l1tLeadingEGammaTree -> Branch("l1tEGamma_id", &(this -> _l1tObjectParticle.id), "l1tEGamma_id/i");
  this -> _l1tLeadingEGammaTree -> Branch("l1tEGamma_pt", &(this -> _l1tObjectParticle.pt), "l1tEGamma_pt/F");
  this -> _l1tLeadingEGammaTree -> Branch("l1tEGamma_eta", &(this -> _l1tObjectParticle.eta), "l1tEGamma_eta/F");
  this -> _l1tLeadingEGammaTree -> Branch("l1tEGamma_phi", &(this -> _l1tObjectParticle.phi), "l1tEGamma_phi/F");
  this -> _l1tLeadingEGammaTree -> Branch("l1tEGamma_qual", &(this -> _l1tObjectParticle.hwQual), "l1tEGamma_qual/i");

  this -> _l1tLeadingTauTree -> Branch("l1tTau_id", &(this -> _l1tObjectParticle.id), "l1tTau_id/i");
  this -> _l1tLeadingTauTree -> Branch("l1tTau_pt", &(this -> _l1tObjectParticle.pt), "l1tTau_pt/F");
  this -> _l1tLeadingTauTree -> Branch("l1tTau_eta", &(this -> _l1tObjectParticle.eta), "l1tTau_eta/F");
  this -> _l1tLeadingTauTree -> Branch("l1tTau_phi", &(this -> _l1tObjectParticle.phi), "l1tTau_phi/F");
  this -> _l1tLeadingTauTree -> Branch("l1tTau_qual", &(this -> _l1tObjectParticle.hwQual), "l1tTau_qual/i");

  this -> _l1tLeadingMuonTree -> Branch("l1tMuon_id", &(this -> _l1tObjectParticle.id), "l1tMuon_id/i");
  this -> _l1tLeadingMuonTree -> Branch("l1tMuon_pt", &(this -> _l1tObjectParticle.pt), "l1tMuon_pt/F");
  this -> _l1tLeadingMuonTree -> Branch("l1tMuon_eta", &(this -> _l1tObjectParticle.eta), "l1tMuon_eta/F");
  this -> _l1tLeadingMuonTree -> Branch("l1tMuon_phi", &(this -> _l1tObjectParticle.phi), "l1tMuon_phi/F");
  this -> _l1tLeadingMuonTree -> Branch("l1tMuon_qual", &(this -> _l1tObjectParticle.hwQual), "l1tMuon_qual/i");  

  this -> _l1tLeadingJetTree -> Branch("l1tJet_id", &(this -> _l1tObjectParticle.id), "l1tJet_id/i");
  this -> _l1tLeadingJetTree -> Branch("l1tJet_pt", &(this -> _l1tObjectParticle.pt), "l1tJet_pt/F");
  this -> _l1tLeadingJetTree -> Branch("l1tJet_eta", &(this -> _l1tObjectParticle.eta), "l1tJet_eta/F");
  this -> _l1tLeadingJetTree -> Branch("l1tJet_phi", &(this -> _l1tObjectParticle.phi), "l1tJet_phi/F");
  this -> _l1tLeadingJetTree -> Branch("l1tJet_qual", &(this -> _l1tObjectParticle.hwQual), "l1tJet_qual/i");

}

void MatchGenJetWithL1Objects::_getTokens(const edm::ParameterSet& iConfig)
{
  
  this -> _genJetCollectionTag = new edm::EDGetTokenT< std::vector< reco::GenJet > >(consumes< std::vector< reco::GenJet > > (iConfig.getParameter< edm::InputTag >("genJetCollectionTag")));
  // Taking the tag of the various L1T object collections
  // If a parameter is omitted that object will not be studied
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

void MatchGenJetWithL1Objects::_freeTokens()
{
  if (this -> _genJetCollectionTag) delete this -> _genJetCollectionTag;
  if (this -> _l1tMuonCollectionTag) delete this -> _l1tMuonCollectionTag;
  if (this -> _l1tJetCollectionTag) delete this -> _l1tJetCollectionTag;
  if (this -> _l1tEGammaCollectionTag) delete this -> _l1tEGammaCollectionTag;
  if (this -> _l1tTauCollectionTag) delete this -> _l1tTauCollectionTag;
}

MatchGenJetWithL1Objects::~MatchGenJetWithL1Objects()
{
  this -> _freeTokens();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MatchGenJetWithL1Objects::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //Retrieving gen and l1t stuff
  edm::Handle < std::vector< reco::GenJet > > genJetCollectionHandle;
  iEvent.getByToken(*(this -> _genJetCollectionTag), genJetCollectionHandle);

  // I want to save for each event the highest momentum l1t(Muon/EGamma/Tau/Jet) for performance purposes
  
  if (this -> _l1tMuonCollectionTag)
  {
    edm::Handle < BXVector< l1t::Muon > > l1tMuonCollectionHandle;
    iEvent.getByToken(*(this -> _l1tMuonCollectionTag), l1tMuonCollectionHandle);
    auto l1tMuonGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tMuonCollectionHandle);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tMuonGenJetTree), l1tMuonGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tMuonIterator = l1tMuonCollectionHandle -> begin(0); l1tMuonIterator != l1tMuonCollectionHandle -> end(0); l1tMuonIterator++ )
    {
      if (l1tMuonIterator -> hwQual() < 8) continue;
      if (l1tMuonIterator -> pt() > maxPt)
      {
        save = true;
        maxPt = l1tMuonIterator -> pt();
        this -> _l1tObjectParticle.hwQual = l1tMuonIterator -> hwQual();      
        this -> _l1tObjectParticle.id = (l1tMuonIterator - l1tMuonCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tMuonIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tMuonIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tMuonIterator -> phi();
      }
    }
    if (save) {
      //std::cout << "I am saving a muon object w/ qual " << this -> _l1tObjectParticle.hwQual << std::endl;      
      this -> _l1tLeadingMuonTree -> Fill();
    }
    
    for (auto l1tMuonIterator = l1tMuonCollectionHandle -> begin(0); l1tMuonIterator != l1tMuonCollectionHandle -> end(0); l1tMuonIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tMuonIterator - l1tMuonCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tMuonIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tMuonIterator -> eta();
      this -> _l1tObjectParticle.phi = l1tMuonIterator -> phi();
      this -> _l1tObjectParticle.hwQual = l1tMuonIterator -> hwQual();      
      this -> _l1tMuonTree -> Fill();
    }
    
  }
  
  if (this -> _l1tJetCollectionTag)
  {
    edm::Handle < BXVector< l1t::Jet > > l1tJetCollectionHandle;
    iEvent.getByToken(*(this -> _l1tJetCollectionTag), l1tJetCollectionHandle);
    auto l1tJetGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tJetCollectionHandle);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tJetGenJetTree), l1tJetGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tJetIterator = l1tJetCollectionHandle -> begin(0); l1tJetIterator != l1tJetCollectionHandle -> end(0); l1tJetIterator++ )
    {
      if ((l1tJetIterator->pt() > maxPt) && (l1tJetIterator->eta() < 1.44) && (l1tJetIterator->eta() > -1.44))
      {
        save = true;        
        maxPt = l1tJetIterator -> pt();
        this -> _l1tObjectParticle.hwQual = l1tJetIterator -> hwQual();      
        this -> _l1tObjectParticle.id = (l1tJetIterator - l1tJetCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tJetIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tJetIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tJetIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingJetTree -> Fill();
    
    for (auto l1tJetIterator = l1tJetCollectionHandle -> begin(0); l1tJetIterator != l1tJetCollectionHandle -> end(0); l1tJetIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tJetIterator - l1tJetCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tJetIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tJetIterator -> eta();
      this -> _l1tObjectParticle.phi = l1tJetIterator -> phi();
      this -> _l1tObjectParticle.hwQual = l1tJetIterator -> hwQual();            
      this -> _l1tJetTree -> Fill();
    }
    
  }
  
  if (this -> _l1tEGammaCollectionTag)
  {
    edm::Handle < BXVector< l1t::EGamma > > l1tEGammaCollectionHandle;
    iEvent.getByToken(*(this -> _l1tEGammaCollectionTag), l1tEGammaCollectionHandle);
    auto l1tEGammaGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tEGammaCollectionHandle);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tEGammaGenJetTree), l1tEGammaGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tEGammaIterator = l1tEGammaCollectionHandle -> begin(0); l1tEGammaIterator != l1tEGammaCollectionHandle -> end(0); l1tEGammaIterator++ )
    {
      if (l1tEGammaIterator -> pt() > maxPt)
      {
        save = true;
        maxPt = l1tEGammaIterator -> pt();
        this -> _l1tObjectParticle.hwQual = l1tEGammaIterator -> hwQual();        
        this -> _l1tObjectParticle.id = (l1tEGammaIterator - l1tEGammaCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tEGammaIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tEGammaIterator -> eta();
        this -> _l1tObjectParticle.phi = l1tEGammaIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingEGammaTree -> Fill();
    
    for (auto l1tEGammaIterator = l1tEGammaCollectionHandle -> begin(0); l1tEGammaIterator != l1tEGammaCollectionHandle -> end(0); l1tEGammaIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tEGammaIterator - l1tEGammaCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tEGammaIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tEGammaIterator -> eta();
      this -> _l1tObjectParticle.phi = l1tEGammaIterator -> phi();
      this -> _l1tObjectParticle.hwQual = l1tEGammaIterator -> hwQual();                  
      this -> _l1tEGammaTree -> Fill();
    }
    
  }
  
  if (this -> _l1tTauCollectionTag)
  {
    edm::Handle < BXVector< l1t::Tau > > l1tTauCollectionHandle;
    iEvent.getByToken(*(this -> _l1tTauCollectionTag), l1tTauCollectionHandle);
    auto l1tTauGenJetPairs = this -> _matchGenJetWithL1Object<>(genJetCollectionHandle, l1tTauCollectionHandle);
    this -> _fillTreeWithMatchedPairs(*(this -> _l1tTauGenJetTree), l1tTauGenJetPairs);
    float maxPt = 0;
    bool save = false;
    for (auto l1tTauIterator = l1tTauCollectionHandle -> begin(0); l1tTauIterator != l1tTauCollectionHandle -> end(0); l1tTauIterator++ )
    {
      if (l1tTauIterator -> pt() > maxPt)
      {
        save = true;
        maxPt = l1tTauIterator -> pt();
        this -> _l1tObjectParticle.id = (l1tTauIterator - l1tTauCollectionHandle->begin(0));
        this -> _l1tObjectParticle.pt = l1tTauIterator -> pt();
        this -> _l1tObjectParticle.eta = l1tTauIterator -> eta();
        this -> _l1tObjectParticle.hwQual = l1tTauIterator -> hwQual();                        
        this -> _l1tObjectParticle.phi = l1tTauIterator -> phi();
      }
    }
    if (save) this -> _l1tLeadingTauTree -> Fill();
    
    for (auto l1tTauIterator = l1tTauCollectionHandle -> begin(0); l1tTauIterator != l1tTauCollectionHandle -> end(0); l1tTauIterator++ )
    {
      this -> _l1tObjectParticle.id = (l1tTauIterator - l1tTauCollectionHandle->begin(0));
      this -> _l1tObjectParticle.pt = l1tTauIterator -> pt();
      this -> _l1tObjectParticle.eta = l1tTauIterator -> eta();
      this -> _l1tObjectParticle.hwQual = l1tTauIterator -> hwQual();                        
      this -> _l1tObjectParticle.phi = l1tTauIterator -> phi();
      this -> _l1tTauTree -> Fill();
    }
    
  }

  // Saving every genJet to get the misidentification probability

  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    this -> _genJetParticle.id = (genJetIterator - genJetCollectionHandle->begin());
    this -> _genJetParticle.pt = genJetIterator -> pt();
    //std::cout << "Jet pt " << genJetIterator -> pt() << std::endl;
    this -> _genJetParticle.eta = genJetIterator -> eta();
    this -> _genJetParticle.phi = genJetIterator -> phi();
    this -> _genJetTree -> Fill();
  }


}

void
MatchGenJetWithL1Objects::_fillTreeWithMatchedPairs
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

// Considers a genJet and looks for every l1tObject within that radius
//If more than 1 object is present the highest momentum one is considered and saved

template <class T> // <3
const std::vector <std::tuple < TriggerObject, Particle, float > > 
MatchGenJetWithL1Objects::_matchGenJetWithL1Object
(
  const edm::Handle < std::vector< reco::GenJet > > & genJetCollectionHandle,
  const edm::Handle < BXVector < T > > & l1tObjectCollectionHandle 
)
{

  std::vector< std::tuple < TriggerObject, Particle , float > > l1tObjectGenJetPairs;
  // for each object in the l1t collection we look for the closest jet in a wide range
  for (auto genJetIterator = genJetCollectionHandle -> begin(); genJetIterator != genJetCollectionHandle -> end(); genJetIterator++ )
  {
    
    bool foundMatch = false;
    float dr2Min = 0.25; // i.e. dr = 0.5, half of the jet size, 0.5
    std::tuple<TriggerObject, Particle, float > l1tObjectGenJetPair;
    
    Particle & matchedGenJet = std::get<1>(l1tObjectGenJetPair);
    matchedGenJet.id = (genJetIterator - genJetCollectionHandle -> begin());
    matchedGenJet.pt = genJetIterator -> pt();
    matchedGenJet.phi = genJetIterator -> phi();
    matchedGenJet.eta = genJetIterator -> eta();

    for (
      typename BXVector<T>::const_iterator bx0Iterator = l1tObjectCollectionHandle -> begin(0);
      bx0Iterator != l1tObjectCollectionHandle -> end(0);
      bx0Iterator++
    )
    {
      float dr2 = reco::deltaR2(*bx0Iterator, *genJetIterator);
      
      if ((dr2 < dr2Min))
      {
        TriggerObject & matchedL1TObject = std::get<0>(l1tObjectGenJetPair);
        matchedL1TObject.id = (bx0Iterator - l1tObjectCollectionHandle -> begin(0));
        matchedL1TObject.pt = bx0Iterator -> pt();
        matchedL1TObject.phi = bx0Iterator -> phi();
        matchedL1TObject.eta = bx0Iterator -> eta();
        matchedL1TObject.hwQual = bx0Iterator -> hwQual();
        
        std::get<2>(l1tObjectGenJetPair) = dr2;
        dr2Min = dr2;
        foundMatch = true;
      }
    }
    
    //if we have found a compatible gen jet we push the object-jet pair in a vector which will be our result
    if (foundMatch)
      l1tObjectGenJetPairs.push_back(l1tObjectGenJetPair);
  }

  return l1tObjectGenJetPairs;
}

// ------------ method called once each job just before starting event loop  ------------
void MatchGenJetWithL1Objects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchGenJetWithL1Objects::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchGenJetWithL1Objects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MatchGenJetWithL1Objects);
