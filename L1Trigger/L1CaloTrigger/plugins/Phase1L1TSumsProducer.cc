// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase1L1TSumsProducer
//
/**\class Phase1L1TSumsProducer Phase1L1TSumsProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase1L1TSumsProducer.cc

Description: Produces jets with a phase-1 like sliding window algorithm using a collection of reco::Candidates in input.
  If flag is enabled, computes MET and MHT.

*** INPUT PARAMETERS ***
  * etaBinning: vdouble with eta binning (allows non-homogeneous binning in eta)
  * nBinsPhi: uint32, number of bins in phi
  * phiLow: double, min phi (typically -pi)
  * phiUp: double, max phi (typically +pi)
  * jetIEtaSize: uint32, jet cluster size in ieta
  * jetIPhiSize: uint32, jet cluster size in iphi
  * seedPtThreshold: double, threshold of the seed tower
  * puSubtraction: bool, runs chunky doughnut pile-up subtraction, 9x9 jet only
  * outputCollectionName: string, tag for the output collection
  * vetoZeroPt: bool, controls whether jets with 0 pt should be save. 
    It matters if PU is ON, as you can get negative or zero pt jets after it.
  * outputSumsCollectionName: string, tag for the output collection containing MET and HT
  * inputCollectionTag: inputtag, collection of reco::candidates used as input to the algo
  * htPtThreshold: double, threshold for computing ht

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "TH2F.h"

#include <algorithm>

class Phase1L1TSumsProducer : public edm::one::EDProducer<edm::one::SharedResources> {
   public:
      explicit Phase1L1TSumsProducer(const edm::ParameterSet&);
      ~Phase1L1TSumsProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
      l1t::EtSum _computeHT(const std::vector<reco::CaloJet>& l1jetVector);
      template <class ParticleCollection>
      l1t::EtSum _computeMET(const ParticleCollection & caloGrid);

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_particleCollectionTag;
      edm::EDGetTokenT<std::vector<reco::CaloJet> > *_jetCollectionTag;

      std::vector<double> _sinPhi;
      std::vector<double> _cosPhi;
      unsigned int _nBinsPhi;
      double _phiLow;
      double _phiUp;
      double _phiStep;
      double _htPtThreshold;

      std::string _outputCollectionName;

};


Phase1L1TSumsProducer::Phase1L1TSumsProducer(const edm::ParameterSet& iConfig):
  // getting configuration settings
  _sinPhi(iConfig.getParameter<std::vector<double> >("sinPhi")),
  _cosPhi(iConfig.getParameter<std::vector<double> >("cosPhi")),
  _nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
  _phiLow(iConfig.getParameter<double>("phiLow")),
  _phiUp(iConfig.getParameter<double>("phiUp")),
  _htPtThreshold(iConfig.getParameter<double>("htPtThreshold")),
  _outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName"))
{
  this -> _particleCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("particleCollectionTag")));  
  this -> _jetCollectionTag = new edm::EDGetTokenT< std::vector<reco::CaloJet> >(consumes< std::vector<reco::CaloJet> > (iConfig.getParameter< edm::InputTag >("jetCollectionTag")));  
  this -> _phiStep = ( this -> _phiUp - this -> _phiLow ) / this -> _nBinsPhi;
  produces< BXVector<l1t::EtSum> >( this -> _outputCollectionName ).setBranchAlias(this -> _outputCollectionName);
}

Phase1L1TSumsProducer::~Phase1L1TSumsProducer()
{
  delete this -> _particleCollectionTag;
  delete this -> _jetCollectionTag;
}

void Phase1L1TSumsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle < edm::View< reco::Candidate > > particleCollectionHandle;
  edm::Handle < std::vector<reco::CaloJet> > jetCollectionHandle;
  iEvent.getByToken(*(this -> _particleCollectionTag), particleCollectionHandle);
  iEvent.getByToken(*(this -> _jetCollectionTag), jetCollectionHandle);
  
  // computing sums
  l1t::EtSum lHT = this -> _computeHT(*jetCollectionHandle);
  l1t::EtSum lMET = this -> _computeMET<>(*particleCollectionHandle);

  //packing in BXVector for event saving
  std::unique_ptr< BXVector<l1t::EtSum> > lSumVectorPtr(new BXVector<l1t::EtSum>(2));
  lSumVectorPtr -> push_back(0, lHT);
  lSumVectorPtr -> push_back(0, lMET);
  // std::cout << "HT-MET sums prod: " << lHT.pt() << "\t" << lMET.pt() << std::endl;
  //saving sums
  iEvent.put(std::move(lSumVectorPtr), this -> _outputCollectionName);

  return;

}

// computes ht, adds jet pt to ht only if the pt of the jet is above the ht calculation threshold
l1t::EtSum Phase1L1TSumsProducer::_computeHT(const std::vector<reco::CaloJet>& l1jetVector) 
{
  double lHT = 0;
  for (const auto & jet: l1jetVector)
  {
    lHT += (jet.pt() >= this -> _htPtThreshold) ? jet.pt() : 0;
  }

  reco::Candidate::PolarLorentzVector lHTVector;
  lHTVector.SetPt(lHT);
  lHTVector.SetEta(0);
  lHTVector.SetPhi(0);
  l1t::EtSum lHTSum(lHTVector, l1t::EtSum::EtSumType::kTotalHt);

  return lHTSum;
}

template <class ParticleCollection>
l1t::EtSum Phase1L1TSumsProducer::_computeMET(const ParticleCollection & particleCollection) 
{

  double lTotalPx = 0;
  double lTotalPy = 0;

  for (const auto & particle : particleCollection)
  {
    double lParticlePhi = particle.phi();
    double lParticlePt = particle.pt();
    if ((lParticlePhi < this -> _phiLow) || (lParticlePhi > this -> _phiUp)) continue;
    unsigned int iPhi = ( lParticlePhi - this -> _phiLow ) / this -> _phiStep;
    double lSinPhi = this -> _sinPhi[iPhi];
    double lCosPhi = this -> _cosPhi[iPhi];
    lTotalPx += (lParticlePt * lCosPhi);
    lTotalPy += (lParticlePt * lSinPhi);
  }

  double lMET = sqrt(lTotalPx * lTotalPx + lTotalPy * lTotalPy);
  //packing in EtSum object
  reco::Candidate::PolarLorentzVector lMETVector;
  lMETVector.SetPt(lMET);
  lMETVector.SetEta(0);
  lMETVector.SetPhi(0);
  l1t::EtSum lMETSum(lMETVector, l1t::EtSum::EtSumType::kMissingEt);

  return lMETSum;

}

void Phase1L1TSumsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Phase1L1TSumsProducer);
