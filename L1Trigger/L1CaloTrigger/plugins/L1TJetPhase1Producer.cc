// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: L1TJetPhase1Producer
//
/**\class L1TJetPhase1Producer L1TJetPhase1Producer.cc L1Trigger/L1CaloTrigger/plugin/L1TJetPhase1Producer.cc

Description: Produces jets with sliding window algorithm using pfcluster and pfcandidates

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"

#include "TH2F.h"

class L1TJetPhase1Producer : public edm::EDProducer {
   public:
      explicit L1TJetPhase1Producer(const edm::ParameterSet&);
      ~L1TJetPhase1Producer();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);
      void _buildJets(TH2F* caloGrid);
      edm::EDGetTokenT<std::vector<l1t::PFCandidate>> *_pfCandidateCollectionTag;
      edm::EDGetTokenT<std::vector<l1t::PFCluster>> *_pfClusterCollectionTag;
      TH2F* _caloGrid;

};

L1TJetPhase1Producer::L1TJetPhase1Producer(const edm::ParameterSet& iConfig)
{
  // Retrieving the pfCandidates and pfClusters if input tag has been provided
  try
  {
    this -> _pfCandidateCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCandidate> >(consumes< std::vector<l1t::PFCandidate> > (iConfig.getParameter< edm::InputTag >("pfCandidateCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfCandidateCollectionTag not found" << std::endl;
    this -> _pfCandidateCollectionTag = NULL;
  }
  
  try
  {
    this -> _pfClusterCollectionTag = new edm::EDGetTokenT< std::vector<l1t::PFCluster> >(consumes< std::vector<l1t::PFCluster> > (iConfig.getParameter< edm::InputTag >("pfClusterCollectionTag")));
  } catch (std::exception const & ex) 
  {
    std::cerr << ">>> pfClusterCollectionTag not found" << std::endl;
    this -> _pfClusterCollectionTag = NULL;
  }

  std::vector<double> etaBinning = iConfig.getParameter<std::vector<double> >("etaBinning");
  size_t nBinsEta = etaBinning.size() - 1;
  unsigned int nBinsPhi = iConfig.getParameter<unsigned int>("nBinsPhi");
  double phiLow = iConfig.getParameter<double>("phiLow");
  double phiUp = iConfig.getParameter<double>("phiUp");
  
  this -> _caloGrid = new TH2F("CaloEtaPhiGrid", "Calorimeter grid", nBinsEta, etaBinning.data(), nBinsPhi, phiLow, phiUp);
  this -> _caloGrid -> GetXaxis() -> SetTitle("#eta");
  this -> _caloGrid -> GetYaxis() -> SetTitle("#phi");

}

L1TJetPhase1Producer::~L1TJetPhase1Producer()
{
  if (this -> _pfCandidateCollectionTag) delete this -> _pfCandidateCollectionTag;
  if (this -> _pfClusterCollectionTag) delete this -> _pfClusterCollectionTag;
  delete this -> _caloGrid;
}

void L1TJetPhase1Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Retrieving the pfCandidates and pfClusters if input tag has been provided

  if (this -> _pfCandidateCollectionTag) {
    edm::Handle < std::vector< l1t::PFCandidate > > pfCandidateCollectionHandle;
    iEvent.getByToken(*(this -> _pfCandidateCollectionTag), pfCandidateCollectionHandle);
    // dumping the data
    std::cout << ">>>>>> DUMPING PFCANDIDATES <<<<<<" << std::endl;
    for (auto pfCandidateIterator = pfCandidateCollectionHandle -> begin(); pfCandidateIterator != pfCandidateCollectionHandle -> end(); pfCandidateIterator++) 
    {
      std::cout << pfCandidateIterator -> pt() << "\t" << pfCandidateIterator -> eta() << "\t" << pfCandidateIterator -> phi() << "\t" << std::endl;
    }
  }

  if (this -> _pfClusterCollectionTag) {
    edm::Handle < std::vector< l1t::PFCluster > > pfClusterCollectionHandle;
    iEvent.getByToken(*(this -> _pfClusterCollectionTag), pfClusterCollectionHandle);
    // dumping the data
    std::cout << ">>>>>> DUMPING PFCLUSTERS <<<<<<" << std::endl;
    for (auto pfClusterIterator = pfClusterCollectionHandle -> begin(); pfClusterIterator != pfClusterCollectionHandle -> end(); pfClusterIterator++) 
    {
      std::cout << pfClusterIterator -> pt() << "\t" << pfClusterIterator -> eta() << "\t" << pfClusterIterator -> phi() << "\t" << std::endl;
    }
  }

  return;

}


void L1TJetPhase1Producer::_buildJets(TH2F* caloGrid) 
{
  int nBinsX = caloGrid -> GetNbinsX();
  int nBinsY = caloGrid -> GetNbinsY();

  // for each point of the grid check if it is a local maximum
  // to do so I take a point, and look if is greater than the points around it (in the 3x3 neighborhood)
  // to prevent mutual exclusion, I check greater or equal for points above and right to the one I am considering (including the top-left point)
  // to prevent mutual exclusion, I check greater for points below and left to the one I am considering (including the bottom-right point)

  for (int y = 1; y <= nBinsY; y++)
  {
    for (int x = 1; x <= nBinsX; x++)
    {
      float centralPt = caloGrid -> GetBinContent(x, y);
      if (centralPt <= 0.) continue;
      bool isLocalMaximum = true;

      // Y gets recomputed in order to take account of periodicity in phi
      int oldY = y;
      if (y == nBinsY) {
        y = 0;
      }
      //top right angle (greater/equal)
      if (x > 1) isLocalMaximum = ((isLocalMaximum) && (centralPt >= caloGrid -> GetBinContent(x - 1, y + 1)));
      isLocalMaximum = ((isLocalMaximum) && (centralPt >= caloGrid -> GetBinContent(x, y + 1)));
      if (x < nBinsX) isLocalMaximum = ((isLocalMaximum) && (centralPt >= caloGrid -> GetBinContent(x + 1, y + 1)));
      y = oldY;

      if (x < nBinsX) isLocalMaximum = ((isLocalMaximum) && (centralPt >= caloGrid -> GetBinContent(x + 1, y)));


      //bottom left angle (greater)
      if (x > 1) isLocalMaximum = ((isLocalMaximum) && (centralPt > caloGrid -> GetBinContent(x - 1, y)));

      // Y gets recomputed in order to take account of periodicity in phi      
      oldY = y;
      if (y == 1) {
        y = nBinsY + 1;
      }
      if (x > 1) isLocalMaximum = ((isLocalMaximum) && (centralPt > caloGrid -> GetBinContent(x - 1 , y - 1)));
      isLocalMaximum = ((isLocalMaximum) && (centralPt > caloGrid -> GetBinContent(x, y - 1)));
      if (x < nBinsX) isLocalMaximum = ((isLocalMaximum) && (centralPt > caloGrid -> GetBinContent(x + 1, y - 1)));
      y = oldY;

      if (isLocalMaximum) {
        std::cout << "I found a local maximum at " << caloGrid -> GetXaxis() -> GetBinCenter(x) << "\t" << caloGrid -> GetXaxis() -> GetBinCenter(y) << std::endl;
      }
    }
  }

  return;
}

DEFINE_FWK_MODULE(L1TJetPhase1Producer);