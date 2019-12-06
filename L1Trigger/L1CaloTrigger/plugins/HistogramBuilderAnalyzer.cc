// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: HistogramBuilderAnalyzer
//
/**\class HistogramBuilderAnalyzer HistogramBuilderAnalyzer.cc L1Trigger/L1CaloTrigger/plugin/HistogramBuilderAnalyzer.cc

Description: analyzes histograms pictures using ROOT and saves them to a folder

*** INPUT PARAMETERS ***
  * etaBinning: vdouble with eta binning (allows non-homogeneous binning in eta)
  * nBinsPhi: uint32, number of bins in phi
  * phiLow: double, min phi (typically -pi)
  * phiUp: double, max phi (typically +pi)
  * jetIEtaSize: uint32, jet cluster size in ieta
  * jetIPhiSize: uint32, jet cluster size in iphi
  * seedPtThreshold: double, threshold of the seed tower
  * puSubtraction: bool, runs chunky doughnut pile-up subtraction, 9x9 jet only
  * outputFolderName: string, tag for the output collection
  * vetoZeroPt: bool, controls whether jets with 0 pt should be save. 
    It matters if PU is ON, as you can get negative or zero pt jets after it.
  * inputCollectionTag: inputtag, collection of reco::candidates used as input to the algo

*/
//
// Original Simone Bologna
// Created: Mon Jul 02 2018
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
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
#include "TCanvas.h"

#include "TH2F.h"

#include <algorithm>

class HistogramBuilderAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit HistogramBuilderAnalyzer(const edm::ParameterSet&);
      ~HistogramBuilderAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      
      // <3 handy method to fill the calogrid with whatever type
      template <class Container>
      void _fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives);

      edm::EDGetTokenT<edm::View<reco::Candidate>> *_inputCollectionTag;
      // histogram containing our clustered inputs
      TH2F* _caloGrid;

      std::vector<double> _etaBinning;
      size_t _nBinsEta;
      unsigned int _nBinsPhi;
      double _phiLow;
      double _phiUp;

      std::string _outputFolderName;

};


HistogramBuilderAnalyzer::HistogramBuilderAnalyzer(const edm::ParameterSet& iConfig)
{

  // getting configuration settings
  this -> _etaBinning = iConfig.getParameter<std::vector<double> >("etaBinning");
  this -> _nBinsEta = this -> _etaBinning.size() - 1;
  this -> _nBinsPhi = iConfig.getParameter<unsigned int>("nBinsPhi");
  this -> _phiLow = iConfig.getParameter<double>("phiLow");
  this -> _phiUp = iConfig.getParameter<double>("phiUp");
  this -> _outputFolderName = iConfig.getParameter<std::string>("outputFolderName");

  this -> _inputCollectionTag = new edm::EDGetTokenT< edm::View<reco::Candidate> >(consumes< edm::View<reco::Candidate> > (iConfig.getParameter< edm::InputTag >("inputCollectionTag")));
  this -> _caloGrid = new TH2F("caloGrid", "Calorimeter grid", this -> _nBinsEta, this-> _etaBinning.data(), this -> _nBinsPhi, this -> _phiLow, this -> _phiUp);
  this -> _caloGrid -> GetXaxis() -> SetTitle("#eta");
  this -> _caloGrid -> GetYaxis() -> SetTitle("#phi");
  
}

HistogramBuilderAnalyzer::~HistogramBuilderAnalyzer()
{
  delete this -> _inputCollectionTag;
  if (this -> _caloGrid) delete this -> _caloGrid;
}

void HistogramBuilderAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle < edm::View< reco::Candidate > > inputCollectionHandle;
  iEvent.getByToken(*(this -> _inputCollectionTag), inputCollectionHandle);
  // dumping the data
  this -> _caloGrid -> Reset();
  this -> _fillCaloGrid<>(*(this -> _caloGrid), *inputCollectionHandle);
  edm::RunNumber_t lRunNumber = iEvent.run();
  edm::LuminosityBlockNumber_t lLumiBlock = iEvent.luminosityBlock();
  edm::EventNumber_t lEventNumber = iEvent.id().event();
  TCanvas lCanvas("canvas", "canvas", 1024, 1024);
  lCanvas.SetCanvasSize(1024, 1024);
  this -> _caloGrid -> SetStats(false);
  this -> _caloGrid -> Draw("COL");
  std::string lFileName;
  std::ostringstream lISS;
  lISS << this -> _outputFolderName << "/event_" << lRunNumber << "_" << lLumiBlock << "_" << lEventNumber << ".png";
  lCanvas.SaveAs(lISS.str().c_str(), "png");
  return;

}

template <class Container>
void HistogramBuilderAnalyzer::_fillCaloGrid(TH2F & caloGrid, const Container & triggerPrimitives)
{
  //Filling the calo grid with the primitives
  for (auto primitiveIterator = triggerPrimitives.begin(); primitiveIterator != triggerPrimitives.end(); primitiveIterator++)
    {
      //if(primitiveIterator->eta() >= 0 && primitiveIterator->eta() < 1.5 && primitiveIterator->phi() >= 0 && primitiveIterator->phi() < 0.7)
      caloGrid.Fill((float) primitiveIterator -> eta(), (float) primitiveIterator -> phi(), (float) primitiveIterator -> pt());
    }
  return;
}

void
HistogramBuilderAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(HistogramBuilderAnalyzer);
