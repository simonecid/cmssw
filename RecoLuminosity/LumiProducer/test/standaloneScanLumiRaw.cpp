#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "LumiRawDataStructures.h"
#include <iostream>
#include <iomanip>
/** This programm scans a given lumi raw data file and print out the content
 **/
int main(int argc, char** argv){
  const char* filename="test.root";
  //default file to read. file name is taken from command argument
  if(argc>1){
    filename=argv[1];
  }
  TFile *myfile=new TFile(filename,"READ");

  HCAL_HLX::RUN_SUMMARY *myRunSummary = new HCAL_HLX::RUN_SUMMARY;
  TTree *runsummaryTree = (TTree *) myfile->Get("RunSummary");
  if(!runsummaryTree) std::cout<<"no run summary data"<<std::endl;
  runsummaryTree->SetBranchAddress("RunSummary.",&myRunSummary);
  size_t runsummaryentries=runsummaryTree->GetEntries();
  std::cout<<"n run summary entries "<<runsummaryentries<<std::endl;
  for(size_t i=0;i<runsummaryentries;++i){
    runsummaryTree->GetEntry(i);
    std::cout<<"Summary for run : "<<myRunSummary->runNumber<<std::endl;
    std::cout<<std::setw(20)<<"timestamp : "<<myRunSummary->timestamp<<" : timestamp micros : "<<myRunSummary->timestamp_micros<<" : start orbit : "<<myRunSummary->startOrbitNumber<<" : end orbit : "<<myRunSummary->endOrbitnumber<<" : fill number : "<<myRunSummary->fillNumber<<" : number CMS LS : "<<myRunSummary->numberCMSLumiSections<<" : number DAQ LS : "<<myRunSummary->numberLumiDAQLumiSections<<std::endl;
  }

  HCAL_HLX::LEVEL1_TRIGGER *myTRG = new HCAL_HLX::LEVEL1_TRIGGER;
  TTree *trgTree = (TTree *) myfile->Get("L1Trigger");
  if(!trgTree) std::cout<<"no trg data"<<std::endl;
  trgTree->SetBranchAddress("L1Trigger.",&myTRG);
  size_t trgentries=trgTree->GetEntries();
  for(size_t i=0;i<trgentries;++i){
    trgTree->GetEntry(i);
    //std::cout<<"trg runnumber "<<myTRG->runNumber<<std::endl;
    std::cout<<"TRG for run : "<< myTRG->runNumber<<" : LS : "<<myTRG->sectionNumber<<" : deadtime : "<< myTRG->deadtimecount<<std::endl;
    for( unsigned int j=0; j<128; ++j){
      std::cout<<std::setw(20)<<"GT Algo  "<<j;
      std::cout<<" : path : "<< myTRG->GTAlgo[j].pathName<<" : counts : "<< myTRG->GTAlgo[j].counts<<" : prescale : "<< myTRG->GTAlgo[j].prescale<<std::endl;
    }
    for( unsigned int k=0; k<64; ++k){
      std::cout<<std::setw(20)<<"GT Tech : "<<k;
      std::cout<<" : path : "<< myTRG->GTTech[k].pathName<<" : counts : "<< myTRG->GTTech[k].counts<<" : prescale : "<< myTRG->GTTech[k].prescale<<std::endl;
    }
  }

  HCAL_HLX::HLTRIGGER *myHLT = new HCAL_HLX::HLTRIGGER;
  TTree *hltTree = (TTree *) myfile->Get("HLTrigger");
  if(!hltTree) std::cout<<"no hlt data"<<std::endl;
  hltTree->SetBranchAddress("HLTrigger.",&myHLT);
  size_t hltentries=hltTree->GetEntries();
  for(size_t i=0;i<hltentries;++i){
    hltTree->GetEntry(i);
    std::cout<<"HLT for run : "<< myHLT->runNumber<<":  LS  : "<<myHLT->sectionNumber<<" : total hlt paths : "<<myHLT->numPaths<<std::endl;
    for( unsigned int j=0; j<myHLT->numPaths; ++j){
      std::cout<<std::setw(20)<<"HLTConfigId : "<<myHLT->HLTPaths[j].HLTConfigId<<"path : "<<myHLT->HLTPaths[j].PathName<<" : L1Pass : "<<myHLT->HLTPaths[j].L1Pass<<" : PSPass : "<<myHLT->HLTPaths[j].PSPass<<" : PAccept : "<<myHLT->HLTPaths[j].PAccept<<" : PExcept : "<<myHLT->HLTPaths[j].PExcept<<" : PReject : "<<myHLT->HLTPaths[j].PReject<<" : PSIndex : "<<myHLT->HLTPaths[j].PSIndex<<" : Prescale : "<<myHLT->HLTPaths[j].Prescale<<std::endl;
    }
  }
  HCAL_HLX::LUMI_SECTION *myLumiSection=new HCAL_HLX::LUMI_SECTION;
  HCAL_HLX::LUMI_SECTION_HEADER *myLumiHeader = &(myLumiSection->hdr);
  HCAL_HLX::LUMI_SUMMARY *myLumiSummary = &(myLumiSection->lumiSummary);
  HCAL_HLX::LUMI_DETAIL *myLumiDetail = &(myLumiSection->lumiDetail);
  
  TTree *hlxTree = (TTree *) myfile->Get("HLXData");
  if(!hlxTree) std::cout<<"no hlx data"<<std::endl;
  hlxTree->SetBranchAddress("Header.",&myLumiHeader);
  hlxTree->SetBranchAddress("Summary.",&myLumiSummary);
  hlxTree->SetBranchAddress("Detail.",&myLumiDetail);
  size_t hlxentries=hlxTree->GetEntries();
  std::cout<<"hlxentries "<<hlxentries<<std::endl;
  for(size_t i=0;i<hlxentries;++i){
    hlxTree->GetEntry(i);
    std::cout<<"Lumi summary for run : "<<myLumiHeader->runNumber<<" : LS : "<<myLumiHeader->sectionNumber<<std::endl;
    std::cout<<std::setw(20)<<"deadtime norm : "<<myLumiSummary->DeadTimeNormalization<<" : LHC norm : "<<myLumiSummary->LHCNormalization<<" : instantlumi : "<<myLumiSummary->InstantLumi<<" : instantlumiErr : "<<myLumiSummary->InstantLumiErr<<" : instantlumiQlty : "<<myLumiSummary->InstantLumiQlty<<std::endl;
    std::cout<<std::setw(20)<<"lumi details : "<<std::endl;
    for(size_t j=0;j<HCAL_HLX_MAX_BUNCHES;++j){
      std::cout<<std::setw(20)<<"    LHCLumi : "<<myLumiDetail->LHCLumi[j]<<" : ETLumi : "<<myLumiDetail->ETLumi[j]<<" : ETLumiErr : "<<myLumiDetail->ETLumiErr[j]<<" : ETLumiQlty : "<<myLumiDetail->ETLumiQlty[j]<<" : ETBXNormalization : "<<myLumiDetail->ETBXNormalization[j]<<std::endl;
    }
  }
}

