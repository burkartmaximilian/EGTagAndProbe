#ifndef NTUPLIZER_H
#define NTUPLIZER_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <TNtuple.h>


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"


class Ntuplizer : public edm::EDAnalyzer {
    public:
        /// Constructor
        explicit Ntuplizer(const edm::ParameterSet&);
        /// Destructor
        virtual ~Ntuplizer();  
        
    private:
        //----edm control---
        virtual void beginJob() ;
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        void Initialize(); 
        
        TTree *_tree;
        std::string _treeName;
        // -------------------------------------
        // variables to be filled in output tree
        ULong64_t       _indexevents;
        Int_t           _runNumber;
        Int_t           _lumi;

        edm::EDGetTokenT<pat::MuonRefVector>  _muonsTag;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _triggerObjects;
        edm::EDGetTokenT<edm::TriggerResults> _triggerBits;

        edm::InputTag _processName;

        std::vector<Float_t> _tauPtVector;
        std::vector<uint> _tauTriggeredVector;

        unsigned int _doubleMediumIsoPFTau32Index;
        unsigned int _doubleMediumIsoPFTau35Index;
        unsigned int _doubleMediumIsoPFTau40Index;

        HLTConfigProvider _hltConfig;

        
};


// ----Constructor and Destructor -----
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
_muonsTag(consumes<pat::MuonRefVector>(iConfig.getParameter<edm::InputTag>("muons"))),
_triggerObjects(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerSet"))),
_triggerBits(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultsLabel")))
{
    this -> _treeName = iConfig.getParameter<std::string>("treeName");
    this -> _processName = iConfig.getParameter<edm::InputTag>("triggerResultsLabel");
    this -> Initialize();
    return;
}

Ntuplizer::~Ntuplizer()
{}

void Ntuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    Bool_t changedConfig = false;
    
    if(!this -> _hltConfig.init(iRun, iSetup, this -> _processName.process(), changedConfig)){
        edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!"; 
        return;
    }

    this -> _doubleMediumIsoPFTau32Index = -1;
    this -> _doubleMediumIsoPFTau35Index = -1;
    this -> _doubleMediumIsoPFTau40Index = -1;

    const edm::TriggerNames::Strings& triggerNames = this -> _hltConfig.triggerNames();
    for(unsigned int j=0; j < triggerNames.size(); j++)
    {
        // if (triggerNames[j].find("HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg_v") != std::string::npos) this -> _doubleMediumIsoPFTau32Index = j;
        if (triggerNames[j].find("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v") != std::string::npos) this -> _doubleMediumIsoPFTau32Index = j;
        
        if (triggerNames[j].find("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v") != std::string::npos) this -> _doubleMediumIsoPFTau35Index = j;
        if (triggerNames[j].find("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v") != std::string::npos) this -> _doubleMediumIsoPFTau40Index = j;
        std::cout << j << " -- " << triggerNames[j] << std::endl;

    } 
    
}

void Ntuplizer::Initialize() {
    this -> _indexevents = 0;
    this -> _runNumber = 0;
    this -> _lumi = 0;

    this -> _tauPtVector.clear();    
}


void Ntuplizer::beginJob()
{
    edm::Service<TFileService> fs;
    this -> _tree = fs -> make<TTree>(this -> _treeName.c_str(), this -> _treeName.c_str());

    //Branches
    this -> _tree -> Branch("EventNumber",&_indexevents,"EventNumber/l");
    this -> _tree -> Branch("RunNumber",&_runNumber,"RunNumber/I");
    this -> _tree -> Branch("lumi",&_lumi,"lumi/I");
    this -> _tree -> Branch("tauPt", &_tauPtVector);
    this -> _tree -> Branch("isTriggered", &_tauTriggeredVector);
    
    return;
}


void Ntuplizer::endJob()
{
    return;
}


void Ntuplizer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    return;
}

void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& eSetup)
{

    this -> Initialize();

    _indexevents = iEvent.id().event();
    _runNumber = iEvent.id().run();
    _lumi = iEvent.luminosityBlock();

    std::auto_ptr<pat::MuonRefVector> resultMuon(new pat::MuonRefVector);

    // search for the tag in the event
    edm::Handle<pat::MuonRefVector> muonHandle;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<edm::TriggerResults> triggerBits;

    iEvent.getByToken(this -> _muonsTag, muonHandle);
    iEvent.getByToken(this -> _triggerObjects, triggerObjects);
    iEvent.getByToken(this -> _triggerBits, triggerBits);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    
    int idx = 0;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects)
    {
        obj.unpackPathNames(names);
        const edm::TriggerNames::Strings& triggerNames = names.triggerNames();

        uint isTriggered = 0;

        // if ((obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau32Index], true, false)) || 
        //     (obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau35Index], true, false)) ||
        //     (obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau40Index], true, false))) isTriggered = 1;
//last , L3
        bool nameTF = obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau32Index], true, false);
        bool nameFT = obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau32Index], false, true);
        bool nameFF = obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau32Index], false, false);
        bool nameTT = obj.hasPathName(triggerNames[this -> _doubleMediumIsoPFTau32Index], true, true);


        std::string found = triggerNames[this -> _doubleMediumIsoPFTau32Index];
        found += std::string("AAA");
        // // std::cout << "fOUND: " << triggerNames[this -> _doubleMediumIsoPFTau32Index] << std::endl;

        // if ( obj.pt() > 40 ) 
        //      std::cout << nameTF << " " << nameFT << " " << nameFF << " " << nameTT << std::endl;

        const std::vector<std::string>& vLabels = obj.filterLabels();
        for (std::string str : vLabels)
        {
            if (str == std::string("hltDoublePFTau35TrackPt1MediumIsolationDz02Reg"))
            {
                std::cout << idx << " trovato " << obj.hasTriggerObjectType(trigger::TriggerTau) << " " << obj.hasTriggerObjectType(trigger::TriggerL1TauJet) << " " <<
                nameTF << " " << nameFT << " " << nameFF << " " << nameTT << " " <<
                obj.pt() << " " << obj.eta() << " " << obj.phi() << " " << obj.energy() << std::endl;
                // break;
            }
        }

        this -> _tauPtVector.push_back(obj.pt());
        this -> _tauTriggeredVector.push_back(isTriggered);
        ++idx;
    }

    /*
    for (size_t imu = 0; imu < muonHandle -> size(); ++imu )
    {
        const pat::MuonRef mu = (*muonHandle)[imu] ;
        std::cout << "##### MUON PT: " << mu -> pt() << std::endl;
        this -> _muonsPtVector.push_back(mu -> pt());
    }
    */ 

    this -> _tree -> Fill();
    
}

#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(Ntuplizer);

#endif //NTUPLIZER_H
