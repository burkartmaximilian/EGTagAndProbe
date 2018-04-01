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
#include <TString.h>
#include <bitset>


#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "tParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"



//Set this variable to decide the number of triggers that you want to check simultaneously
#define NUMBER_OF_MAXIMUM_TRIGGERS 64


/*
██████  ███████  ██████ ██       █████  ██████   █████  ████████ ██  ██████  ███    ██
██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ████   ██
██   ██ █████   ██      ██      ███████ ██████  ███████    ██    ██ ██    ██ ██ ██  ██
██   ██ ██      ██      ██      ██   ██ ██   ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
██████  ███████  ██████ ███████ ██   ██ ██   ██ ██   ██    ██    ██  ██████  ██   ████
*/

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
        bool hasFilters(const pat::TriggerObjectStandAlone&  obj , const std::vector<std::string>& filtersToLookFor);
        bool matchToTruth(const edm::Ptr<reco::GsfElectron> ele, 
			  const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

        TTree *_tree;
        TTree *_triggerNamesTreeTag;
        TTree *_triggerNamesTreeProbe;

        std::string _treeName;
        // -------------------------------------
        // variables to be filled in output tree
        ULong64_t       _indexevents;
        Int_t           _runNumber;
        Int_t           _lumi;
        unsigned long _eleProbeTriggerBits;
        unsigned long _eleTagTriggerBits;
        float _eleProbePt;
        float _eleProbeEta;
        float _eleProbePhi;
        float _eleProbeSclEt;
        int _eleProbeCharge;

        float _tauProbePt;
        float _tauProbeEta;
        float _tauProbePhi;       
        int _tauProbeCharge;
        int _tauProbeDM;
        bool _tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT;
        bool _tauProbeByLooseIsolationMVArun2v1DBoldDMwLT;
        bool _tauProbeByMediumIsolationMVArun2v1DBoldDMwLT;
        bool _tauProbeByTightIsolationMVArun2v1DBoldDMwLT;
        bool _tauProbeByVTightIsolationMVArun2v1DBoldDMwLT;
        bool _tauProbeAgainstMuonLoose3;
        bool _tauProbeAgainstMuonTight3;
        bool _tauProbeAgainstElectronVLooseMVA6;
        bool _tauProbeAgainstElectronLooseMVA6;
        bool _tauProbeAgainstElectronMediumMVA6;
        bool _tauProbeAgainstElectronTightMVA6;
        bool _tauProbeAgainstElectronVTightMVA6;

        vector<float> _hltPt;
        vector<float> _hltEta;
        vector<float> _hltPhi;
        int _l1tQual;
        float _l1tPt;
        float _l1tEta;
        float _l1tPhi;
        int _l1tIso;
        int _l1tEmuQual;
        float _l1tEmuPt;
        float _l1tEmuEta;
        float _l1tEmuPhi;
        int _l1tEmuIso;
        int _l1tEmuNTT;
        int _l1tEmuTowerIEta;
        int _l1tEmuTowerIPhi;
        int _l1tEmuRawEt;
        int _l1tEmuIsoEt;
        Bool_t _isTagHLTmatched;
        Bool_t _isProbeHLTmatched;
        Bool_t _isOS;
        int _foundTag;
        float _eleTagPt;
        float _eleTagEta;
        float _eleTagPhi;
        int _eleTagCharge;
        float _Mee;
        int _Nvtx;
        
        int _hasL1[100];
        int _hasL1_iso[100];

        edm::EDGetTokenT<edm::View<reco::GsfElectron> >  _electronsTag;
        edm::EDGetTokenT<pat::TauRefVector>   _tauTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle> > _genParticlesTag;
        edm::EDGetTokenT<edm::ValueMap<bool> > _eleLooseIdMapTag;
        edm::EDGetTokenT<edm::ValueMap<bool> > _eleMediumIdMapTag;
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _triggerObjects;
        edm::EDGetTokenT<edm::TriggerResults> _triggerBits;
        edm::EDGetTokenT<l1t::EGammaBxCollection> _L1EGTag  ;
        edm::EDGetTokenT<l1t::EGammaBxCollection> _L1EmuEGTag  ;
        edm::EDGetTokenT<std::vector<reco::Vertex>> _VtxTag;

        bool _useGenMatch;
        bool _useHLTMatch;

        //!Contains the parameters
        tVParameterSet _parametersTag;
        tVParameterSet _parametersProbe;

        edm::InputTag _processName;
        //! Maximum
        std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> _eleProbeTriggerBitSet;
        std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> _eleTagTriggerBitSet;



        HLTConfigProvider _hltConfig;


};

/*
██ ███    ███ ██████  ██      ███████ ███    ███ ███████ ███    ██ ████████  █████  ████████ ██  ██████  ███    ██
██ ████  ████ ██   ██ ██      ██      ████  ████ ██      ████   ██    ██    ██   ██    ██    ██ ██    ██ ████   ██
██ ██ ████ ██ ██████  ██      █████   ██ ████ ██ █████   ██ ██  ██    ██    ███████    ██    ██ ██    ██ ██ ██  ██
██ ██  ██  ██ ██      ██      ██      ██  ██  ██ ██      ██  ██ ██    ██    ██   ██    ██    ██ ██    ██ ██  ██ ██
██ ██      ██ ██      ███████ ███████ ██      ██ ███████ ██   ████    ██    ██   ██    ██    ██  ██████  ██   ████
*/

// ----Constructor and Destructor -----
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
_electronsTag       (consumes<edm::View<reco::GsfElectron> >                     (iConfig.getParameter<edm::InputTag>("electrons"))),
_tauTag         (consumes<pat::TauRefVector>                      (iConfig.getParameter<edm::InputTag>("taus"))),
_genParticlesTag (consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genParticles"))),
_eleLooseIdMapTag  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
_eleMediumIdMapTag  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
_triggerObjects (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerSet"))),
_triggerBits    (consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerResultsLabel"))),
_L1EGTag       (consumes<l1t::EGammaBxCollection>                   (iConfig.getParameter<edm::InputTag>("L1EG"))),
_L1EmuEGTag    (consumes<l1t::EGammaBxCollection>                   (iConfig.getParameter<edm::InputTag>("L1EmuEG"))),
_VtxTag         (consumes<std::vector<reco::Vertex>>              (iConfig.getParameter<edm::InputTag>("Vertices")))
{
    _treeName = iConfig.getParameter<std::string>("treeName");
    _processName = iConfig.getParameter<edm::InputTag>("triggerResultsLabel");
    _useGenMatch = iConfig.getParameter<bool>("useGenMatch");
    _useHLTMatch = iConfig.getParameter<bool>("useHLTMatch");

    TString triggerNameTag;
    edm::Service<TFileService> fs;
    _triggerNamesTreeTag = fs -> make<TTree>("triggerNamesTag", "triggerNamesTag");
    _triggerNamesTreeTag -> Branch("triggerNamesTag",&triggerNameTag);

    //Building the trigger arrays
    const std::vector<edm::ParameterSet>& HLTListTag = iConfig.getParameter <std::vector<edm::ParameterSet> > ("triggerListTag");
    for (const edm::ParameterSet& parameterSet : HLTListTag) {
        tParameterSet pSet;
        pSet.hltPath = parameterSet.getParameter<std::string>("HLT");
        triggerNameTag = pSet.hltPath;
        pSet.hltFilters1 = parameterSet.getParameter<std::vector<std::string> >("path1");
        pSet.hltFilters2 = parameterSet.getParameter<std::vector<std::string> >("path2");
        pSet.leg1 = parameterSet.getParameter<int>("leg1");
        pSet.leg2 = parameterSet.getParameter<int>("leg2");
        _parametersTag.push_back(pSet);

        _triggerNamesTreeTag -> Fill();
    }


    TString triggerNameProbe;
    _triggerNamesTreeProbe = fs -> make<TTree>("triggerNamesProbe", "triggerNamesProbe");
    _triggerNamesTreeProbe -> Branch("triggerNamesProbe",&triggerNameProbe);

    //Building the trigger arrays
    const std::vector<edm::ParameterSet>& HLTListProbe = iConfig.getParameter <std::vector<edm::ParameterSet> > ("triggerListProbe");
    for (const edm::ParameterSet& parameterSet : HLTListProbe) {
        tParameterSet pSet;
        pSet.hltPath = parameterSet.getParameter<std::string>("HLT");
        triggerNameProbe = pSet.hltPath;
        pSet.hltFilters1 = parameterSet.getParameter<std::vector<std::string> >("path1");
        pSet.hltFilters2 = parameterSet.getParameter<std::vector<std::string> >("path2");
        pSet.leg1 = parameterSet.getParameter<int>("leg1");
        pSet.leg2 = parameterSet.getParameter<int>("leg2");
        _parametersProbe.push_back(pSet);

        _triggerNamesTreeProbe -> Fill();
    }


    Initialize();
    return;
}

Ntuplizer::~Ntuplizer()
{}

void Ntuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    Bool_t changedConfig = false;

    if(!_hltConfig.init(iRun, iSetup, _processName.process(), changedConfig)){
        edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
        return;
    }

    const edm::TriggerNames::Strings& triggerNames = _hltConfig.triggerNames();
    //std::cout << " ===== LOOKING FOR THE PATH INDEXES =====" << std::endl;
    for (tParameterSet& parameter : _parametersTag){
        const std::string& hltPath = parameter.hltPath;
        bool found = false;
        for(unsigned int j=0; j < triggerNames.size(); j++)
        {
            //std::cout << triggerNames[j] << std::endl;
            if (triggerNames[j].find(hltPath) != std::string::npos) {
                found = true;
                parameter.hltPathIndex = j;

                //std::cout << "### FOUND AT INDEX #" << j << " --> " << triggerNames[j] << std::endl;
            }
        }
        if (!found) parameter.hltPathIndex = -1;
    }


    for (tParameterSet& parameter : _parametersProbe){
        const std::string& hltPath = parameter.hltPath;
        bool found = false;
        for(unsigned int j=0; j < triggerNames.size(); j++)
        {
            //std::cout << triggerNames[j] << std::endl;
            if (triggerNames[j].find(hltPath) != std::string::npos) {
                found = true;
                parameter.hltPathIndex = j;

                //std::cout << "### FOUND AT INDEX #" << j << " --> " << triggerNames[j] << std::endl;
            }
        }
        if (!found) parameter.hltPathIndex = -1;
    }

}

void Ntuplizer::Initialize() {
    _indexevents = 0;
    _runNumber = 0;
    _lumi = 0;
    _eleProbePt = -1.;
    _eleProbeEta = -1.;
    _eleProbePhi = -1.;
    _eleProbeSclEt = -1.;
    _eleProbeCharge = 0;

    _tauProbePt = -1.;
    _tauProbeEta = -1.;
    _tauProbePhi = -1.;
    _tauProbeCharge = 0;
    _tauProbeDM = -1;
    _tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT = 0;
    _tauProbeByLooseIsolationMVArun2v1DBoldDMwLT = 0;
    _tauProbeByMediumIsolationMVArun2v1DBoldDMwLT = 0;
    _tauProbeByTightIsolationMVArun2v1DBoldDMwLT = 0;
    _tauProbeByVTightIsolationMVArun2v1DBoldDMwLT = 0;
    _tauProbeAgainstMuonLoose3 = 0;
    _tauProbeAgainstMuonTight3 = 0;
    _tauProbeAgainstElectronVLooseMVA6 = 0;
    _tauProbeAgainstElectronLooseMVA6 = 0;
    _tauProbeAgainstElectronMediumMVA6 = 0;
    _tauProbeAgainstElectronTightMVA6 = 0;
    _tauProbeAgainstElectronVTightMVA6 = 0;

    _eleTagPt = -1.;
    _eleTagEta = -1.;
    _eleTagPhi = -1.;
    _eleTagCharge = 0;
    _Mee = 0;
    _isTagHLTmatched = false;
    _isProbeHLTmatched = false;    
    _hltPt.assign(NUMBER_OF_MAXIMUM_TRIGGERS,-1);
    _hltEta.assign(NUMBER_OF_MAXIMUM_TRIGGERS,666);
    _hltPhi.assign(NUMBER_OF_MAXIMUM_TRIGGERS,666);
    _l1tPt = -1;
    _l1tEta = 666;
    _l1tPhi = 666;
    _l1tQual = -1;
    _l1tIso = -1;
    _l1tEmuPt = -1;
    _l1tEmuEta = 666;
    _l1tEmuPhi = 666;
    _l1tEmuQual = -1;
    _l1tEmuIso = -1;
    _l1tEmuNTT = -1;
    _l1tEmuTowerIEta = -1;
    _l1tEmuTowerIPhi = -1;
    _l1tEmuRawEt = -1;
    _l1tEmuIsoEt = -1;
    _foundTag = 0;

    for(unsigned int i=0;i<100;i++){
      _hasL1[i] = -1;
      _hasL1_iso[i] = -1;
    }

}


void Ntuplizer::beginJob()
{
    edm::Service<TFileService> fs;
    _tree = fs -> make<TTree>(_treeName.c_str(), _treeName.c_str());

    //Branches
    _tree -> Branch("EventNumber",&_indexevents,"EventNumber/l");
    _tree -> Branch("RunNumber",&_runNumber,"RunNumber/I");
    _tree -> Branch("lumi",&_lumi,"lumi/I");
    _tree -> Branch("eleProbeTriggerBits", &_eleProbeTriggerBits, "eleProbeTriggerBits/l");
    _tree -> Branch("eleTagTriggerBits", &_eleTagTriggerBits, "eleTagTriggerBits/l");

    _tree -> Branch("eleProbePt",  &_eleProbePt,  "eleProbePt/F");
    _tree -> Branch("eleProbeEta", &_eleProbeEta, "eleProbeEta/F");
    _tree -> Branch("eleProbePhi", &_eleProbePhi, "eleProbePhi/F");
    _tree -> Branch("eleProbeSclEt",  &_eleProbeSclEt,  "eleProbeSclEt/F");
    _tree -> Branch("eleProbeCharge",  &_eleProbeCharge,  "eleProbeCharge/I");

    _tree -> Branch("tauProbePt",  &_tauProbePt,  "tauProbePt/F");
    _tree -> Branch("tauProbeEta", &_tauProbeEta, "tauProbeEta/F");
    _tree -> Branch("tauProbePhi", &_tauProbePhi, "tauProbePhi/F");
    _tree -> Branch("tauProbeCharge",  &_tauProbeCharge,  "tauProbeCharge/I");
    _tree -> Branch("tauProbeDM",  &_tauProbeDM,  "tauProbeDM/I");
    _tree -> Branch("tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT",  &_tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT,  "tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT/O");
    _tree -> Branch("tauProbeByLooseIsolationMVArun2v1DBoldDMwLT",  &_tauProbeByLooseIsolationMVArun2v1DBoldDMwLT,  "tauProbeByLooseIsolationMVArun2v1DBoldDMwLT/O");
    _tree -> Branch("tauProbeByMediumIsolationMVArun2v1DBoldDMwLT",  &_tauProbeByMediumIsolationMVArun2v1DBoldDMwLT,  "tauProbeByMediumIsolationMVArun2v1DBoldDMwLT/O");
    _tree -> Branch("tauProbeByTightIsolationMVArun2v1DBoldDMwLT",  &_tauProbeByTightIsolationMVArun2v1DBoldDMwLT,  "tauProbeByTightIsolationMVArun2v1DBoldDMwLT/O");
    _tree -> Branch("tauProbeByVTightIsolationMVArun2v1DBoldDMwLT",  &_tauProbeByVTightIsolationMVArun2v1DBoldDMwLT,  "tauProbeByVTightIsolationMVArun2v1DBoldDMwLT/O");
    _tree -> Branch("tauProbeAgainstMuonLoose3", &_tauProbeAgainstMuonLoose3, "tauProbeAgainstMuonLoose3/F");
    _tree -> Branch("tauProbeAgainstMuonTight3", &_tauProbeAgainstMuonTight3, "tauProbeAgainstMuonTIght3/F");
    _tree -> Branch("tauProbeAgainstElectronVLooseMVA6", &_tauProbeAgainstElectronVLooseMVA6, "tauProbeAgainstElectronVLooseMVA6/F");
    _tree -> Branch("tauProbeAgainstElectronLooseMVA6", &_tauProbeAgainstElectronLooseMVA6, "tauProbeAgainstElectronLooseMVA6/F");
    _tree -> Branch("tauProbeAgainstElectronMediumMVA6", &_tauProbeAgainstElectronMediumMVA6, "tauProbeAgainstElectronMediumMVA6/F");
    _tree -> Branch("tauProbeAgainstElectronTightMVA6", &_tauProbeAgainstElectronTightMVA6, "tauProbeAgainstElectronTightMVA6/F");
    _tree -> Branch("tauProbeAgainstElectronVTightMVA6", &_tauProbeAgainstElectronVTightMVA6, "tauProbeAgainstElectronVTightMVA6/F");


    _tree -> Branch("eleTagPt",  &_eleTagPt,  "eleTagPt/F");
    _tree -> Branch("eleTagEta", &_eleTagEta, "eleTagEta/F");
    _tree -> Branch("eleTagPhi", &_eleTagPhi, "eleTagPhi/F");
    _tree -> Branch("eleTagCharge",  &_eleTagCharge,  "eleTagCharge/I");
    _tree -> Branch("Mee",  &_Mee,  "Mee/F");

    _tree -> Branch("hltPt",  &_hltPt);
    _tree -> Branch("hltEta", &_hltEta);
    _tree -> Branch("hltPhi", &_hltPhi);
    _tree -> Branch("l1tPt",  &_l1tPt,  "l1tPt/F");
    _tree -> Branch("l1tEta", &_l1tEta, "l1tEta/F");
    _tree -> Branch("l1tPhi", &_l1tPhi, "l1tPhi/F");
    _tree -> Branch("l1tQual", &_l1tQual, "l1tQual/I");
    _tree -> Branch("l1tIso", &_l1tIso, "l1tIso/I");
    _tree -> Branch("l1tEmuPt",  &_l1tEmuPt,  "l1tEmuPt/F");
    _tree -> Branch("l1tEmuEta", &_l1tEmuEta, "l1tEmuEta/F");
    _tree -> Branch("l1tEmuPhi", &_l1tEmuPhi, "l1tEmuPhi/F");
    _tree -> Branch("l1tEmuQual", &_l1tEmuQual, "l1tEmuQual/I");
    _tree -> Branch("l1tEmuIso", &_l1tEmuIso, "l1tEmuIso/I");
    _tree -> Branch("l1tEmuNTT", &_l1tEmuNTT, "l1tEmuNTT/I");
    _tree -> Branch("l1tEmuTowerIEta", &_l1tEmuTowerIEta, "l1tEmuTowerIEta/I");
    _tree -> Branch("l1tEmuTowerIPhi", &_l1tEmuTowerIPhi, "l1tEmuTowerIPhi/I");
    _tree -> Branch("l1tEmuRawEt", &_l1tEmuRawEt, "l1tEmuRawEt/I");
    _tree -> Branch("l1tEmuIsoEt", &_l1tEmuIsoEt, "l1tEmuIsoEt/I");


    _tree -> Branch("isTagHLTmatched", &_isTagHLTmatched, "isTagHLTmatched/O");
    _tree -> Branch("isProbeHLTmatched", &_isProbeHLTmatched, "isProbeHLTmatched/O");
    _tree -> Branch("isOS", &_isOS, "isOS/O");
    _tree -> Branch("foundTag", &_foundTag, "foundTag/I");
    _tree -> Branch("Nvtx", &_Nvtx, "Nvtx/I");

    for(unsigned int i=0;i<100;i++){  
      TString name = Form("hasL1_%i",i);
      _tree -> Branch(name,  &_hasL1[i],  name+"/I");
      TString name_iso = Form("hasL1_iso_%i",i);
      _tree -> Branch(name_iso,  &_hasL1_iso[i],  name_iso+"/I");
    }

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
    Initialize();

    _indexevents = iEvent.id().event();
    _runNumber = iEvent.id().run();
    _lumi = iEvent.luminosityBlock();


    // search for the tag in the event
    edm::Handle<edm::View<reco::GsfElectron> > electrons;
    edm::Handle<pat::TauRefVector>  taus;
    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<reco::Vertex> >  vertices;

    iEvent.getByToken(_electronsTag, electrons);
    iEvent.getByToken(_tauTag, taus);
    iEvent.getByToken(_eleLooseIdMapTag, loose_id_decisions);
    iEvent.getByToken(_eleMediumIdMapTag, medium_id_decisions);
    if(_useHLTMatch)
      iEvent.getByToken(_triggerObjects, triggerObjects);
    iEvent.getByToken(_triggerBits, triggerBits);
    iEvent.getByToken(_VtxTag,vertices);

    if(_useGenMatch)
      iEvent.getByToken(_genParticlesTag, genParticles);


    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    for (unsigned int i = 0; i< electrons->size(); ++i){
      
      const auto eleTag = electrons->ptrAt(i);
      int isTagIDMedium = (*medium_id_decisions)[eleTag];
      if(!isTagIDMedium || eleTag->p4().Pt()<30) continue;
      

      for (unsigned int j = 0; j< electrons->size(); ++j){

	if(i==j) continue;	

	const auto eleProbe = electrons->ptrAt(j);
	int isProbeLoose = (*loose_id_decisions)[eleProbe];
	float eleProbeEta = eleProbe->p4().Eta();
	if(!isProbeLoose || (abs(eleProbeEta)>1.4442 && abs(eleProbeEta)<1.566)) continue;

	bool isOS = (eleTag->charge() / eleProbe->charge() < 0) ? true : false;
	if(!isOS) continue;

	float Mee = (eleTag->p4() + eleProbe->p4()).M();
	if(!(Mee>60 && Mee<120)) continue;

	if(_useGenMatch){
	  if(!matchToTruth(eleProbe,genParticles))
	    continue;
	}


	_isOS = isOS;
	_Mee = Mee;

	//! TagAndProbe on HLT eles

	_eleProbeTriggerBitSet.reset();
	_eleTagTriggerBitSet.reset();

	if(_useHLTMatch){

	  for (pat::TriggerObjectStandAlone  obj : *triggerObjects)
	    {

	      //if(!obj.hasCollection("hltEgammaCandidates::HLT")) continue;
	      
	      const float dR_tag = deltaR (*eleTag, obj);
	      if ( dR_tag < 0.3)
		{		  

		  obj.unpackPathNames(names);

		  const edm::TriggerNames::Strings& triggerNames = names.triggerNames();

		  //Looking for the path
		  unsigned int x = 0;
		  bool foundTrigger = false;	

		  for (const tParameterSet& parameter : _parametersTag)
		    {
		      if ((parameter.hltPathIndex >= 0)&&(obj.hasPathName(triggerNames[parameter.hltPathIndex], true, false)))
			{
			  foundTrigger = true;			  
			  //Retrieving filter list for the event
			  const std::vector<std::string>& filters = (parameter.hltFilters1);
			  if (hasFilters(obj, filters))
			    {
			      //std::cout << "#### FOUND ELE WITH HLT PATH " << x << " ####" << std::endl;			     
			      _eleTagTriggerBitSet[x] = true;			      
			    }
			}
		      x++;
		    }
		  if (foundTrigger){
		    _isTagHLTmatched = true;
		    _foundTag++;
		  }
		}
	      
	      
	      const float dR_probe = deltaR (*eleProbe, obj);
	      if ( dR_probe < 0.3)
		{
		  _isProbeHLTmatched = true;
		  
		  obj.unpackPathNames(names);
		  const edm::TriggerNames::Strings& triggerNames = names.triggerNames();
		  //Looking for the path
		  unsigned int x = 0;
		  bool foundTrigger = false;
		  for (const tParameterSet& parameter : _parametersProbe)
		    {
		      if ((parameter.hltPathIndex >= 0)&&(obj.hasPathName(triggerNames[parameter.hltPathIndex], true, false)))
			{
			  foundTrigger = true;			  
			  const std::vector<std::string>& filters = (parameter.hltFilters1);
			  if (hasFilters(obj, filters))
			  {
			    //std::cout << "#### FOUND ELE WITH HLT PATH " << x << " ####" << std::endl;
			    _hltPt[x] = obj.pt();
			    _hltEta[x] = obj.eta();
			    _hltPhi[x] = obj.phi();
			    _eleProbeTriggerBitSet[x] = true;
			  }
		      }
		    x++;
		  }
		if (foundTrigger) _isProbeHLTmatched = true;
	      }

	    }

	  if(!(_isTagHLTmatched)) continue;

	}      


	// Tau matching
	
	for (pat::TauRef  tau : *taus){

	  const float dR = deltaR(*eleProbe,*tau);
	  if( dR<0.3 ){
	    
	    _tauProbePt = tau->pt();
	    _tauProbeEta = tau->eta();
	    _tauProbePhi = tau->phi();
	    _tauProbeCharge = tau->charge();
	    _tauProbeDM = tau->decayMode();
	    _tauProbeByVLooseIsolationMVArun2v1DBoldDMwLT = tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
	    _tauProbeByLooseIsolationMVArun2v1DBoldDMwLT = tau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
	    _tauProbeByMediumIsolationMVArun2v1DBoldDMwLT = tau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
	    _tauProbeByTightIsolationMVArun2v1DBoldDMwLT = tau->tauID("byTightIsolationMVArun2v1DBoldDMwLT");
	    _tauProbeByVTightIsolationMVArun2v1DBoldDMwLT = tau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT");	    
	    _tauProbeAgainstMuonLoose3 = tau->tauID("againstMuonLoose3");
	    _tauProbeAgainstMuonTight3 = tau->tauID("againstMuonTight3");
	    _tauProbeAgainstElectronVLooseMVA6 = tau->tauID("againstElectronVLooseMVA6");
	    _tauProbeAgainstElectronLooseMVA6 = tau->tauID("againstElectronLooseMVA6");
	    _tauProbeAgainstElectronMediumMVA6 = tau->tauID("againstElectronMediumMVA6");
	    _tauProbeAgainstElectronTightMVA6 = tau->tauID("againstElectronTightMVA6");
	    _tauProbeAgainstElectronVTightMVA6 = tau->tauID("againstElectronVTightMVA6");

	  }

	}





	//! TagAndProbe on L1T EG

	edm::Handle< BXVector<l1t::EGamma> >  L1EGHandle;
	iEvent.getByToken(_L1EGTag, L1EGHandle);

	float minDR = 0.3; //Uncomment for new match algo

	for (l1t::EGammaBxCollection::const_iterator bx0EGIt = L1EGHandle->begin(0); bx0EGIt != L1EGHandle->end(0) ; bx0EGIt++)
	  {
	    const float dR = deltaR(*eleProbe, *bx0EGIt);
	    const l1t::EGamma& l1tEG = *bx0EGIt;	   

	    if (dR < minDR) //Uncomment for new match algo
	      {
		minDR = dR; //Uncomment for new match algo
		_l1tPt = l1tEG.pt();
		_l1tEta = l1tEG.eta();
		_l1tPhi = l1tEG.phi();
		_l1tIso = l1tEG.hwIso();
		_l1tQual = l1tEG.hwQual();
	      }
	  }

	for(unsigned int i=0;i<100;i++){
	  _hasL1[i] = (_l1tPt)>=i;
	  _hasL1_iso[i] = ((_l1tIso) && (_l1tPt)>=i);
	}

	edm::Handle< BXVector<l1t::EGamma> >  L1EmuEGHandle;
	try {iEvent.getByToken(_L1EmuEGTag, L1EmuEGHandle);} catch (...) {;}

	if (L1EmuEGHandle.isValid())
	  {
	    minDR = 0.3;
	
	    for (l1t::EGammaBxCollection::const_iterator bx0EmuEGIt = L1EmuEGHandle->begin(0); bx0EmuEGIt != L1EmuEGHandle->end(0) ; bx0EmuEGIt++)
	      {
		const float dR = deltaR(*eleProbe, *bx0EmuEGIt);
		const l1t::EGamma& l1tEmuEG = *bx0EmuEGIt;
		
		//cout<<"Emul EG, pT = "<<l1tEmuEG.pt()<<", eta = "<<l1tEmuEG.eta()<<", phi = "<<l1tEmuEG.phi()<<endl;
		
		if (dR < minDR) //Uncomment for new match algo
		  {
		    minDR = dR; //Uncomment for new match algo
		    _l1tEmuPt        = l1tEmuEG.pt();
		    _l1tEmuEta       = l1tEmuEG.eta();
		    _l1tEmuPhi       = l1tEmuEG.phi();
		    _l1tEmuIso       = l1tEmuEG.hwIso();
		    _l1tEmuNTT       = l1tEmuEG.nTT();
		    _l1tEmuQual      = l1tEmuEG.hwQual();
		    _l1tEmuTowerIEta = l1tEmuEG.towerIEta();
		    _l1tEmuTowerIPhi = l1tEmuEG.towerIPhi();
		    _l1tEmuRawEt     = l1tEmuEG.rawEt();
		    _l1tEmuIsoEt     = l1tEmuEG.isoEt();
		    
		  }
	      }
	  }
	
	_eleProbePt = eleProbe->pt();
	_eleProbeEta = eleProbe->eta();
	_eleProbePhi = eleProbe->phi();
	_eleProbeSclEt = (eleProbe->superCluster()->energy()) / cosh(eleProbe->superCluster()->eta()) ;
	_eleProbeCharge = eleProbe->charge();


	_eleTagPt = eleTag->pt();
	_eleTagEta = eleTag->eta();
	_eleTagPhi = eleTag->phi();
	_eleTagCharge = eleTag->charge();

	_Nvtx = vertices->size();

 
	_eleProbeTriggerBits = _eleProbeTriggerBitSet.to_ulong();
	_eleTagTriggerBits = _eleTagTriggerBitSet.to_ulong();
	//std::cout << "++++++++++ FILL ++++++++++" << std::endl;
	_tree -> Fill();

      }

    }


}




bool Ntuplizer::hasFilters(const pat::TriggerObjectStandAlone&  obj , const std::vector<std::string>& filtersToLookFor) {

    const std::vector<std::string>& eventLabels = obj.filterLabels();
    for (const std::string& filter : filtersToLookFor)
    {
        //Looking for matching filters
        bool found = false;
        for (const std::string& label : eventLabels)
        {
            //if (label == std::string("hltOverlapFilterIsoMu17MediumIsoPFTau40Reg"))
            if (label == filter)
            {

                //std::cout << "#### FOUND FILTER " << label << " == " << filter << " ####" << std::endl;
                found = true;
            }
        }
        if(!found) return false;
    }

    return true;
}







bool Ntuplizer::matchToTruth(const edm::Ptr<reco::GsfElectron> ele, 
			     const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 || particle->pt()<5)
      continue;
    //
    double dRtmp = deltaR( ele->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return false;
  }

  return true;

}





#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(Ntuplizer);

#endif //NTUPLIZER_H
