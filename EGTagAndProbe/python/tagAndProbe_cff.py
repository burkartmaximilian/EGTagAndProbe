import FWCore.ParameterSet.Config as cms

print "Running on data"

# filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt



HLTLISTTAG = cms.VPSet(
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle35noerWPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
)


HLTLISTPROBE = cms.VPSet(
    cms.PSet (
        HLT = cms.string("HLT_Ele32_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle32WPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_Ele35_WPTight_Gsf_v"),
        path1 = cms.vstring ("hltEle35noerWPTightGsfTrackIsoFilter"), #FIXME: to check
        path2 = cms.vstring (""),
        leg1 = cms.int32(11),
        leg2 = cms.int32(999)
        ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
)



# filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = ['HLT_Ele32_WPTight_Gsf_v*','HLT_Ele35_WPTight_Gsf_v*'],
    andOr = cms.bool(True), # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(True) #if True: throws exception if a trigger path is invalid
)


## good taus - apply analysis selection
goodTaus = cms.EDFilter("PATTauRefSelector",
        src = cms.InputTag("slimmedTaus"),
        cut = cms.string(
        #        'pt > 5 && abs(eta) < 2.1 ' #kinematics
                'pt > 20 && abs(eta) < 2.1 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 ' #sometimes 2 prongs have charge != 1
                '&& tauID("decayModeFinding") > 0.5 ' # tau ID
                '&& tauID("byMediumIsolationMVArun2v1DBoldDMwLT") > 0.5 '
                '&& tauID("againstMuonTight3") > 0.5 ' # anti Muon tight
        ),
        filter = cms.bool(True)
)


patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
                                    triggerResults = cms.InputTag('TriggerResults', '', "HLT"),
                                    unpackFilterLabels = cms.bool(True)
                                    )


Ntuplizer = cms.EDAnalyzer("Ntuplizer",
    treeName = cms.string("TagAndProbe"),
    electrons = cms.InputTag("gedGsfElectrons"),
    taus = cms.InputTag("goodTaus"),
    genParticles = cms.InputTag("genParticles"),                       
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
    triggerSet = cms.InputTag("patTriggerUnpacker"),
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),   
    L1EG = cms.InputTag("caloStage2Digis", "EGamma", "RECO"),
    L1EmuEG = cms.InputTag("simCaloStage2Digis", "MP"),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    triggerListTag = HLTLISTTAG,
    triggerListProbe = HLTLISTPROBE,
    useGenMatch = cms.bool(False),
    useHLTMatch = cms.bool(True)
)

NtupleSeq = cms.Sequence(
    hltFilter        +
    goodTaus         +
    patTriggerUnpacker +
    Ntuplizer
)
