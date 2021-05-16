import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.VarParsing as VarParsing
from PhysicsTools.NanoAOD.common_cff import *

process = cms.Process("GenAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.ntuplizer_seq = cms.Sequence()


#1. setting GT
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v4')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            'file:/afs/cern.ch/user/s/soffi/public/10C23D4F-94BD-E811-9588-E0071B7B2320.root'
#                                '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/60000/854B1DC0-2F71-694D-A3F5-8DC1CDE1EF18.root'
#                                    'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/GenProduction/BPH-RunIISummer20UL18MiniAODv2-00008.root'
#                                    '/store/cmst3/group/bpark/BToKmumu_1000Events_MINIAOD.root'
#                                '/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'
#                                'file:/afs/cern.ch/work/m/mcampana/public/Tetraquark/BPH-RunIISummer20UL18MiniAODv2-00008.root'
#                                    'file:/afs/cern.ch/user/m/mcampana/public/MiniAOD_Y.root'
#                                    '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_2.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_3.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_4.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_5.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_6.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_7.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_8.root',
#          '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_step2/privateBPH_2021Apr28/210428_062206/0000/BPH-RunIISummer20UL18MiniAOD-14GeV_9.root'
                                    '/store/user/mcampana/BPH_Production/private_XToYYTo2mu2e_26GeV_pseudoscalar/privateBPH_2021Apr30/210430_141958/0000/BPH-RunIISummer20UL18MiniAOD-26GeV_22.root'

              )
                            )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple_tetraquarks_wLowPt.root")
)


#2. setting input variables
options = VarParsing.VarParsing('analysis')
options.register ('isMC',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,           # string, int, float, bool
                  "Bool isMC")
options.register ('isSignal',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,           # string, int, float, bool
                  "Bool isSignal")
options.register ('isUpsilon',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,           # string, int, float, bool
                  "Bool isSignal")
options.parseArguments()

if options.isMC and options.isSignal :
    print'Sample is MC Signal'
    xsec=1.

if options.isMC and not options.isSignal : print'Sample is MC Background'
if options.isMC == 0 : print'Sample is Data'

if options.isMC and options.isSignal : index = 0
if options.isMC and not options.isSignal : index = -100
if options.isMC ==0 : index==100

year=2018

if options.isUpsilon :
    massRef=9.46
else:
    massRef=3.09

#3. setting up post reco tools from EGM to have all IDs and SS
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2018-Prompt')  



#4. setting json file
if (options.isMC==False):
    print "applying json"
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    #change json below once you run on new data
    JSONfile = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
    process.source.lumisToProcess.extend(myLumis)
#    print myLumis



#5. setting regression
process.GlobalTag.toGet = cms.VPSet(
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalOnly_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalOnly_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_mean"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_mean_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_eb_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_eb_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")),
cms.PSet(record = cms.string("GBRDWrapperRcd"),
         label = cms.untracked.string("gsfElectron_ee_ecalTrk_05To50_sigma"),
         tag = cms.string("gsfElectron_ee_ecalTrk_05To50_sigma_2018V1"),
         connect = cms.string("sqlite_file:lowPtEleReg_2018_02062020_nv.db")))


#from RecoEgamma.EgammaTools.regressionModifierNoPhoReg_cfi import regressionModifier106XUL
from RecoEgamma.EgammaTools.regressionModifier_cfi import regressionModifier106XUL
_gsfRegressionModifier = regressionModifier106XUL.clone(
    modifierName = 'EGRegressionModifierV3',
    rhoTag = 'fixedGridRhoFastjetAll', # this is ...Tmp in 10_2_X
    eleRegs = dict(
        ecalOnlyMean = dict(
            ebLowEtForestName = "gsfElectron_eb_ecalOnly_05To50_mean",
            ebHighEtForestName = "gsfElectron_eb_ecalOnly_05To50_mean",
            eeLowEtForestName = "gsfElectron_ee_ecalOnly_05To50_mean",
            eeHighEtForestName = "gsfElectron_ee_ecalOnly_05To50_mean",
            ),
        ecalOnlySigma = dict(
            ebLowEtForestName = "gsfElectron_eb_ecalOnly_05To50_sigma",
            ebHighEtForestName = "gsfElectron_eb_ecalOnly_05To50_sigma",
            eeLowEtForestName = "gsfElectron_ee_ecalOnly_05To50_sigma",
            eeHighEtForestName = "gsfElectron_ee_ecalOnly_05To50_sigma",
            ),
        epComb = dict(
            ecalTrkRegressionConfig = dict(
                ebLowEtForestName = "gsfElectron_eb_ecalTrk_05To50_mean",
                ebHighEtForestName = "gsfElectron_eb_ecalTrk_05To50_mean",
                eeLowEtForestName = "gsfElectron_ee_ecalTrk_05To50_mean",
                eeHighEtForestName = "gsfElectron_ee_ecalTrk_05To50_mean",
                ),
            ecalTrkRegressionUncertConfig = dict(
                ebLowEtForestName = "gsfElectron_eb_ecalTrk_05To50_sigma",
                ebHighEtForestName = "gsfElectron_eb_ecalTrk_05To50_sigma",
                eeLowEtForestName = "gsfElectron_ee_ecalTrk_05To50_sigma",
                eeHighEtForestName = "gsfElectron_ee_ecalTrk_05To50_sigma",
                ),
        ),
    ),
)

process.regressionForEle = cms.EDProducer(
    'ElectronRegresser',
    pfSrc = cms.InputTag('slimmedElectrons'),
    gsfRegressionConfig = _gsfRegressionModifier,
)

#5b retrained ID for PF low pt electrons

mvaConfigsForEleProducer = cms.VPSet( )
# Import and add all desired MVAs
from mvaElectronID_BParkRetrain_cff \
    import mvaEleID_BParkRetrain_producer_config
mvaConfigsForEleProducer.append( mvaEleID_BParkRetrain_producer_config )

process.electronMVAValueMapProducer = cms.EDProducer(
  'ElectronMVAValueMapProducer',
  # AOD case
  src = cms.InputTag('gedGsfElectrons'),  
  # miniAOD case
  #srcMiniAOD = cms.InputTag('slimmedElectrons',processName=cms.InputTag.skipCurrentProcess()),
  srcMiniAOD = cms.InputTag('regressionForEle:regressedElectrons'),

  # MVA configurations
  mvaConfigurations = mvaConfigsForEleProducer
)

process.egmGsfElectronIDs = cms.EDProducer(
    "VersionedGsfElectronIdProducer",
    physicsObjectSrc = cms.InputTag('gedGsfElectrons'),
    physicsObjectIDs = cms.VPSet( )
)

process.egmGsfElectronIDTask = cms.Task(
    process.electronMVAValueMapProducer,
    process.egmGsfElectronIDs,
)

process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDTask)


#6. setting Analyzer
process.GenAnalysis = cms.EDAnalyzer('TQGenAnalyzer',
                                     generatorInfo= cms.InputTag("generator"),
                                     prunedGenParticles    = cms.InputTag("prunedGenParticles"),
                                     patMuons = cms.InputTag("slimmedMuons"),
                                     pfSrc = cms.InputTag("regressionForEle:regressedElectrons"),                 # MINIAOD
                                     lowptSrc = cms.InputTag("slimmedLowPtElectrons"), #change when running on our old signals
                                     vtx=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     rho= cms.InputTag('fixedGridRhoAll'),
                                     PileUp = cms.InputTag('slimmedAddPileupInfo'),
                                     bits         = cms.InputTag("TriggerResults::HLT"),
                                     flags        = cms.InputTag("TriggerResults::SIM"),
                                     year  = cms.untracked.int32(year),
                                     sampleIndex  = cms.untracked.int32(index),
                                     sampleXsec  = cms.untracked.double(xsec),
                                     massRef  = cms.untracked.double(massRef),
                                     drForCleaning = cms.double(0.03),
                                     dzForCleaning = cms.double(0.5), ##keep tighter dZ to check overlap of pfEle with lowPt (?)
                                     mvaValuePF = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2BParkRetrainRawValues'),



                                )


#8. setting full path
process.electronsBParkSequence = cms.Sequence(
  process.regressionForEle+
  process.egmGsfElectronIDSequence
)
#process.Seq = cms.Sequence(process.electronsBParkSequence+process.GenAnalysis)
#process.p = cms.Path((process.egammaPostRecoSeq+process.electronsBParkSequence) + process.GenAnalysis)
process.p = cms.Path((process.electronsBParkSequence) + process.GenAnalysis)

