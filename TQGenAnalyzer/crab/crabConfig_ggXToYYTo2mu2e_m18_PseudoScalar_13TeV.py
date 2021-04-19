from CRABClient.UserUtilities import config
config = config()

## Name of output directory ##
config.General.requestName = 'ggXToYYTo2mu2e_m18_PseudoScalar_13TeV'
config.General.workArea    = 'crab_projects'

## Input analyzer pset ## 
config.JobType.pluginName  = 'analysis'
config.JobType.psetName    = 'ConfFile_cfg.py'
#config.JobType.pyCfgParams = ['globalTag=92X_dataRun2_Prompt_v7','useOOTPhotons=True']
config.JobType.inputFiles  = ['lowPtEleReg_2017UL_25112020.db',
                              'pileup_ALL.root',
                             'pileup_2018.root',
                             'pileup_2017.root',
                             'pileup_2016.root',
                             'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',
                             'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
                              'Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt']
config.JobType.allowUndistributedCMSSW = True

## Input Data ##
#config.Data.inputDataset   = ''
config.Data.userInputFiles = open('ggXToYYTo2mu2e_m18_PseudoScalar_13TeV.txt').readlines()
config.Data.unitsPerJob    = 1
config.Data.splitting      = 'FileBased' 

## Where to run ##
#config.Site.whitelist     = ['T1_US_FNAL']

## Output Data ##
config.Data.outputPrimaryDataset = 'ggXToYYTo2mu2e_m18_PseudoScalar_13TeV_OUTPUT'
config.Data.publication   = False
config.Site.storageSite   = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_egamma/soffi/ggXToYYTo2mu2e/'
