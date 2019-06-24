from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_pPb8TeV_AOD_Pbp_RW1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiWithAOD_cfg_prod.py'
config.JobType.maxMemoryMB = 2500
config.JobType.outputFiles = ['Chi_c_pPb8TeV_AOD.root']

config.section_("Data")
config.Data.inputDataset = '/PADoubleMuon/PARun2016C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '285479-285832'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/Chi_c/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'