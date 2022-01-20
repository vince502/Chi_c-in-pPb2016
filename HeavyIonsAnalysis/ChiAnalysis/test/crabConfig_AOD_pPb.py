from CRABClient.UserUtilities import config

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_pPb8TeV_AOD_pPb_RW7'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True #Allow SLC7
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiWithAOD_cfg_prod.py'
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 800
config.JobType.outputFiles = ['Chi_c_pPb8TeV.root']

config.section_("Data")
config.Data.inputDataset = '/PADoubleMuon/PARun2016C-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 60
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/HI/Cert_285952-286496_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '285952-286496'
#config.Data.runRange = '285993-285993'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/okukral/Chi_c/%s' % (config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'