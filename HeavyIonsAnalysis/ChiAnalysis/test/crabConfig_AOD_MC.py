from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_pPb8TeV_AOD_MC1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiWithAOD_cfg_MC.py'
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 800
config.JobType.outputFiles = ['Chi_c_pPb8TeV_MC.root']

config.section_("Data")
config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO-eb0de96e274499c444c51980f0cf37bd/USER'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob   = 10 #3 was suggested by dryrun, but turned out to be too small
config.Data.totalUnits   = -1
config.Data.splitting     = 'FileBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/Chi_c/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'