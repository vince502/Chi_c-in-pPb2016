from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_pPb8TeV_MC_2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiConAOD.py'
config.JobType.maxMemoryMB = 2500
config.JobType.outputFiles = ['Chi_c_pPb8TeV_AOD.root']

config.section_("Data")
config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO-eb0de96e274499c444c51980f0cf37bd/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/pPb/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
config.Site.blacklist = ['T2_US_Vanderbilt']
