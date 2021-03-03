from CRABClient.UserUtilities import config, getUsernameFromCRIC

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_pPb8TeV_AOD_MC6v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.allowUndistributedCMSSW = True #Allow SLC7
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiWithAOD_cfg_MC.py'
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 800
config.JobType.outputFiles = ['Chi_c_pPb8TeV_MC.root']

config.section_("Data")
# config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO-eb0de96e274499c444c51980f0cf37bd/USER' #v1
# config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v2-61416874db099c53202c8cb2d81ec4a3/USER' #v2
# config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v3-61416874db099c53202c8cb2d81ec4a3/USER' #v3
# config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v4-61416874db099c53202c8cb2d81ec4a3/USER' #v4
#config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v5-61416874db099c53202c8cb2d81ec4a3/USER' #v5
config.Data.inputDataset = '/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v6-61416874db099c53202c8cb2d81ec4a3/USER' #v6
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob   = 5 #10 #3 was suggested by dryrun, but turned out to be too small
config.Data.totalUnits   = -1
config.Data.splitting     = 'FileBased'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/Chi_c/%s' % (getUsernameFromCRIC(), config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'