from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.section_("General")
config.General.requestName = 'Chi_c_XeXe'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runChiWithAOD_cfg_XeXe.py'
config.JobType.maxMemoryMB = 2500
config.JobType.outputFiles = ['Chi_c_XeXe.root']

config.section_("Data")
config.Data.inputDataset = ['/HIMinimumBias/XeXeRun2017-13Dec2017-v1/AOD','/HIMinimumBias1/XeXeRun2017-13Dec2017-v1/AOD','/HIMinimumBias2/XeXeRun2017-13Dec2017-v1/AOD','/HIMinimumBias3/XeXeRun2017-13Dec2017-v1/AOD','/HIMinimumBias4/XeXeRun2017-13Dec2017-v1/AOD']
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/HI/Cert_304899-304907_5TeV_PromptReco_XeXe_Collisions17_JSON.txt'
config.Data.runRange = '285952-286496'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/%s/Chi_c/%s' % (getUsernameFromSiteDB(), config.General.requestName)
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
