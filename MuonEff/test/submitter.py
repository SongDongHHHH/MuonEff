import os

pythonCfg = 'MuonDetail_cfg.py'
joblst = ['nopu', 'pu35', 'pu140', 'pu200']
datadir = os.environ["CMSSW_BASE"]+'/src/MuonEff/MuonEff/doc/'
for jobName in joblst:
    fileList = datadir+jobName+'.txt'
    createbatch = "create-batch --cfg %s --jobName %s --fileList %s --maxFiles 5"%(pythonCfg, jobName, fileList)
    print createbatch
    #os.system(createbatch)
