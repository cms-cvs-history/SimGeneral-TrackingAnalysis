import FWCore.ParameterSet.Config as cms

process = cms.Process('TrackingTruthPlayback')

# TrackingTruth
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")

# Output definition
process.output = cms.OutputModule(
  'PoolOutputModule',
  fileName = cms.untracked.string('TrackingTruthPlayback.root'),
  outputCommands = cms.untracked.vstring(
    'keep edmHepMCProduct_source_*_*',
    'keep *_mergedtruth__*',
    'keep *_mergedtruth_*_*'
  )
)

process.path = cms.Path(process.trackingParticles)
process.outpath = cms.EndPath(process.output)

# Input definition
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
readFiles = cms.untracked.vstring('rfio:/castor/cern.ch/user/b/benedet/Data219SingleEle/Single_e_pt4_nosmear_5.root')
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

