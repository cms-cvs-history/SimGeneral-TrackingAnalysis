import FWCore.ParameterSet.Config as cms

process = cms.Process('TrackingTruthTest')

# Playback
from Configuration.StandardSequences.Services_cff import *
del RandomNumberGeneratorService.theSource
RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')
from SimGeneral.MixingModule.mixNoPU_cfi import *
mix.playback = cms.untracked.bool(True)

