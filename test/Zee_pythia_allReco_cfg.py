import FWCore.ParameterSet.Config as cms

process = cms.Process('TrackingTruthTest')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/Geometry_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/Simulation_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/VtxSmearedGauss_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')


process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.6 $'),
    annotation = cms.untracked.string('RelVal Z to ee'),
    name = cms.untracked.string('Source: /cvs_server/repositories/CMSSW/CMSSW/SimGeneral/TrackingAnalysis/test/Zee_pythia_allReco.cfg,v $')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

process.source = cms.Source("PythiaSource",
  pythiaHepMCVerbosity = cms.untracked.bool(False),
  PythiaParameters = cms.PSet(
    parameterSets = cms.vstring(
      "pythiaZee"
    ),
    pythiaZee = cms.vstring(
      "MSEL = 11 ",           
      "MDME( 174,1) = 0            !Z decay into d dbar",
      "MDME( 175,1) = 0            !Z decay into u ubar",
      "MDME( 176,1) = 0            !Z decay into s sbar",
      "MDME( 177,1) = 0            !Z decay into c cbar",
      "MDME( 178,1) = 0            !Z decay into b bbar",
      "MDME( 179,1) = 0            !Z decay into t tbar",
      "MDME( 182,1) = 1            !Z decay into e- e+",
      "MDME( 183,1) = 0            !Z decay into nu_e nu_ebar",
      "MDME( 184,1) = 0            !Z decay into mu- mu+",
      "MDME( 185,1) = 0            !Z decay into nu_mu nu_mubar",
      "MDME( 186,1) = 0            !Z decay into tau- tau+",
      "MDME( 187,1) = 0            !Z decay into nu_tau nu_taubar",
      "MSTJ( 11) = 3               !Choice of the fragmentation function",
      "MSTP( 2) = 1                !which order running alphaS",
      "MSTP( 33) = 0               !(D=0) ",
      "MSTP( 51) = 7               !structure function chosen",
      "MSTP( 81) = 1               !multiple parton interactions 1 is Pythia default",
      "MSTP( 82) = 4               !Defines the multi-parton model",
      "PARJ( 71) = 10.             !for which ctau  10 mm",
      "PARP( 82) = 1.9             !pt cutoff for multiparton interactions",
      "PARP( 89) = 1000.           !sqrts for which PARP82 is set",
      "PARP( 83) = 0.5             !Multiple interactions: matter distrbn parameter Registered by Chris.Seez@cern.ch",
      "PARP( 84) = 0.4             !Multiple interactions: matter distribution parameter Registered by Chris.Seez@cern.ch",
      "PARP( 90) = 0.16            !Multiple interactions: rescaling power Registered by Chris.Seez@cern.ch",
      "CKIN( 1) = 40.              !(D=2. GeV)",
      "CKIN( 2) = -1.              !(D=-1. GeV)"
    )
  )
)

# Output definition
process.output = cms.OutputModule(
  "PoolOutputModule",
  fileName = cms.untracked.string('Zee_pythia_withTkElectrons.root'),
  outputCommands = cms.untracked.vstring(
    "keep *",
    "drop *_mergedtruth__*"
  )
)

# Path and EndPath definitions
process.p1 = cms.Path(process.psim)
process.p2 = cms.Path(process.pdigi)
process.p3 = cms.Path(process.globalreco_plusRS_plusGSF)
process.outpath = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.p1,process.p2,process.outpath)

