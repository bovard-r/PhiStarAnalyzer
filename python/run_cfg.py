import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("PhiStarAnalyzer")

options = VarParsing.VarParsing('analysis')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("analysis.PhiStarAnalyzer.phistaranalyzer_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

options.register('fileIdx', 0,
        VarParsing.VarParsing.multiplicity.singleton,
        VarParsing.VarParsing.varType.int,
        "Index of the filelist chunk to run over")

options.parseArguments()

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        phiStarAnalyzer = cms.PSet(
            initialSeed = cms.untracked.uint32(options.fileIdx + 5150),
            engineName = cms.untracked.string('TRandom3')
        )
)

with open(f"analysis/PhiStarAnalyzer/python/filelist_{options.fileIdx}.txt") as f:
    file_list = [line.strip() for line in f if line.strip()]

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(file_list),
        duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
        )

process.phiStarAnalyzer = cms.EDAnalyzer("PhiStarAnalyzer")

process.phiStarAnalyzer.genParticles = cms.InputTag('genParticles')
process.phiStarAnalyzer.ak4GenJetsNoNu = cms.InputTag('ak4GenJetsNoNu')
process.phiStarAnalyzer.lheInfo = cms.InputTag('externalLHEProducer')
process.phiStarAnalyzer.genEventInfo = cms.InputTag('generator')
process.phiStarAnalyzer.variations = process.variations
for idx in range(len(process.variations)):
    process.phiStarAnalyzer.variations[idx].outputDir = cms.string(process.phiStarAnalyzer.variations[idx].outputDir.value() + str(options.fileIdx))

process.p = cms.Path(process.phiStarAnalyzer)
