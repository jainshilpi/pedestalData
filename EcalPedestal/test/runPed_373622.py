import FWCore.ParameterSet.Config as cms

process = cms.Process("Noise")

process.load("FWCore.MessageService.test.Services_cff")

process.MessageLogger = cms.Service("MessageLogger",
 # destinations = cms.untracked.vstring('messages.txt'),
)

#######################################################################################
# input
### Case 1:  PoolSource (from DBS, castor, etc.)
### Case 2:  Run on raw .dat files
process.source = cms.Source("NewEventStreamFileReader",
#  skipEvents = cms.untracked.uint32(1),
#  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),   # only with PoolSource!
  fileNames =  cms.untracked.vstring(
  ##  'file:run366748_ls0010_streamDQM_pid49016.dat'
      'root://eoscms//eos/cms/store/group/dpg_ecal/comm_ecal/fromP5/run373622/run373622_ls0010_streamDQM_pid12996.dat'
  )
)

'''
process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound')
)
'''

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(100)
)
#######################################################################################

 
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
 
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
# process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

# Conditions (Global Tag is used here):
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.prefer("GlobalTag")
# process.GlobalTag.globaltag = 'GR_P_V56'
process.GlobalTag.globaltag = '120X_dataRun3_HLT_v3'

process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")

process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi")
# unpack raw data from global run
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi")
process.ecalEBunpacker.InputLabel = cms.InputTag('rawDataCollector')

process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("RecoLocalCalo.EcalRecProducers.ecalRatioUncalibRecHit_cfi")
process.load("RecoLocalCalo.EcalRecProducers.ecalGlobalUncalibRecHit_cfi")
process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff")
#process.load("RecoLocalCalo.EcalRecProducers.ecalDetIdToBeRecovered_cfi")
#
# UncalibRecHit producer
import RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi
process.ecalWeightUncalibRecHit = RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi.ecalWeightUncalibRecHit.clone()
import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
process.ecalFixedAlphaBetaFitUncalibRecHit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()
import RecoLocalCalo.EcalRecProducers.ecalMaxSampleUncalibRecHit_cfi
process.ecalMaxSampleUncalibRecHit = RecoLocalCalo.EcalRecProducers.ecalMaxSampleUncalibRecHit_cfi.ecalMaxSampleUncalibRecHit.clone()
process.ecalWeightUncalibRecHit.EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")
process.ecalWeightUncalibRecHit.EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")
process.ecalFixedAlphaBetaFitUncalibRecHit.EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")
process.ecalFixedAlphaBetaFitUncalibRecHit.EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")


process.ecalGlobalUncalibRecHit.EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")
process.ecalGlobalUncalibRecHit.EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")
process.ecalMultiFitUncalibRecHit.cpu.EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")
process.ecalMultiFitUncalibRecHit.cpu.EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")

process.TFileService = cms.Service("TFileService", fileName = cms.string('timeSampleTree_ped.root'))



process.Pedestal = cms.EDAnalyzer("EcalPedestal",
  digiProducer = cms.string("ecalEBunpacker"),
  EBdigiCollection = cms.string("ebDigis"),
  EEdigiCollection = cms.string("eeDigis"),

  ebRecHitCollection = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
  eeRecHitCollection = cms.InputTag("ecalRecHit:EcalRecHitsEE"),

  EEuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit:EcalUncalibRecHitsEE"),
  EBuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit:EcalUncalibRecHitsEB"),

  digiCollection = cms.string("ecalEBunpacker"),

  hitProducer = cms.string("ecal2006TBWeightUncalibRecHit"),
  hitCollection = cms.string("EcalUncalibRecHitsEB"),

  eventHeaderProducer = cms.string("ecalEBunpacker"),
  eventHeaderCollection = cms.string(""),

  runnumber = cms.untracked.int32(373622),
#  ECALType = cms.string("EA"),
  ECALType = cms.string("EB"),
  runType = cms.string("Pedes"),
  startevent = cms.untracked.uint32(1),
  xtalnumber = cms.untracked.int32(0),
  readPedestals = cms.bool(False),
  patternNoise = cms.bool(False),
  SampleCorr = cms.bool(False),
  NormalSequence = cms.bool(False),
#  b904 = cms.bool(False),
  gain12events = cms.untracked.int32(300),
  gain6events = cms.untracked.int32(600),
  gain1events = cms.untracked.int32(900),

                                  
)

'''
process.InterimOutput = cms.OutputModule("PoolOutputModule",
                                         fileName = cms.untracked.string('myOutputFile.root'), 
                                          SelectEvents = cms.untracked.PSet(
                                              SelectEvents = cms.vstring("p")
                                              ),
                                         outputCommands = cms.untracked.vstring('keep *')
                                         
                                     )
'''

process.ecalDetIdToBeRecovered.ebSrFlagCollection = cms.InputTag("ecalEBunpacker")
process.ecalDetIdToBeRecovered.eeSrFlagCollection = cms.InputTag("ecalEBunpacker")
process.ecalDetIdToBeRecovered.ebIntegrityChIdErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityChIdErrors")
process.ecalDetIdToBeRecovered.ebIntegrityGainErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityGainErrors")
process.ecalDetIdToBeRecovered.ebIntegrityGainSwitchErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityGainSwitchErrors")


process.ecalDetIdToBeRecovered.eeIntegrityChIdErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityChIdErrors")
process.ecalDetIdToBeRecovered.eeIntegrityGainErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityGainErrors")
process.ecalDetIdToBeRecovered.eeIntegrityGainSwitchErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityGainSwitchErrors")


process.ecalDetIdToBeRecovered.integrityTTIdErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityTTIdErrors")
process.ecalDetIdToBeRecovered.integrityBlockSizeErrors = cms.InputTag("ecalEBunpacker","EcalIntegrityBlockSizeErrors")


process.p = cms.Path(process.ecalEBunpacker*
                     process.ecalWeightUncalibRecHit*
                     process.ecalFixedAlphaBetaFitUncalibRecHit*
                     process.ecalMaxSampleUncalibRecHit*
                     process.bunchSpacingProducer *
                     process.ecalLocalRecoSequence *
                     process.Pedestal)


#process.e = cms.EndPath(process.InterimOutput)  

print(process.dumpPython())
