#ifndef ECALPEDESTAL_H
#define ECALPEDESTAL_H

/**\class EcalPedestal

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// $Id: EcalPedestal.h,v 0.0 2006/11/17 19:24:46 meridian Exp $
//



// system include files

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "CaloOnlineTools/EcalTools/interface/EcalFedMap.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <TTree.h>

#include <string>
//#include "TTree.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include<fstream>
#include<map>
//#include<stl_pair>

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

using namespace edm;
using namespace cms;
using namespace std;

namespace {
  class EcalPedestal : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
  public:
    explicit EcalPedestal( const edm::ParameterSet& );
    ~EcalPedestal();


    void beginJob() override;
    void beginRun( const edm::Run&, const edm::EventSetup&) override;
    void analyze( const edm::Event&, const edm::EventSetup& ) override;
    void endRun( const edm::Run&, const edm::EventSetup&) override;
    void endJob() override;
    enum { kChannels = 61200,kEBChannels = 61200, kEEChannels = 14648};
    enum { kGains = 3, kFirstGainId = 1};
  private:
    //  const EcalElectronicsMapping ecalElectronicsMap_;
    

    const edm::ESGetToken<EcalElectronicsMapping, EcalMappingRcd> EcalElectronicsMappingToken;

        ///SJ
    edm::InputTag ebcalibRecHitCollectionTag_; // secondary name given to collection of EB uncalib rechits
    edm::InputTag eecalibRecHitCollectionTag_; // secondary name
    edm::EDGetTokenT<EcalRecHitCollection>           ebRecHitCollection_;
    edm::EDGetTokenT<EcalRecHitCollection>           eeRecHitCollection_;
    ///SJ

    const edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> ChannelStatusToken;
    EcalElectronicsId elecId[kEEChannels];
    EcalFedMap* fedMap_;
    edm::InputTag EcalTBEventHeader_;
    edm::InputTag EcalRawDataCollection_;
    //  edm::InputTag EBDigiCollection_;
    edm::InputTag EBuncalibRecHitCollection_; // secondary name given to collection of EB uncalib rechits
    edm::InputTag EEuncalibRecHitCollection_; // secondary name given to collection of EE uncalib rechits
    edm::InputTag EBDigisTag = edm::InputTag("ecalEBunpacker", "ebDigis");
    edm::InputTag EEDigisTag = edm::InputTag("ecalEBunpacker", "eeDigis");
    edm::EDGetTokenT<EBDigiCollection> EBDigisToken;
    edm::EDGetTokenT<EEDigiCollection> EEDigisToken;
    edm::EDGetTokenT<EBRecHitCollection> EBUncalibRecHitToken;
    edm::EDGetTokenT<EBRecHitCollection> EBUncalibRecHitMaxToken;
    edm::EDGetTokenT<EBRecHitCollection> EBUncalibRecHitFitToken;
    edm::EDGetTokenT<EERecHitCollection> EEUncalibRecHitToken;
    std::string rootfile_;
    std::string hitCollection_;
    std::string hitProducer_;
    std::string tdcRecInfoCollection_;
    std::string tdcRecInfoProducer_;
    std::string eventHeaderCollection_;
    std::string eventHeaderProducer_;
    int runnumber_;
    int firstReadEvent_;
    std::string ECALType_; // EB or EE
    std::string runType_; // Pedes or Other
    int gain12events_;
    int gain6events_;
    int gain1events_;
    unsigned int startevent_; 
    int xtalnumber_;
    bool readPedestals_;
    bool patternNoise_;
    bool SampleCorr_;
    bool NormalSequence_;
    unsigned int cnt_evt_;
    unsigned int kEventGains;
    int kEta;
    int maxIdRef_;
    int nevent[10];
    int iEvent_;
    Timestamp begtime[3], endtime[3];
    std::vector<int> maskedChannels_;
    std::vector<int> maskedEEChannels_;
    std::vector<int> badChannels_;
    std::vector<int> badEEChannels_;
    std::ofstream fLargeAmp;

    std::vector<int>    nEntriesStandard_;
    std::vector<double> pedestalStandard_;
    std::vector<double> pedestal2Standard_;
    std::vector<double> PedestalRead_;

    std::vector<int>    nEntriesHit_;
    std::vector<double> pedestalHit_;
    std::vector<double> pedestal2Hit_;

    std::vector<int>    nEntriesMinMax_;
    std::vector<double> pedestalMinMax_;
    std::vector<double> pedestal2MinMax_;

    std::vector<double> avgSamples_;
    std::vector<double> avgSample2_;
    std::vector<int>    nEntriesAvg_;
    std::vector<double> pedavg_;
    std::vector<double> pedavg2_;

    std::vector<int>    nEntriesSingle_;
    std::vector<double> pedestalSingle_;
    std::vector<double> pedestal2Single_;
    std::vector<double> pedestal_5_;
    std::vector<double> pedestal2_5_;
    std::vector<double> pedestal_3p5_;
    std::vector<double> pedestal2_3p5_;
    std::vector<double> pedestal_2p1_;
    std::vector<double> pedestal2_2p1_;
    std::vector<double> pedestal_3p1_;
    std::vector<double> pedestal2_3p1_;

    std::vector<double> pedestal_3g1_;
    std::vector<double> pedestal2_3g1_;
    std::vector<double> pedestal_3g2_;
    std::vector<double> pedestal2_3g2_;
    std::vector<double> pedestal_3g3_;
    std::vector<double> pedestal2_3g3_;
    std::vector<double> pedestal_3g4_;
    std::vector<double> pedestal2_3g4_;
    std::vector<double> pedestal_3g5_;
    std::vector<double> pedestal2_3g5_;
    std::vector<double> pedestal_3g6_;
    std::vector<double> pedestal2_3g6_;
    std::vector<double> pedestal_3g7_;
    std::vector<double> pedestal2_3g7_;

    std::vector<double> pedestal_4p1_;
    std::vector<double> pedestal2_4p1_;
    std::vector<double> pedestal_5p1_;
    std::vector<double> pedestal2_5p1_;
    std::vector<double> pedestal_6p1_;
    std::vector<double> pedestal2_6p1_;
    std::vector<double> pedestal_7p1_;
    std::vector<double> pedestal2_7p1_;
    std::vector<double> pedestal_8p1_;
    std::vector<double> pedestal2_8p1_;
    std::vector<double> pedestal_9p1_;
    std::vector<double> pedestal2_9p1_;
    std::vector<double> noise_1_;
    std::vector<double> noise_25_;
    std::vector<double> noisesum_25_;
    std::vector<int>    nXtals_25_;
    std::vector<double> noise_9_;
    std::vector<double> noisesum_9_;
    std::vector<int>    nXtals_9_;

    std::vector<double> noiseSingle_25_;
    //  std::vector<double> noisesumSingle_25_;
    std::vector<double> noiseSingle_9_;
    //  std::vector<double> noisesumSingle_9_;

    std::vector<double> noise_5_25_;
    std::vector<double> noise_5_9_;

    std::vector<double> noise_3p5_;
    std::vector<double> noise_3p5_25_;
    std::vector<double> noise_3p5_9_;

    std::vector<double> noise_2p1_25_;
    std::vector<double> noise_2p1_9_;

    std::vector<double> noise_3p1_25_;
    std::vector<double> noise_3p1_9_;

    std::vector<double> noise_4p1_25_;
    std::vector<double> noise_4p1_9_;

    std::vector<double> noise_5p1_25_;
    std::vector<double> noise_5p1_9_;

    std::vector<double> noise_6p1_25_;
    std::vector<double> noise_6p1_9_;

    std::vector<double> noise_7p1_25_;
    std::vector<double> noise_7p1_9_;

    std::vector<double> noise_8p1_25_;
    std::vector<double> noise_8p1_9_;

    std::vector<double> noise_9p1_25_;
    std::vector<double> noise_9p1_9_;
 
    // Reconstructed energies
    std::map< int, std::map< int, TH1F*> > hPedestal;
    std::map< int, std::map< int, TH1F*> > hPedestal_9;
    std::map< int, std::map< int, TH1F*> > hPedestal_25;

    // endcaps
    std::vector<int>    nEntriesStandardEE_;
    std::vector<double> pedestalStandardEE_;
    std::vector<double> pedestal2StandardEE_;
    std::vector<double> PedestalReadEE_;

    std::vector<double> avgSamplesEE_;
    std::vector<double> avgSample2EE_;
    std::vector<int>    nEntriesAvgEE_;
    std::vector<double> pedavgEE_;
    std::vector<double> pedavg2EE_;

    std::vector<int>    nEntriesHitEE_;
    std::vector<double> pedestalHitEE_;
    std::vector<double> pedestal2HitEE_;

    float meanEBped[kGains][kEBChannels][1000], meanEEped[kGains][kEEChannels][1000];

    //SJ
    const edm::ESGetToken<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd> tok_ecalPFRecHitThresholds_;
    
    std::vector<std::vector<double>> hitsAmplitudes_;
    std::vector<double> hitsEnergy_;
    std::vector<double> hitsThr_;
    
    int nCrys_;
    TTree   *tree;



  };
}
#endif
