/**\class EcalPedestal
*/
//
// $Id: EcalEcalPedestal.cc,v 0.0 2006/11/17 19:01:00 fay Exp $
//

//#include "Pit/Pedestal/interface/EcalPedestal.h"
//#include "Pit/H4CovMat/interface/H4CovMat.h"

#include "pedestalData/EcalPedestal/interface/EcalPedestal.h"
#include "pedestalData/H4CovMat/interface/H4CovMat.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCRecInfo.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBEventHeader.h"
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCSample.h" 
#include "TBDataFormats/EcalTBObjects/interface/EcalTBTDCRawInfo.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include <DataFormats/EcalRawData/interface/EcalRawDataCollections.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"


//#include<fstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"

//
// constants, enums and typedefs
//
const Int_t kSample=10;
TProfile* hEntries;
TProfile** hPedsample = new TProfile*[3];
TH1F* hGains;
//
// static data member definitions
//
int SMCal[36] = {24,22,13,31,26,16, 2,11, 5, 0,29,28,14,33,32, 3, 7,19,
		 12,17,10, 1, 8, 4,27,20,23,25, 6,34,35,15,18,30,21, 9};

vector<H4CovMat*>* vCovMatEB;
vector<H4CovMat*>* vCovMatEE;

int gainValues[3] = {12, 6, 1};
TH1F** hMaxMin = new TH1F*[3];
TH2F** hEBEntry = new TH2F*[3];
TProfile** hGain = new TProfile*[3];
TProfile** hGainEB = new TProfile*[3];
TProfile** hGainEE = new TProfile*[3];
TH2F** hEBDiff = new TH2F*[3];
TH2F*** hEEDiff = new TH2F**[3];

//
// constructors and destructor
//

//===================================================================
EcalPedestal::EcalPedestal( const edm::ParameterSet& iConfig )
    : EcalElectronicsMappingToken(esConsumes()),
      //ebRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRecHitCollection"))),
      //eeRecHitCollection_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRecHitCollection"))),
      ebRecHitCollection_(consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit", "EcalRecHitsEB"))),
      eeRecHitCollection_(consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit", "EcalRecHitsEE"))),
      geometryToken_(esConsumes()), 
      ChannelStatusToken(esConsumes()),
      tok_ecalPFRecHitThresholds_(esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>())


//====================================================================
{
   //now do what ever initialization is needed   
  ///SJ
  //ebcalibRecHitCollectionTag_ = iConfig.getParameter<edm::InputTag>("ebRecHitCollection");
  //eecalibRecHitCollectionTag_ = iConfig.getParameter<edm::InputTag>("eeRecHitCollection");
  //ebRecHitCollection_ = EDGetTokenT<EcalRecHitCollection>(consumes<EcalRecHitCollection>(ebcalibRecHitCollectionTag_));
  //eeRecHitCollection_ = EDGetTokenT<EcalRecHitCollection>(consumes<EcalRecHitCollection>(eecalibRecHitCollectionTag_));

  
  hitCollection_             = iConfig.getParameter<std::string>("hitCollection");
  EEuncalibRecHitCollection_ = iConfig.getParameter<edm::InputTag>("EEuncalibRecHitCollection");
  EBuncalibRecHitCollection_ = iConfig.getParameter<edm::InputTag>("EBuncalibRecHitCollection");
  hitProducer_               = iConfig.getParameter<std::string>("hitProducer");
//  tdcRecInfoCollection_      = iConfig.getParameter<std::string>("tdcRecInfoCollection");
//  tdcRecInfoProducer_        = iConfig.getParameter<std::string>("tdcRecInfoProducer");
  eventHeaderCollection_     = iConfig.getParameter<std::string>("eventHeaderCollection");
  eventHeaderProducer_       = iConfig.getParameter<std::string>("eventHeaderProducer");
  runnumber_                 = iConfig.getUntrackedParameter<int>("runnumber",-1);
  ECALType_                  = iConfig.getParameter<std::string>("ECALType");
  runType_                   = iConfig.getParameter<std::string>("runType");
  startevent_                = iConfig.getUntrackedParameter<unsigned int>("startevent", 1);
  xtalnumber_                = iConfig.getUntrackedParameter<int>("xtalnumber",-1);
  readPedestals_             = iConfig.getParameter<bool>("readPedestals");
  patternNoise_              = iConfig.getParameter<bool>("patternNoise");
  SampleCorr_                = iConfig.getParameter<bool>("SampleCorr");
  NormalSequence_            = iConfig.getParameter<bool>("NormalSequence");
  gain12events_              = iConfig.getUntrackedParameter<int>("gain12events",-1);
  gain6events_               = iConfig.getUntrackedParameter<int>("gain6events",-1);
  gain1events_               = iConfig.getUntrackedParameter<int>("gain1events",-1);

  EBDigisToken = EDGetTokenT<EBDigiCollection>(consumes<EBDigiCollection>(EBDigisTag));
  EEDigisToken = EDGetTokenT<EEDigiCollection>(consumes<EEDigiCollection>(EEDigisTag));

  
  /*  
  EBUncalibRecHitToken = EDGetTokenT<EBRecHitCollection>(consumes<EBDigiCollection>(EBuncalibRecHitCollection_));
  edm::InputTag EBuncalibRecHitMaxCollection = edm::InputTag("ecalMaxSampleUncalibRecHit:EcalUncalibRecHitsEB");
  EBUncalibRecHitMaxToken = EDGetTokenT<EBRecHitCollection>(consumes<EBDigiCollection>(EBuncalibRecHitMaxCollection));
  edm::InputTag EBuncalibRecHitFitCollection = edm::InputTag("ecalFixedAlphaBetaFitUncalibRecHit:EcalUncalibRecHitsEB");
  EBUncalibRecHitFitToken = EDGetTokenT<EBRecHitCollection>(consumes<EBDigiCollection>(EBuncalibRecHitFitCollection));
  EEUncalibRecHitToken = EDGetTokenT<EERecHitCollection>(consumes<EEDigiCollection>(EEuncalibRecHitCollection_));
  */
  vector<int> listDefaults;
  listDefaults.push_back(-1);  

  cnt_evt_ = 0;

   cout << "Exiting constructor" << endl;
}//constructor


//========================================================================
EcalPedestal::~EcalPedestal()
//========================================================================
{
  cout << "ANALYSIS FINISHED" << endl;
}//destructor

//========================================================================
void EcalPedestal::beginRun(edm::Run const & run, edm::EventSetup const & c)
//========================================================================
{
  cout << "Entering beginRun" << endl;
  /*
  //  edm::ESHandle< EcalElectronicsMapping > handle;
  //  const edm::ESGetToken<EcalElectronicsMapping, EcalMappingRcd> EcalElectronicsMappingToken;
  //  c.get< EcalMappingRcd >().get(handle);
  //  ecalElectronicsMap_ = handle.product();
  const auto& ecalElectronicsMap_ = c.getData(EcalElectronicsMappingToken);
  cout << "et un " << endl;
  for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
    EEDetId myEEDetId = EEDetId::unhashIndex(iChannel);
    elecId[iChannel] = ecalElectronicsMap_.getElectronicsId(myEEDetId);
    cout << " channel " << iChannel << " DCC " << elecId[iChannel].dccId() << endl;
  }
  */
  fedMap_ = new EcalFedMap();

  //  edm::ESHandle<EcalChannelStatus> pChannelStatus;
  //  if(auto handle = c.getHandle(
  //  c.get<EcalChannelStatusRcd>().get(pChannelStatus);
  //  const EcalChannelStatus* chStatus = pChannelStatus.product();
  /*
  const auto& chStatus = c.getData(ChannelStatusToken);
  cout << "tr " << endl;
  EcalChannelStatusMap::const_iterator chit;
  cout << "qu " << endl;
  for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
    EBDetId id = EBDetId::unhashIndex(iChannel);
    chit = chStatus.getMap().find(id.rawId());
    if( chit != chStatus.getMap().end() ) {
      EcalChannelStatusCode ch_code = (*chit);
      uint16_t statusCode = ch_code.getStatusCode() & 31;
      //      if(statusCode == 1 || (statusCode > 7 && statusCode < 13))
      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616: status 13 or 14 means bad data!
	maskedChannels_.push_back(iChannel);
      //      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || (statusCode > 7 && statusCode < 13))
      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	badChannels_.push_back(iChannel);
    }
    for (int igain = 0; igain < kGains; igain++) 
      for (int iev = 0; iev < 1000; iev++) 
	meanEBped[igain][iChannel][iev] = -999.;
  }
  // channel 46619 is very noisy but not in the tag (??)
  badChannels_.push_back(46619);
  for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
    EEDetId id = EEDetId::unhashIndex(iChannel);
    chit = chStatus.getMap().find(id.rawId());
    if( chit != chStatus.getMap().end() ) {
      EcalChannelStatusCode ch_code = (*chit);
      uint16_t statusCode = ch_code.getStatusCode() & 31;
      //      if(statusCode == 1 || (statusCode > 7 && statusCode < 13))
      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	maskedEEChannels_.push_back(iChannel);
      //      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || (statusCode > 7 && statusCode < 13))
      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	badEEChannels_.push_back(iChannel);
    }
    for (int igain = 0; igain < kGains; igain++) 
      for (int iev = 0; iev < 1000; iev++) 
	meanEEped[igain][iChannel][iev] = -999.;
  }
  */
  cout << "Exiting beginRun" << endl;
}//beginRun

//========================================================================
void EcalPedestal::endRun(edm::Run const & run, edm::EventSetup const & c) {
//========================================================================
  cout << "Entering endRun" << endl;

  cout << "Exiting endRun" << endl;
}//endRun

//========================================================================
void EcalPedestal::beginJob() {
///========================================================================

  cout << "Entering beginJob" << endl;

  ///SJ
  edm::Service<TFileService> fs;
  tree    = fs->make<TTree>("EventTreeEle", "Event data");
  tree->Branch("hitsAmplitudes",         &hitsAmplitudes_);
  tree->Branch("hitsEnergy",         &hitsEnergy_);
  tree->Branch("hitsThr",         &hitsThr_);
  tree->Branch("hitsEta",         &hitsEta_);
  tree->Branch("hitsPhi",         &hitsPhi_);
  tree->Branch("nCrys",         &nCrys_);
  
  nEntriesStandard_.resize(kEBChannels*kGains, 0);
  pedestalStandard_.resize(kEBChannels*kGains, 0.);
  pedestal2Standard_.resize(kEBChannels*kGains, 0.);

  avgSamples_.resize(kEBChannels*kGains, 0.);
  avgSample2_.resize(kEBChannels*kGains, 0.);
  nEntriesAvg_.resize(kEBChannels*kGains, 0);
  pedavg_.resize(kEBChannels*kGains, 0.);
  pedavg2_.resize(kEBChannels*kGains, 0.);

  nEntriesHit_.resize(kEBChannels*kGains, 0);
  pedestalHit_.resize(kEBChannels*kGains, 0.);
  pedestal2Hit_.resize(kEBChannels*kGains, 0.);

  nEntriesMinMax_.resize(kEBChannels*kGains, 0);
  pedestalMinMax_.resize(kEBChannels*kGains, 0.);
  pedestal2MinMax_.resize(kEBChannels*kGains, 0.);

  nEntriesSingle_.resize(kEBChannels*kGains, 0);
  pedestalSingle_.resize(kEBChannels*kGains, 0.);
  pedestal2Single_.resize(kEBChannels*kGains, 0.);


  // endcaps
  nEntriesStandardEE_.resize(kEEChannels*kGains, 0);
  pedestalStandardEE_.resize(kEEChannels*kGains, 0.);
  pedestal2StandardEE_.resize(kEEChannels*kGains, 0.);

  avgSamplesEE_.resize(kEEChannels*kGains, 0.);
  avgSample2EE_.resize(kEEChannels*kGains, 0.);
  nEntriesAvgEE_.resize(kEEChannels*kGains, 0);
  pedavgEE_.resize(kEEChannels*kGains, 0.);
  pedavg2EE_.resize(kEEChannels*kGains, 0.);

  nEntriesHitEE_.resize(kEEChannels*kGains, 0);
  pedestalHitEE_.resize(kEEChannels*kGains, 0.);
  pedestal2HitEE_.resize(kEEChannels*kGains, 0.);

  hEntries = new TProfile("Entries","entries per SM", 36,1. ,37);
  // Reconstructed energies
  hGains = new TH1F("hGains","Gain vs event #", 1000,0 ,1000.);

  for (int gainId = 0; gainId < kGains; gainId++) {
    hEBEntry[gainId] = new TH2F(Form("EBEntry_%i",gainId),
				Form("Digits gain %i",gainValues[gainId]),360,1.,361.,171,-85.,86.);
    hMaxMin[gainId] = new TH1F(Form("MaxMin_%i",gainId),
			       Form("max - min gain %i",gainValues[gainId]),100, 0., 100.);
    hPedsample[gainId] = new TProfile(Form("Pedsample_%i",gainId),Form("Pedestal samples gain %i",gainValues[gainId]), 
				      10,0 ,10,-10,10);
    int nev = gain12events_;
    if(!NormalSequence_) {
      if(gainId == 1) nev = gain6events_ -  gain12events_;
      else if(gainId == 2)  nev = gain1events_ -  gain6events_;
    }
    float xnev = (float)nev + 1.;
    cout << "hGain for gain " << gainValues[gainId] << " nb of bins " << nev << " last bin " << xnev << endl;
    hGain[gainId] = new TProfile(Form("Gain%i", gainValues[gainId]), Form("gain%i pedestal evolution", gainValues[gainId]), nev, 1. , xnev);
    hGainEB[gainId] = new TProfile(Form("GainEB%i", gainValues[gainId]), Form("EB gain%i pedestal evolution", gainValues[gainId]), nev, 1. , xnev);
    hGainEE[gainId] = new TProfile(Form("GainEE%i", gainValues[gainId]), Form("EE gain%i pedestal evolution", gainValues[gainId]), nev, 1. , xnev);
    hEBDiff[gainId] = new TH2F(Form("EBDiff_%i",gainId),
			       Form("EB pedestal difference gain %i",gainValues[gainId]),360,1.,361.,171,-85.,86.);
    hEEDiff[gainId] = new TH2F*[2];
    for (int Zside = 0; Zside <2; Zside++) {
      int izz = Zside;
      if(Zside == 0) izz = -1;
      hEEDiff[gainId][Zside] = new TH2F(Form("EEDiff_%i_%i", gainId, Zside),
					Form("Endcap pedestal difference gain %i side %i", gainValues[gainId], izz ),100, 1.,101., 100, 1., 101);
    }
    begtime[gainId] = Timestamp(std::numeric_limits<TimeValue_t>::max());
    endtime[gainId] = Timestamp(1);
  }

  for(int i = 0; i < 10; i++)
    nevent[i] = 0;

  vCovMatEB = new vector<H4CovMat*>(kEBChannels);
  vCovMatEE = new vector<H4CovMat*>(kEEChannels);
  if(SampleCorr_) {
    for (unsigned int i = 0; i < vCovMatEB->size(); i++) {
      (*vCovMatEB)[i] = new H4CovMat();
      (*vCovMatEB)[i]->init(runnumber_);
      (*vCovMatEB)[i]->SetMode(true);
    }
    for (unsigned int i = 0; i < vCovMatEE->size(); i++) {
      (*vCovMatEE)[i] = new H4CovMat();
      (*vCovMatEE)[i]->init(runnumber_);
      (*vCovMatEE)[i]->SetMode(true);
    }
  }

  cout << "Exiting beginJob" << endl;
}//beginJob

//========================================================================
void
EcalPedestal::endJob() {
//========================================================================

  cout << "Entering endJob" << endl;

  TFile f(Form("ecalPedestal_%d.root",runnumber_),"RECREATE");
  float noiseup[kGains] = {5., 2., 1.};
  //  float noisecut[kGains] = {1.5, 1.0, 0.7};
  //  float noisecut[kGains] = {1.8, 1.1, 0.8};
  //  float noisecut[kGains] = {1.9, 1.2, 0.9};       // 161017
  //  float noisecut[kGains] = {2.1, 1.3, 0.9};       // 170821
  //  float noisecut[kGains] = {2.1, 1.4, 0.9};       // 180716
  float noisecut[kGains] = {3.0, 1.5, 0.9};         // 220615
  //  float noisecutEE[kGains] = {2.5, 1.5, 1.0};
  float noisecutEE[kGains] = {3.0, 1.5, 0.9};       // 220615
  float hitup[kGains] = {5., 5., 10.};
  float sum25up[kGains] = {15., 10., 5.};
  float sum9up[kGains] = {10., 5., 2.5};

  cout << " Nb EB events         " << nevent[0] << endl
       << " Nb EB events gain 12 " << nevent[1] << endl
       << " Nb EB events gain 6  " << nevent[2] << endl
       << " Nb EB events gain 1  " << nevent[3] << endl
       << " Nb EE events         " << nevent[5] << endl
       << " Nb EE events gain 12 " << nevent[6] << endl
       << " Nb EE events gain 6  " << nevent[7] << endl
       << " Nb EE events gain 1  " << nevent[8] << endl;

  char buf[256];
  cout << "\n\n  ***************   Time    ****************" << endl;
  /*   does not work,   why?
  for(int ig = 0; ig < 3; ig++) {
    unsigned long long bbb = begtime[ig].value(), eee = endtime[ig].value();
    unsigned long ddd = eee - bbb;
    unsigned long mic = static_cast<unsigned int>(0xFFFFFFFF & ddd);
    unsigned long sec = static_cast<unsigned int>(ddd >> 32);
    cout << "gain " << ig << " beg " << bbb << " end " << eee << " duration " << ddd 
	 << " " << sec << "." << mic << endl;
  }
  */
  double micro, btOff[3], etOff[3];
  // always 6 digits in microsecondOffset!
  for(int ig = 0; ig < 3; ig++) {
    btOff[ig] = begtime[ig].microsecondOffset();
    etOff[ig] = endtime[ig].microsecondOffset();
    for(int di = 0; di < 5; di ++) {
      if(btOff[ig] > 100000) break;
      else {
	cout << " gain " << ig << " begin microsecondOffset was " << btOff[ig];
	btOff[ig] =btOff[ig] * 10;
	cout << " and corrected to " << btOff[ig] << endl;
      }
    }
    for(int di = 0; di < 5; di ++) {
      if(etOff[ig] > 100000) break;
      else {
	cout << " gain " << ig << " end microsecondOffset was " << etOff[ig];
	etOff[ig] = etOff[ig] * 10;
 	cout << " and corrected to " << etOff[ig] << endl;
      }
    }
  }
  time_t unixb = begtime[2].unixTime();
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixb));
  cout << " gain 12 begins at " << buf << "." << begtime[2].microsecondOffset();
  time_t unixe = endtime[2].unixTime();
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixe));
  if(etOff[2] > btOff[2])
    micro = double(etOff[2] - btOff[2])/1000000.;
  else 
    micro = -double(btOff[2] - etOff[2])/1000000.;
  cout << " and ends at " << buf << "." << etOff[2] 
       << " duration " << unixe - unixb + micro << endl;

  unixb = begtime[1].unixTime();
  if(etOff[2] > btOff[1])
    micro = -double(etOff[2] - btOff[1])/1000000.;
  else 
    micro = double(btOff[1] - etOff[2])/1000000.;
  cout << " time between gain 12 and 6 " << unixb - unixe + micro << endl;
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixb));

  cout << " gain 6 begins at " << buf << "." << btOff[1];
  unixe = endtime[1].unixTime() + etOff[1]/1000000.;;
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixe));
  if(etOff[1] > btOff[1])
    micro = double(etOff[1] - btOff[1])/1000000.;
  else 
    micro = -double(btOff[1] - etOff[1])/1000000.;
  cout << " and ends at " << buf << "." << etOff[1]
       << " duration " << unixe - unixb + micro << endl;

  unixb = begtime[0].unixTime();
  if(btOff[0] > etOff[1])
    micro = double(btOff[0] - etOff[1])/1000000.;
  else 
    micro = -double(etOff[1] - btOff[0])/1000000.;
  cout << " time between gain 6 and 1 " << unixb - unixe + micro << endl;

  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixb));
  cout << " gain 1 begins at " << buf << "." << btOff[0];
  unixe = endtime[0].unixTime();
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unixe));
  if(etOff[0] > btOff[0]) {
    micro = double(etOff[0] - btOff[0])/1000000.;
  }
  else {
    micro = -double(btOff[0] - etOff[0])/1000000.;
  }
  cout << " and ends at " << buf << "." << etOff[0]
       << " duration " << unixe - unixb + micro << endl;
  cout << "  ***************************************\n\n" << endl;
  // Create output text file for standard pedestals
  std::ofstream fPedestal;
  fPedestal.open(Form("pedestal_%d",runnumber_));
  std::ofstream fDQM;
  fDQM.open(Form("DQMlike_%d",runnumber_));
  std::ofstream fDB;
  fDB.open(Form("DBinput_%d",runnumber_));

  if(nevent[0] > 0) {  // barrel data present
    hEntries->Write();
    hGains->Write();

    TH2F** hEBPed = new TH2F*[kGains];
    TH2F** hEBNoise = new TH2F*[kGains];
    TH2F** hEBNoiseLF = new TH2F*[kGains];
    TH2F** hEBNoiseHF = new TH2F*[kGains];

    TH2F** hEBNoiseHit = new TH2F*[kGains];
    TH2F** hEBNoiseMinMax = new TH2F*[kGains];

    for (int gainId = 0; gainId < kGains; gainId++) {
      hPedsample[gainId]->Write();
      if(hEBEntry[gainId]->GetEntries() > 0.) hEBEntry[gainId]->Write();
      hEBPed[gainId] = new TH2F(Form("EBPed_%d", gainId),
				Form("pedestal gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);
      hEBNoise[gainId] = new TH2F(Form("EBNoise_%d", gainId),
				  Form("Total noise gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);
      hEBNoiseLF[gainId] = new TH2F(Form("EBNoiseLF_%d", gainId),
				    Form("Low frequency noise gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);
      hEBNoiseHF[gainId] = new TH2F(Form("EBNoiseHF_%d", gainId),
				    Form("High frequency noise gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);

      hEBNoiseHit[gainId] = new TH2F(Form("EBNoiseHit_%d", gainId),
				     Form("Weights Hit noise gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);
      hEBNoiseMinMax[gainId] = new TH2F(Form("EBNoiseMax_%d", gainId),
					Form("Min MaxHit noise gain %i", gainValues[gainId]),360,1.,361.,171,-85.,86.);

    }   // loop on gains

    TH2F*** hTotalNoiseVsPos = new TH2F**[36];
    TH1F** hPedestal = new TH1F*[kGains];
    TH1F** hTotNoise = new TH1F*[kGains];
    TH1F** hLFNoise  = new TH1F*[kGains];
    TH1F** hHFNoise  = new TH1F*[kGains];

    TH1F** hPedestalHit = new TH1F*[kGains];
    TH1F** hTotNoiseHit = new TH1F*[kGains];

    TH1F** hPedestalMinMax = new TH1F*[kGains];
    TH1F** hTotNoiseMinMax = new TH1F*[kGains];

    for (int gainId = 0; gainId < kGains; gainId++) {
      std::string hTitle = Form("Pedestal gain x%d", gainValues[gainId]);
      hPedestal[gainId] = new TH1F(Form("Pedestal_%d", gainId),
				   hTitle.c_str(), 100, 150., 250.);
      hTitle = Form("Tot Noise gain x%d", gainValues[gainId]);
      hTotNoise[gainId] = new TH1F(Form("TotNoise_%d", gainId),
				   hTitle.c_str(), 100, 0., noiseup[gainId]);
      hTitle = Form("LF noise gain x%d", gainValues[gainId]);
      hLFNoise[gainId] = new TH1F(Form("LFNoise_%d",  gainId),
				  hTitle.c_str(), 100, 0., noiseup[gainId]);
      hTitle = Form("HF noise gain x%d", gainValues[gainId]);
      hHFNoise[gainId] = new TH1F(Form("HFNoise_%d", gainId),
				  hTitle.c_str(), 100, 0., noiseup[gainId]);

      hTitle = Form("Pedestal for weights hits gain x%d", gainValues[gainId]);
      hPedestalHit[gainId] = new TH1F(Form("PedestalHit_%d", gainId),
				      hTitle.c_str(),  100, -2., 2.);
      hTitle = Form("Tot Noise forweights hits gain x%d", gainValues[gainId]);
      hTotNoiseHit[gainId] = new TH1F(Form("TotNoiseHit_%d", gainId),
				      hTitle.c_str(), 100, 0., hitup[gainId]);
      hTitle = Form("Pedestal for Min Max hits gain x%d", gainValues[gainId]);
      hPedestalMinMax[gainId] = new TH1F(Form("PedestalMinMax_%d", gainId),
					 hTitle.c_str(),  100, -10., 10.);
      hTitle = Form("Tot Noise for  Min Max hits gain x%d", gainValues[gainId]);
      hTotNoiseMinMax[gainId] = new TH1F(Form("TotNoiseMinMax_%d", gainId),
					 hTitle.c_str(), 100, 0., 10.);
    }  // loop on gain

    for(int SM = 0; SM <36; SM++) {
      hTotalNoiseVsPos[SM] = new TH2F*[kGains];
    }
    TH1F** hNoise_25 = new TH1F*[kGains];
    TH1F** hNoise_9 = new TH1F*[kGains];

    TH1F** hNoise_25Xt = new TH1F*[kGains];
    TH1F** hNoise_9Xt = new TH1F*[kGains];
    TH1F** hTotCorr25 = new TH1F*[kGains];
    TH1F** hTotCorr9 = new TH1F*[kGains];
    TH1F** hTotCorrp25 = new TH1F*[kGains];
    TH1F** hTotCorrp9 = new TH1F*[kGains];

    TH1F** hPedestalSingle = new TH1F*[kGains];
    TH1F** hTotNoiseSingle = new TH1F*[kGains];
    TH1F** hNoiseSingle_9 = new TH1F*[kGains];
    TH1F** hNoiseSingle_25 = new TH1F*[kGains];

    for (int gainId = 0; gainId < kGains; gainId++) {
      for(int SM = 0; SM <36; SM++) {
	string hTitle = Form("Noise vs (#eta, #phi) gain x%d", gainValues[gainId]);
	hTotalNoiseVsPos[SM][gainId] = new TH2F(Form("TotalNoiseVsPos_SM%d_%d", SM + 1, gainId),
						hTitle.c_str(), 85, 0., 85, 20, 0., 20.);
      }  // end loop on SM
      string hTitle = Form("Noise sum 25 Xtals, gain x%d", gainValues[gainId]);
      hNoise_25[gainId] = new TH1F(Form("Noise_25_%d", gainId),  
				   hTitle.c_str(), 500, 0., sum25up[gainId]);
      hTitle = Form("Noise sum 9 Xtals, gain x%d", gainValues[gainId]);
      hNoise_9[gainId] = new TH1F(Form("Noise_9_%d", gainId),  
				  hTitle.c_str(), 500, 0., sum9up[gainId]);

      hTitle = Form("Noise sum ind. 25 Xtals, gain x%d", gainValues[gainId]);
      hNoise_25Xt[gainId] = new TH1F(Form("Noise_25Xt_%d", gainId),
					 hTitle.c_str(), 500, 0., sum25up[gainId]);
      hTitle = Form("Noise sum ind. 9 Xtals, gain x%d", gainValues[gainId]);
      hNoise_9Xt[gainId] = new TH1F(Form("Noise_9Xt_%d", gainId),
					hTitle.c_str(), 500, 0., sum9up[gainId]);
      hTitle = Form("Corr. Noise 25 xtals, gain x%d", gainValues[gainId]);
      hTotCorr25[gainId] = new TH1F(Form("TotCorr25_%d", gainId),
					hTitle.c_str(), 500, 0., 2.);
      hTitle = Form("Corr. Noise perc.25 xtals, gain x%d", gainValues[gainId]);
      hTotCorrp25[gainId] = new TH1F(Form("TotCorrp25_%d", gainId),
					 hTitle.c_str(), 250, 0., 1.);
      hTitle = Form("Corr. Noise 9 xtals, gain x%d", gainValues[gainId]);
      hTotCorr9[gainId] = new TH1F(Form("TotCorr9_%d", gainId),
				       hTitle.c_str(), 500, 0., 2.);
      hTitle = Form("Corr. Noise perc.9 xtals, gain x%d", gainValues[gainId]);
      hTotCorrp9[gainId] = new TH1F(Form("TotCorrp9_%d", gainId),
					hTitle.c_str(), 250, 0., 1.);

      hTitle = Form("Noise Single 9 gain x%d", gainValues[gainId]);
      hNoiseSingle_9[gainId] = new TH1F(Form("NoiseSingle_9_%d", gainId),
					    hTitle.c_str(), 100, 0., sum9up[gainId]);
      hTitle = Form("Noise Single 25 gain x%d", gainValues[gainId]);
      hNoiseSingle_25[gainId] = new TH1F(Form("NoiseSingle_25_%d", gainId),
					     hTitle.c_str(), 100, 0., sum25up[gainId]);
      hTitle = Form("PedestalSingle gain x%d", gainValues[gainId]);
      hPedestalSingle[gainId] = new TH1F(Form("PedestalSingle_%d", gainId),
					     hTitle.c_str(), 100, 150., 250.);
      hTitle = Form("Tot Noise Single gain x%d", gainValues[gainId]);
      hTotNoiseSingle[gainId] = new TH1F(Form("TotNoiseSingle_%d", gainId),
					     hTitle.c_str(), 100, 0., noiseup[gainId]);

    }  // loop on gain
    //    TH2F** hCorrEB = new TH2F*[kEBChannels];

    float Totmean = 0.;
    float Totsig = 0.;
    //    float LFmean = 0.;
    float LFsig = 0.;
    //    float HFmean = 0.;
    float HFsig = 0.;
    noise_1_.resize(kEBChannels*kGains, 0.);
    noise_25_.resize(kEBChannels*kGains, 0.);
    noisesum_25_.resize(kEBChannels*kGains, 0.);
    nXtals_25_.resize(kEBChannels*kGains, 0);
    noise_9_.resize(kEBChannels*kGains, 0.);
    noisesum_9_.resize(kEBChannels*kGains, 0.);
    nXtals_9_.resize(kEBChannels*kGains, 0);

    noiseSingle_25_.resize(kEBChannels*kGains, 0.);
    noiseSingle_9_.resize(kEBChannels*kGains, 0.);

    int Noentry[36][68][3];
    int NoentrySM[36][3];
    //    int NoentryHit[36][68][3];
    for (int gainId = 0; gainId < kGains; gainId++) {
      for (int SM = 0; SM <36; SM++) {
	NoentrySM[SM][gainId] = 0;
	for (int tower = 0; tower < 68; tower++) {
	  Noentry[SM][tower][gainId] = 0;
	  //	  NoentryHit[SM][tower][gainId] = 0;
	}
      }
    }

    Int_t maxGain, minGain;
    if(runType_ == "Pedes") {
      minGain = 0;
      maxGain = kGains;
    }
    else if(runType_ == "Ped_6") {
      minGain = 1;
      maxGain = 2;
    }
    else if(runType_ == "Other") {
      minGain = 0;
      maxGain = 1; // only gain 12
    }
    else if(runType_ == "Calib") {
      minGain = 0;
      maxGain = 1; // only gain 12
    }
    else {
      cout << " unknown run Type" << runType_ << endl;
      minGain = 0;
      maxGain = 1; // only gain 12
    }

    vector<int>::iterator result;
    int GoodEBChannels = 0, GoodEBChannelsEta[2][17];
    Double_t corrMean[10][10], corrMeanEta[2][17][10][10];
    Double_t corrRSM[10][10], corrRMSEta[2][17][10][10];
    for(int is1 = 0; is1 < 10; is1++) {
      for(int is2 = 0; is2 <= is1; is2++) {
	corrMean[is1][is2] = 0;
	corrRSM[is1][is2] = 0;
	for(int side = 0; side < 2; side++) {
	  for(int eta = 0; eta < 17; eta++) {
	    corrMeanEta[side][eta][is1][is2] = 0;
	    corrRMSEta[side][eta][is1][is2] = 0;
	  }
	}
      }
    }
    for(int side = 0; side < 2; side++)
      for(int eta = 0; eta < 17; eta++)
	GoodEBChannelsEta[side][eta] = 0;

    // Loop over Ecal barrel
    for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
      EBDetId myEBDetId = EBDetId::unhashIndex(iChannel);
      int ieta = myEBDetId.ieta();  // -85:-1,1:85
      int iphi = myEBDetId.iphi();  // 1:360
      int iChperSM = myEBDetId.ic(); // 1:1700
      int ietaSM = myEBDetId.ietaSM();  // 1:85
      int iphiSM = myEBDetId.iphiSM();  // 1:20
      int SM = myEBDetId.ism();  // 1:36
      int towerID = myEBDetId.tower().iTT(); // 1:68
      int SMm1 = SM - 1;  // 0:35 for arrays
      if(SM < 1 || SM > 36)  cout << "**** Aie ***  channel " << iChannel
				  << " SM " << SM << endl;
      int ietamin = ieta - 2;
      int ietamax = ieta + 2;
      if(ieta > 0) {  // take care there is no eta = 0 !!!
	ietamin = ieta - 3;
	ietamax = ieta + 1;
      }
      bool alreadyKnownDead = false;
      result = find(maskedChannels_.begin(), maskedChannels_.end(), iChannel);
      if (result != maskedChannels_.end()) alreadyKnownDead = true;
      // Loop over different gains
      for (int gainId = 0; gainId < 3; gainId++) { 
	if (result == maskedChannels_.end()) {    // do not use bad channels
	  if(meanEBped[gainId][iChannel][firstReadEvent_ - 1] == -999.) {  // first event does not have Ecal data :-(
	    if(iChannel == 0) cout << " no data for event " << firstReadEvent_ << " for gain " << gainId << endl;
	    meanEBped[gainId][iChannel][firstReadEvent_ - 1] = meanEBped[gainId][iChannel][firstReadEvent_];
	  }
	  for (int ievt = 0; ievt < 300; ievt++) { 
	    float diffped = meanEBped[gainId][iChannel][ievt] - meanEBped[gainId][iChannel][0];
	    if(meanEBped[gainId][iChannel][ievt] == -999.) {
	      //	      if(!alreadyKnownDead) cout << "EB event " << ievt << " channel " << iChannel << " gain " << gainValues[gainId] << " mean -999."  << endl;
	    }
	    else {
	      hGain[gainId]->Fill(ievt, diffped);
	      hGainEB[gainId]->Fill(ievt, diffped);
	      //	  if(iChannel%10000 == 0) cout << "EB event " << xEvent << " channel " << iChannel << " gain " << expectedGain
	      //				      << " 1st mean " << meanEBped[iChannel] << " diffped " << diffped << endl;
	      if(abs(diffped) > 10.) {
		//	    cout << "******** EB event " << xEvent << " gain " << expectedGain << " channel " << iChannel
		//		 << " diffped " << diffped << " phi " << iphi << " eta " << ieta << endl;
		hEBDiff[gainId]->Fill(iphi, ieta, 1.);
	      }
	    }  // bad pedestal
	  }   // loop over events
	}
	int arrayId = gainId*kEBChannels + iChannel;

	// Standard way...	
	if (nEntriesStandard_[arrayId]) {
	  double xmean = pedestalStandard_[arrayId] / double(nEntriesStandard_[arrayId]);
	  Totmean = xmean;
	  double sigmean = pedestal2Standard_[arrayId] / double(nEntriesStandard_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  //	  if(iChannel > 7814 && iChannel < 7820) 
	  //	    cout << " final mean " << xmean << " sum ped " << pedestalStandard_[arrayId]
	  //		 << " nb entries " << nEntriesStandard_[arrayId] << endl;

	  Totsig = sqrt(Totsig2);
	  float noisecutEta = noisecut[gainId];
	  if(gainId == 0 && ietaSM > 40) {                                           // 2018
	    if(ietaSM > 60) noisecutEta = noisecut[gainId] + 0.5;
	    else noisecutEta = noisecut[gainId] + 0.2;
	  }
	  //	  if (Totsig == 0. || Totsig > noisecut[gainId])
	  if (Totsig == 0. || Totsig > noisecutEta)
	    cout << " Channel " <<  iChannel
		 << " SM " << SM << " SMCal " << SMCal[SMm1] 
		 << " towerID " << towerID
		 << " eta " << ietaSM << " phi " << iphiSM
		 << " nb in SM " << iChperSM
		 << " gain x" <<  gainValues[gainId]
		 << " Tot noise " << Totsig
		 << endl;
	  if(!alreadyKnownDead && (Totmean < 170. || Totmean > 300.)) 
	    cout << "++++++++++++++++++++++++++\n"
		 << "+ new dead channel " << iChannel << " +\n"
		 << "+       ped " <<  Totmean << "  gain  " << gainValues[gainId] << "       +\n"
		 << "++++++++++++++++++++++++++\n" << endl;
	  if(!alreadyKnownDead && Totsig == 0)
	    cout << "++++++++++++++++++++++++++\n"
		 << "+ new dead channel " << iChannel << " +\n"
		 << "+          RMS 0         +\n"
		 << "++++++++++++++++++++++++++\n" << endl;
	  if(gainId == 0) {
	    fPedestal << std::setw(5) << iChannel << " "
		      << std::setw(3) << iphi << " "
		      << std::setw(3) << ieta << " ";
	    fDB << std::setw(5) << iChannel + 1;
	  }
	  fPedestal << std::setw(7) << std::setprecision(5) << xmean << " " 
		    << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  fDB << std::setw(7) << std::setprecision(5) << xmean << " " 
	      << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  hEBPed[gainId]->Fill(iphi,ieta,xmean);
	  //	float eta = ieta;
	  //	if(ieta > 0) eta = ieta - 1;
	  hEBNoise[gainId]->Fill(iphi,ieta,Totsig);
	  hPedestal[gainId]->Fill(Totmean);
	  hTotNoise[gainId]->Fill(Totsig);
	  hTotalNoiseVsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, Totsig);
	  // Compute sliding windows Sum9 and sum25
	  for (int ip = iphi-2; ip <= iphi+2; ++ip) {
	    if(ip > 0 && ip <= 360) {
	      for (int ie = ietamin; ie <= ietamax; ++ie) {
		if (ie >= -85  && ie < 85 ) {
		  int id =  gainId * kEBChannels + (85 + ie) * 360 + ip -1;
		  noise_25_[id] += Totsig2;
		  noisesum_25_[id] += Totsig;
		  nXtals_25_[id]++;
		  if(ip != iphi-2 && ip != iphi+2 && ie != ietamin && ie != ietamax) {
		    noise_9_[id] += Totsig2;
		    noisesum_9_[id] += Totsig;
		    nXtals_9_[id]++;
		  }
		}
	      }
	    }
	  } // End sliding windows computation
	  // Low Frequency noise
	  xmean = avgSamples_[arrayId] / double(nEntriesAvg_[arrayId]);
	  //	  LFmean = xmean;
	  sigmean = avgSample2_[arrayId] / double(nEntriesAvg_[arrayId]);
	  Totsig2 = sigmean - xmean * xmean;
	  if(Totsig2 >= 0.)
	    LFsig = sqrt(sigmean - xmean * xmean);
	  else {
	    cout << "EB channel " << iChannel << " pb with LF variance " << Totsig2 << " sigmean " << sigmean
		 << " xmean " << xmean * xmean << endl;
	    LFsig = 0.0;
	  }
	  //	hLFNoise[SMm1][gainId]->Fill(LFsig);
	  //	hLFNoiseSM[SMm1][gainId]->Fill(iChperSM, LFsig);
	  hLFNoise[gainId]->Fill(LFsig);
	  fPedestal << std::setw(7) << std::setprecision(4) << LFsig << " ";
	  hEBNoiseLF[gainId]->Fill(iphi,ieta,LFsig);
	  // High Frequency noise
	  xmean = pedavg_[arrayId]/double(nEntriesStandard_[arrayId]);
	  //	  HFmean = xmean;
	  sigmean =  pedavg2_[arrayId]/double(nEntriesStandard_[arrayId]);
	  HFsig = sqrt(sigmean - xmean * xmean);
	  hHFNoise[gainId]->Fill(HFsig);
	  fPedestal << std::setw(7) << std::setprecision(4) << HFsig;
	  hEBNoiseHF[gainId]->Fill(iphi,ieta,HFsig);
	}  // entries > 0
	else {
	  Noentry[SMm1][towerID - 1][gainId]++;
	  NoentrySM[SMm1][gainId]++;
	  double xmean = -999.;
	  double Totsig = -999.;
	  double LFsig = -999.;
	  double HFsig = -999.;
	  if(gainId == 0) {
	    fPedestal << std::setw(5) << iChannel << " "
		      << std::setw(3) << iphi << " "
		      << std::setw(3) << ieta << " ";
	    fDB << std::setw(5) << iChannel + 1;
	  }
	  fPedestal << std::setw(7) << std::setprecision(5) << xmean << " " 
		    << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  fPedestal << std::setw(7) << std::setprecision(4) << LFsig << " ";
	  fPedestal << std::setw(7) << std::setprecision(4) << HFsig;
	  fDB << std::setw(7) << std::setprecision(5) << xmean << " " 
	      << std::setw(7) << std::setprecision(4) << Totsig << " ";
	}  // nEntriesStandard

	//  Now uncalib hits
	if (nEntriesHit_[arrayId]) {
	  double xmean = pedestalHit_[arrayId] / double(nEntriesHit_[arrayId]);
	  Totmean = xmean;
	  double sigmean = pedestal2Hit_[arrayId] / double(nEntriesHit_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  Totsig = sqrt(Totsig2);
	  noise_1_[arrayId] = Totsig;
	  hEBNoiseHit[gainId]->Fill(iphi,ieta,Totsig);
	  hPedestalHit[gainId]->Fill(Totmean);
	  hTotNoiseHit[gainId]->Fill(Totsig);
	}  // nEntries not 0
	//  Now MinMax hits
	if (nEntriesMinMax_[arrayId]) {
	  double xmean = pedestalMinMax_[arrayId] / double(nEntriesMinMax_[arrayId]);
	  Totmean = xmean;
	  double sigmean = pedestal2MinMax_[arrayId] / double(nEntriesMinMax_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  Totsig = sqrt(Totsig2);
	  hEBNoiseMinMax[gainId]->Fill(iphi,ieta,Totsig);
	  hPedestalMinMax[gainId]->Fill(Totmean);
	  hTotNoiseMinMax[gainId]->Fill(Totsig);
	}  // nEntries not 0

      }  // End loop over gains
      fPedestal << std::endl;
      fDB << std::endl;
      if (!alreadyKnownDead && nEntriesStandard_[iChannel]) {
	GoodEBChannels++;
	if(SampleCorr_) (*vCovMatEB)[iChannel]->ComputeEstimators();
	int eta = (myEBDetId.ietaSM() - 1)/5;    // from 1-85 to 0-16
	if(eta < 0 || eta > 16) {
	  cout << " problem with channel " << iChannel << " eta = " << eta << endl;
	  exit(-1);
	}
	int side = 0;
	if(ieta > 0) side = 1;
	GoodEBChannelsEta[side][eta]++;
	if(SampleCorr_) {
	  for(int is1 = 0; is1 < 10; is1++) {
	    for(int is2 = 0; is2 <= is1; is2++) {
	      float corr = (*vCovMatEB)[iChannel]->getCorrelation(is1, is2);
	      corrMean[is1][is2] += corr;
	      corrRSM[is1][is2] += corr * corr;
	      corrMeanEta[side][eta][is1][is2] += corr;
	      corrRMSEta[side][eta][is1][is2] += corr * corr;
	      //	    if(iChannel%10000 == 0) {
	      //	      hCorrEB[iChannel]->Fill(is1, is2, corr);
	      //	      if(is2 !=is1) hCorrEB[iChannel]->Fill(is2, is1, corr);
	      //	    }
	    }
	    //	cout << endl;
	  }  // end loop over samples
	  //	     if(iChannel%10000 == 0) hCorrEB[iChannel]->Write();
	}   // study sample correlation?
      } // good channel
      for (int gainId = 2; gainId > -1; gainId--) { 
	int arrayId = gainId*kEBChannels + iChannel;
	if (nEntriesStandard_[arrayId]) {
	  if(gainId == 2) {
	    if(SM <10) fDQM << "55555  10110" << SM;
	    else fDQM << "55555  1011" << SM;
	    if(iChperSM < 10) fDQM << "000" << iChperSM;
	    else if(iChperSM < 100) fDQM << "00" << iChperSM;
	    else if(iChperSM < 1000) fDQM << "0" << iChperSM;
	    else fDQM  << iChperSM;
	  }
	  double xmean = pedestalStandard_[arrayId] / double(nEntriesStandard_[arrayId]);
	  double sigmean = pedestal2Standard_[arrayId] / double(nEntriesStandard_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  Totsig = sqrt(Totsig2);
	  fDQM << std::setw(7) << std::setprecision(5) << xmean << " " 
	       << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  if(gainId == 0) fDQM << " 0" << std::endl;
	}  // only if read out
      }   // End loop over gains
    }    // End loop over channels
    // print missing SM and towers
    /*
    for (int gainId = minGain; gainId < maxGain; gainId++) { 
      
      for (int SM = 0; SM <36; SM++) {
	if(NoentrySM[SM][gainId] == 1700)
	  cout << " No entry for SM " << SM + 1 << " SMCal " << SMCal[SM] 
	       << " gain " <<  gainValues[gainId]
	       << endl;
	else {
	  for (int tower = 0; tower < 68; tower++) {
	    if(Noentry[SM][tower][gainId] > 0)
	      cout << " No entry for " << Noentry[SM][tower][gainId] 
		   << " channels in tower " << tower + 1
		   << " SM " << SM + 1 << " SMCal " << SMCal[SM] 
		   << " gain " <<  gainValues[gainId]
		   << endl;
	    
	  }
	}
      }
    }  // loop over gains
    */

    // check sliding windows computation
    // Loop over Ecal barrel channels
    for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
      //      EBDetId myEBDetId = EBDetId::unhashIndex(iChannel);
      //    int ietaSM = myEBDetId.ietaSM();
      //    int iphiSM = myEBDetId.iphiSM();
      //    int SM = myEBDetId.ism();  // 1:36
      //    int SMm1 = SM - 1;
      for (int gainId = 0; gainId < kGains; gainId++) {
	int arrayId = gainId*kEBChannels + iChannel;
	//  hits
	if(nXtals_9_[arrayId] > 0 && nXtals_25_[arrayId] > 0) {
	  float noise9 = sqrt(noise_9_[arrayId] / nXtals_9_[arrayId]);
	  float noise25 = sqrt(noise_25_[arrayId] / nXtals_25_[arrayId]);
	  float noisesum25 = noisesum_25_[arrayId];
	  float noisesum9 = noisesum_9_[arrayId];
	  //      float noise25ratio = 1 - noisesum25/noise25;
	  //      float noise9ratio = 1 - noisesum9/noise9;
	  float corr25 = sqrt(noise25*noise25-noisesum25*noisesum25);
	  float corrp25 = corr25/noise25;
	  float corr9 = sqrt(noise9*noise9-noisesum9*noisesum9);
	  float corrp9 = corr9/noise9;
	  //      hTotNoise25[SMm1][gainId]->Fill(noise25ratio);
	  //      hTotNoise9[SMm1][gainId]->Fill(noise9ratio);
	  //      hTotNoise25VsPos[SMm1][gainId]->SetBinContent(ieta + 1, iphi + 1, noise25ratio);
	  //      hTotNoise9VsPos[SMm1][gainId]->SetBinContent(ieta + 1, iphi + 1, noise9ratio);
	  //	hNoise_25VsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, noise25);
	  //	hNoise_9VsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, noise9);
	  hNoise_25Xt[gainId]->Fill(noisesum25);
	  hNoise_9Xt[gainId]->Fill(noisesum9);
	  //	hNoise_25XtVsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, noisesum25);
	//	hNoise_9XtVsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, noisesum9);
	  hTotCorr25[gainId]->Fill(corr25);
	  hTotCorrp25[gainId]->Fill(corrp25);
	  hTotCorr9[gainId]->Fill(corr9);
	  hTotCorrp9[gainId]->Fill(corrp9);
	//	hTotCorr25VsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, corr25);
	//	hTotCorr9VsPos[SMm1][gainId]->SetBinContent(ietaSM, iphiSM, corr9);
	}
	// hits

	// Sliding windows a la Alex
	//      if(ieta > 0 && ieta < 84 && iphi > 0  && iphi < 19) {
	//      if(ieta > -85  && ieta < 85 && iphi > 1  && iphi < 360) {
      }  // End loop over gains
      //    fPedestal << endl;
    }  // End loop over channels

    for (int gainId = minGain; gainId < maxGain; gainId++) { 
      hEBPed[gainId]->Write();
      hEBNoise[gainId]->Write();
      hEBNoiseLF[gainId]->Write();
      hEBNoiseHF[gainId]->Write();
      hEBNoiseHit[gainId]->Write();
      hEBNoiseMinMax[gainId]->Write();

      hPedestal[gainId]->Write();
      hTotNoise[gainId]->Write();
      hLFNoise[gainId]->Write();
      hHFNoise[gainId]->Write();

      // hits
      hPedestalHit[gainId]->Write();
      hTotNoiseHit[gainId]->Write();
      hPedestalMinMax[gainId]->Write();
      hTotNoiseMinMax[gainId]->Write();

      hNoise_25[gainId]->Write();
      hNoise_9[gainId]->Write();
      //    hNoise_25VsPos[gainId]->Write();
      //    hNoise_9VsPos[gainId]->Write();
      hNoise_25Xt[gainId]->Write();
      hNoise_9Xt[gainId]->Write();
      //    hNoise_25XtVsPos[gainId]->Write();
      //    hNoise_9XtVsPos[gainId]->Write();
      hTotCorr25[gainId]->Write();
      hTotCorrp25[gainId]->Write();
      hTotCorr9[gainId]->Write();
      hTotCorrp9[gainId]->Write();
      //    hTotCorr25VsPos[gainId]->Write();
      //    hTotCorr9VsPos[gainId]->Write();

      for(int SM = 1; SM < 37; SM++) {
	int SMm1 = SM - 1;
	hTotalNoiseVsPos[SMm1][gainId]->Write();
      }

      hNoise_25[gainId]->Write();
      hNoise_9[gainId]->Write();
      hNoise_25Xt[gainId]->Write();
      hNoise_9Xt[gainId]->Write();
      hTotCorr25[gainId]->Write();
      hTotCorrp25[gainId]->Write();
      hTotCorr9[gainId]->Write();
      hTotCorrp9[gainId]->Write();
      
    }
    if(SampleCorr_) {
      TH2F* hCorrMeanEB = new TH2F("CorrMeanEB","mean Sample correlation",
				   10 , 0., 10. , 10, 0., 10.);
      TH2F* hCorrRMSEB = new TH2F("CorrRMSEB","RMS Sample correlation",
				  10 , 0., 10. , 10, 0., 10.);
      TH2F*** hCorrMeanEta = new TH2F**[2];
      TH2F*** hCorrRMSEta = new TH2F**[2];
      for(int side = 0; side < 2; side++) {
	hCorrMeanEta[side] = new TH2F*[17];
	hCorrRMSEta[side] = new TH2F*[17];
	int zs = 1;
	if(side == 0) zs = -1;
	for(int eta = 0; eta < 17; eta++) {
	  hCorrMeanEta[side][eta] = new TH2F(Form("CorrMeanEta_%i_%i",side,eta),
					     Form("mean Sample correlation side %i ring %i",zs, eta),
					     10 , 0., 10. , 10, 0., 10.);
	  hCorrRMSEta [side][eta]= new TH2F(Form("CorrRMSEta_%i_%i",side,eta),
					    Form("RMS Sample correlation side %i ring %i",zs, eta),
					    10 , 0., 10. , 10, 0., 10.);
	}
      }
      for(int is1 = 0; is1 < 10; is1++) {
	for(int is2 = 0; is2 <= is1; is2++) {
	  double Mean = corrMean[is1][is2] / double(GoodEBChannels);
	  hCorrMeanEB->Fill(is1, is2, Mean);
	  if(is2 !=is1) hCorrMeanEB->Fill(is2, is1, Mean);
	  double RMS = corrRSM[is1][is2] / double(GoodEBChannels);
	  RMS = RMS - Mean * Mean;
	  RMS = sqrt(RMS);
	  hCorrRMSEB->Fill(is1, is2, RMS);
	  if(is2 !=is1) hCorrRMSEB->Fill(is2, is1, RMS);

	  for(int side = 0; side < 2; side++) {
	    for(int eta = 0; eta < 17; eta++) {
	      double Mean = corrMeanEta[side][eta][is1][is2] / double(GoodEBChannelsEta[side][eta]);
	      hCorrMeanEta[side][eta]->Fill(is1, is2, Mean);
	      if(is2 !=is1) hCorrMeanEta[side][eta]->Fill(is2, is1, Mean);
	      double RMS = corrRMSEta[side][eta][is1][is2] / double(GoodEBChannelsEta[side][eta]);
	      RMS = RMS - Mean * Mean;
	      RMS = sqrt(RMS);
	      hCorrRMSEta[side][eta]->Fill(is1, is2, RMS);
	      if(is2 !=is1) hCorrRMSEta[side][eta]->Fill(is2, is1, RMS);
	    }
	  }
	}
      }
      hCorrMeanEB->Write();
      hCorrRMSEB->Write();
      for(int side = 0; side < 2; side++) {
	for(int eta = 0; eta < 17; eta++) {
	  hCorrMeanEta[side][eta]->Write();
	  hCorrRMSEta [side][eta]->Write();
	}
      }
    }  // barrel data present
  }   // sample correlation study

  if(nevent[5] > 0) {  // End cap data present
    TH1F** hPedestalEE = new TH1F*[kGains];
    TH1F** hTotNoiseEE = new TH1F*[kGains];
    TH1F** hLFNoiseEE  = new TH1F*[kGains];
    TH1F** hHFNoiseEE  = new TH1F*[kGains];

    TH1F** hPedestalHitEE = new TH1F*[kGains];
    TH1F** hTotNoiseHitEE = new TH1F*[kGains];

    TH2F*** hEENoise   = new TH2F**[kGains];
    TH2F*** hEENoiseLF = new TH2F**[kGains];
    TH2F*** hEENoiseHF = new TH2F**[kGains];

    TH2F** hCorrEE = new TH2F*[kEEChannels];

    for (int gainId = 0; gainId < 3; gainId++) {
      std::string hTitle = Form("Endcap Pedestal gain x%d", gainValues[gainId]);
      hPedestalEE[gainId] = new TH1F(Form("PedestalEE_%d", gainId),
				   hTitle.c_str(), 100, 150., 250.);
      hTitle = Form("Endcap Tot Noise gain x%d", gainValues[gainId]);
      hTotNoiseEE[gainId] = new TH1F(Form("TotNoiseEE_%d", gainId),
				   hTitle.c_str(), 100, 0., noiseup[gainId]);
      hTitle = Form("Endcap LF noise gain x%d", gainValues[gainId]);
      hLFNoiseEE[gainId] = new TH1F(Form("LFNoiseEE_%d",  gainId),
				  hTitle.c_str(), 100, 0., noiseup[gainId]);
      hTitle = Form("Endcap HF noise gain x%d", gainValues[gainId]);
      hHFNoiseEE[gainId] = new TH1F(Form("HFNoiseEE_%d", gainId),
				  hTitle.c_str(), 100, 0., noiseup[gainId]);

      hTitle = Form("Endcap Pedestal for weights hits gain x%d", gainValues[gainId]);
      hPedestalHitEE[gainId] = new TH1F(Form("PedestalHitEE_%d", gainId),
				      hTitle.c_str(),  100, -2., 2.);
      hTitle = Form("Endcap Tot Noise forweights hits gain x%d", gainValues[gainId]);
      hTotNoiseHitEE[gainId] = new TH1F(Form("TotNoiseHitEE_%d", gainId),
				      hTitle.c_str(), 100, 0., hitup[gainId]);

      hEENoise[gainId] = new TH2F*[2];
      hEENoiseLF[gainId] = new TH2F*[2];
      hEENoiseHF[gainId] = new TH2F*[2];
      for (int Zside = 0; Zside <2; Zside++) {
	int izz = Zside;
	if(Zside == 0) izz = -1;
	hEENoise[gainId][Zside] = new TH2F(Form("EENoise_%i_%i", gainId, Zside),
				  Form("Endcap Total noise gain %i side %i", gainValues[gainId], izz ),100, 1.,101., 100, 1., 101);
	hEENoiseLF[gainId][Zside] = new TH2F(Form("EENoiseLF_%i_%i", gainId, Zside),
				  Form("Endcap Low Frequency noise gain %i side %i", gainValues[gainId], izz ),100, 1., 101.,100, 1., 101);
	hEENoiseHF[gainId][Zside] = new TH2F(Form("EENoiseHF_%i_%i", gainId, Zside),
				  Form("Endcap High Frequency noise gain %i side %i", gainValues[gainId], izz ),100, 1., 101.,100, 1., 101);
      }

    }

    cout << " kEEChannels " << kEEChannels << endl;
    Int_t NoentryEE[2][316][3];
    int NoentryDCCEE[2][9][3];
    int NoentryChannelEE[2][316][3][25];
    for (int gainId = 0; gainId < kGains; gainId++)
      for (int Zside = 0; Zside <2; Zside++) {
	for (int iSC = 0; iSC < 316; iSC++)
	  NoentryEE[Zside][iSC][gainId] = 0;
	for (int DCC = 0; DCC <9; DCC++)
	  NoentryDCCEE[Zside][DCC][gainId] = 0;
      }

    int GoodEEChannels = 0;
    Double_t corrMean[10][10];
    Double_t corrRSM[10][10];
    for(int is1 = 0; is1 < 10; is1++) {
      for(int is2 = 0; is2 <= is1; is2++) {
	corrMean[is1][is2] = 0;
	corrRSM[is1][is2] = 0;
      }
    }

    vector<int>::iterator result;
    for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
      EEDetId myEEDetId = EEDetId::unhashIndex(iChannel);
      int iz = myEEDetId.zside();
      int ix = myEEDetId.ix();
      int iy = myEEDetId.iy();
      int iSC = myEEDetId.isc();
      int izz = iz;
      if(iz == -1) izz = 0;
      //      EcalElectronicsId elecId = ecalElectronicsMap_.getElectronicsId(myEEDetId);
      int DCCid = elecId[iChannel].dccId() - 1;  // 0-8
      cout << " EE channel " << iChannel << " id " << DCCid << endl;
      if(DCCid > 44) {
	if(iz == -1) cout << " Endcap strange EE+ for DCCId " << DCCid + 1 << endl;
	DCCid -= 45;
      }
      else if(iz == 1) cout << " Endcap strange EE- for DCCId " << DCCid + 1 << endl;
      result = find(maskedEEChannels_.begin(), maskedEEChannels_.end(), iChannel);
      for (int gainId = 0; gainId < 3; gainId++) { 
	int arrayId = gainId*kEEChannels + iChannel;
	if (result == maskedEEChannels_.end()) {    // do not use bad channels
	  if(meanEEped[gainId][iChannel][firstReadEvent_ - 1] == -999.) {  // first event does not have Ecal data :-(
	    if(iChannel == 0) cout << " Endcap : no event "  << firstReadEvent_<< " for gain " << gainId << endl;
	    meanEEped[gainId][iChannel][firstReadEvent_ - 1] = meanEEped[gainId][iChannel][firstReadEvent_];
	  }
	  for (int ievt = 0; ievt < 300; ievt++) { 
	    float diffped = meanEEped[gainId][iChannel][ievt] - meanEEped[gainId][iChannel][0];
	    if(meanEEped[gainId][iChannel][ievt] == -999.) {
	      // cout << "event " << ievt << " EE channel " << iChannel << " gain " << gainValues[gainId] << " mean -999."  << endl;
	    }
	    else {
	      hGain[gainId]->Fill(ievt, diffped);
	      hGainEE[gainId]->Fill(ievt, diffped);
	      //	  if(iChannel%10000 == 0) cout << "EB event " << xEvent << " channel " << iChannel << " gain " << expectedGain
	      //				      << " 1st mean " << meanEBped[iChannel] << " diffped " << diffped << endl;
	      if(abs(diffped) > 10.) {
		//	    cout << "******** EB event " << xEvent << " gain " << expectedGain << " channel " << iChannel
		//		 << " diffped " << diffped << " phi " << iphi << " eta " << ieta << endl;
		hEEDiff[gainId][izz]->Fill(ix, iy, 1.);
	      }
	    }  // bad pedestal
	  }   // loop over events
	}

	// Standard way...	
	if (nEntriesStandardEE_[arrayId]) {
	  NoentryDCCEE[izz][DCCid][gainId]++;
	  double xmean = pedestalStandardEE_[arrayId] / double(nEntriesStandardEE_[arrayId]);
	  double Totmean = xmean;
	  double sigmean = pedestal2StandardEE_[arrayId] / double(nEntriesStandardEE_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  double Totsig = sqrt(Totsig2);
	  if (Totsig == 0. || Totsig > noisecutEE[gainId])
	    cout << " Channel " <<  iChannel
		 << " gain x" <<  gainValues[gainId]
		 << " Tot noise " << Totsig
		 << endl;
	  if(gainId == 0) {
	    fPedestal << std::setw(5) << iChannel << " "
		      << std::setw(3) << iz << " "
		      << std::setw(3) << ix << " "
		      << std::setw(3) << iy << " ";
	    fDB << std::setw(5) << iChannel + 61201;
	  }
	  fPedestal << std::setw(7) << std::setprecision(5) << xmean << " " 
		    << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  fDB << std::setw(7) << std::setprecision(5) << xmean << " " 
	      << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  hPedestalEE[gainId]->Fill(Totmean);
	  hTotNoiseEE[gainId]->Fill(Totsig);
	  hEENoise[gainId][izz]->Fill(ix,iy,Totsig);
	  // Low Frequency noise
	  xmean = avgSamplesEE_[arrayId] / double(nEntriesAvgEE_[arrayId]);
	  sigmean = avgSample2EE_[arrayId] / double(nEntriesAvgEE_[arrayId]);
	  sigmean = sqrt(sigmean - xmean * xmean);
	  double LFsig = sigmean;
	  hLFNoiseEE[gainId]->Fill(LFsig);
	  hEENoiseLF[gainId][izz]->Fill(ix,iy,LFsig);
	  fPedestal << std::setw(7) << std::setprecision(4) << LFsig << " ";
	  // High Frequency noise
	  xmean = pedavgEE_[arrayId]/double(nEntriesStandardEE_[arrayId]);
	  sigmean =  pedavg2EE_[arrayId]/double(nEntriesStandardEE_[arrayId]);
	  double HFsig = sqrt(sigmean - xmean * xmean);
	  hHFNoiseEE[gainId]->Fill(HFsig);
	  hEENoiseHF[gainId][izz]->Fill(ix,iy,HFsig);
	  fPedestal << std::setw(7) << std::setprecision(4) << HFsig;
	}  // nEntriesStandard
	else {
	  int index = NoentryEE[izz][iSC - 1][gainId];
	  NoentryChannelEE[izz][iSC - 1][gainId][index] = iChannel;
	  if(iSC > 316)
	    cout << " EE Channel " <<  iChannel << " Super Crystal " << iSC << endl;
	  NoentryEE[izz][iSC - 1][gainId]++;   // iSC runs from 1 to 316
	  if(gainId == 0) {
	    fPedestal << std::setw(5) << iChannel << " "
		      << std::setw(3) << iz << " "
		      << std::setw(3) << ix << " "
		      << std::setw(3) << iy << " ";
	    fDB << std::setw(5) << iChannel + 61201;
	  }
	  fPedestal << "   -999    -999    -999    -999";
	  fDB << "   -999    -999 ";
	}

	//  Now uncalib hits
	if (nEntriesHitEE_[arrayId]) {
	  double xmean = pedestalHitEE_[arrayId] / double(nEntriesHitEE_[arrayId]);
	  double Totmean = xmean;
	  double sigmean = pedestal2HitEE_[arrayId] / double(nEntriesHitEE_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  double Totsig = sqrt(Totsig2);
	  hPedestalHitEE[gainId]->Fill(Totmean);
	  hTotNoiseHitEE[gainId]->Fill(Totsig);
	}  // nEntriesHit not 0
      }  // loop over gain
      fPedestal << std::endl;
      fDB << std::endl;
      for (int gainId = 2; gainId > -1; gainId--) { 
	int arrayId = gainId*kEEChannels + iChannel;
	if (nEntriesStandardEE_[arrayId]) {
	  if(gainId == 2) {
	    if(iz == -1) fDQM << "55555  2010";
	    else fDQM << "55555  2012";
	    if(ix < 10) fDQM << "00" << ix;
	    else if(ix < 100) fDQM << "0" << ix;
	    else fDQM  << ix;
	    if(iy < 10) fDQM << "00" << iy;
	    else if(iy < 100) fDQM << "0" << iy;
	    else fDQM  << iy;
	  }
	  double xmean = pedestalStandardEE_[arrayId] / double(nEntriesStandardEE_[arrayId]);
	  double sigmean = pedestal2StandardEE_[arrayId] / double(nEntriesStandardEE_[arrayId]);
	  double Totsig2 = sigmean - xmean * xmean;
	  double Totsig = sqrt(Totsig2);
	  fDQM << std::setw(7) << std::setprecision(5) << xmean << " " 
	       << std::setw(7) << std::setprecision(4) << Totsig << " ";
	  if(gainId == 0) fDQM << " 0" << std::endl;
	  GoodEEChannels++;
	  if(SampleCorr_) {
	    (*vCovMatEE)[iChannel]->ComputeEstimators();
	    if(iChannel%2000 == 0) {
	      hCorrEE[iChannel] = new TH2F(Form("CorrEE_%d", iChannel),
					   Form("EE Time correlation xl %d", iChannel),
					   10 , 0., 10. , 10, 0., 10.);
	    } // choosen channel
	    for(int is1 = 0; is1 < 10; is1++) {
	      for(int is2 = 0; is2 <= is1; is2++) {
		float corr = (*vCovMatEE)[iChannel]->getCorrelation(is1, is2);
		corrMean[is1][is2] += corr;
		corrRSM[is1][is2] += corr * corr;
		if(iChannel%2000 == 0) {
		  hCorrEE[iChannel]->Fill(is1, is2, corr);
		  if(is2 !=is1) hCorrEE[iChannel]->Fill(is2, is1, corr);
		}
	      }
	    }  // end loop over samples
	    if(iChannel%2000 == 0) hCorrEE[iChannel]->Write();
	  }    // study sample correlation ?
	} // only if read out
      }  // loop over gain
    }   // loop over channels

    // print missing channel, SC, slice
    for (int gainId = 0; gainId < 3; gainId++) { 
      cout << " Results for gain " << gainValues[gainId] << endl;
      for (int Zside = 0; Zside <2; Zside++) {
	string  side = "-";
	int Fedid = 601;
	if(Zside == 1) {
	  side = "+";
	  Fedid = 646;
	}
	for (int DCC = 0; DCC < 9; DCC++) {
	  cout << " DCC " << DCC << endl; 
	  if(NoentryDCCEE[Zside][DCC][gainId] == 0) {
	    string sliceName = fedMap_->getSliceFromFed(Fedid + DCC);
  	    cout << " No entry for " << sliceName << endl;
	  }
	  cout << " noentry " << NoentryDCCEE[Zside][DCC][gainId] << endl;
	}
	for (int iSC = 0; iSC < 316; iSC++) {
	  cout << " sc " << iSC << endl;
	  cout << " sc " << iSC << " noentry " << NoentryEE[Zside][iSC][gainId] << endl;
	  // check if it is not part of am empty slice (already printed)
	  if(NoentryEE[Zside][iSC][gainId] > 0) {
	    int ich = NoentryChannelEE[Zside][iSC][gainId][0];
	    //	    EEDetId myEEDetId = EEDetId::unhashIndex(ich);
	    //	    EcalElectronicsId elecId = ecalElectronicsMap_.getElectronicsId(myEEDetId);
	    int DCCid = elecId[ich].dccId() - 1;  // 0-8
	    if(DCCid > 44)  DCCid -= 45;
	    if(NoentryDCCEE[Zside][DCCid][gainId] != 0) { // not yet printed
	      if(NoentryEE[Zside][iSC][gainId] == 25)
		cout << " No entry for " << NoentryEE[Zside][iSC][gainId] 
		     << " channels in SC " << iSC + 1 // iSC runs from 1 to 316
		     << " Z side " << Zside
		     << " gain " <<  gainValues[gainId]
		     << endl;
	      else if(NoentryEE[Zside][iSC][gainId] > 0) {
		int index = NoentryEE[Zside][iSC][gainId];
		cout << " No entry for " << index
		     << " channels in SC " << iSC + 1
		     << " Z side " << Zside
		     << " gain " <<  gainValues[gainId] << " channels";
		for (int ich = 0; ich < index; ich++)
		  cout << " " << NoentryChannelEE[Zside][iSC][gainId][ich];
		cout << endl;
	      }
	    } // NoentryDCCEE > 0
	  } // NoentryEE > 0
	} // loop on SC
      }  // loop on Zside
    }  // loop on gain
    for (int gainId = 0; gainId < 3; gainId++) {
      if(hPedestalEE[gainId]->GetEntries() > 0.) { // save only when filled!
	hGain[gainId]->Write();
	hGainEB[gainId]->Write();
	hGainEE[gainId]->Write();
	if(hEBDiff[gainId]->GetEntries() > 0.) hEBDiff[gainId]->Write();
	for(int side = 0; side < 2; side++)
	  if(hEEDiff[gainId][side]->GetEntries() > 0.) hEEDiff[gainId][side]->Write();
	hPedestalEE[gainId]->Write();
	hTotNoiseEE[gainId]->Write();
	hLFNoiseEE[gainId]->Write();
	hHFNoiseEE[gainId]->Write();
	for (int Zside = 0; Zside <2; Zside++) {
	  if(hEENoise[gainId][Zside]->GetEntries() > 0.) {
	    hEENoise[gainId][Zside]->Write();
	    hEENoiseLF[gainId][Zside]->Write();
	    hEENoiseHF[gainId][Zside]->Write();
	  }
	}
	hPedestalHitEE[gainId]->Write();
	hTotNoiseHitEE[gainId]->Write();
      }
    }
    if(SampleCorr_) {
      TH2F* hCorrMeanEE = new TH2F("CorrMeanEE","EE Time mean correlation",
				   10 , 0., 10. , 10, 0., 10.);
      TH2F* hCorrRMSEE = new TH2F("CorrRMSEE","EE Time RMS correlation",
				  10 , 0., 10. , 10, 0., 10.);
      for(int is1 = 0; is1 < 10; is1++) {
	for(int is2 = 0; is2 <= is1; is2++) {
	  double Mean = corrMean[is1][is2] / double(GoodEEChannels);
	  hCorrMeanEE->Fill(is1, is2, Mean);
	  if(is2 !=is1) hCorrMeanEE->Fill(is2, is1, Mean);
	  double RMS = corrRSM[is1][is2] / double(GoodEEChannels);
	  RMS = RMS - Mean * Mean;
	  RMS = sqrt(RMS);
	  hCorrRMSEE->Fill(is1, is2, RMS);
	  if(is2 !=is1) hCorrRMSEE->Fill(is2, is1, RMS);
	}
      }
      hCorrMeanEE->Write();
      hCorrRMSEE->Write();
    }  //  endcap data present
  }    // study sample correlation ?

  fPedestal.close();
  fDQM.close();
  fDB.close();
  f.Close();

  cout << "Exiting endJob" << endl;
}//endJob

//
// member functions
//

//========================================================================
void
EcalPedestal::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
//========================================================================
  ///SJ

  //std::cout<<"Entering analyse"<<std::endl;

  hitsAmplitudes_.clear();
  hitsEnergy_.clear();
  hitsThr_.clear();
  hitsEta_.clear();
  hitsPhi_.clear();
  nCrys_ = 0;
  

  // get geometry
  const CaloGeometry* geo = &iSetup.getData(geometryToken_);
  
   edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;

   //iEvent.getByToken(ebRecHitCollection_,barrelRecHitsHandle);
   //iEvent.getByToken(eeRecHitCollection_,endcapRecHitsHandle);

   try{
   iEvent.getByToken(ebRecHitCollection_,barrelRecHitsHandle);
   iEvent.getByToken(eeRecHitCollection_,endcapRecHitsHandle);
   }
   catch(std::exception& ex){
     //cout<<"No ECALRECHI collection found!!! "<<endl;
   }

   const EcalRecHitCollection* EBRecHits = nullptr;
   const EcalRecHitCollection* EERecHits = nullptr;
   //const EcalRecHitCollection* ESRecHits = nullptr;
  
   if ( !barrelRecHitsHandle.isValid() ){
     std::cout << "Error! EB rechits can't get product!" << std::endl;
   } else{
     EBRecHits = barrelRecHitsHandle.product();
   }

   if ( !endcapRecHitsHandle.isValid() ){
     std::cout << "Error! EE rechits can't get product!" << std::endl;
   } else{
     EERecHits = endcapRecHitsHandle.product();
   }

   //https://github.com/swagata87/OldLocalCovMiniAOD/blob/main/plugins/OldLocalCovMiniAOD.cc
   //edm::ESHandle<EcalPFRecHitThresholds> pThresholds;
   //iSetup.get<EcalPFRecHitThresholdsRcd>().get(pThresholds);
   //const EcalPFRecHitThresholds* thresholds = pThresholds.product();


   ////https://cmssdt.cern.ch/lxr/source/Calibration/HcalAlCaRecoProducers/plugins/AlCaHcalHBHEMuonProducer.cc#0150
   const EcalPFRecHitThresholds* thresholds = &iSetup.getData(tok_ecalPFRecHitThresholds_);

   //std::cout<<"Got handle to thresholds"<<std::endl;
  
  int ievt = iEvent.id().event();
  int lumi = iEvent.id().luminosityBlock();
  Timestamp evtime = iEvent.time();

  if(cnt_evt_ == 0) {
    firstReadEvent_ = ievt; 
    if(ECALType_ == "EB" || ECALType_ == "EA") {
      cout << " Barrel data : nb channels " << kEBChannels << endl;
    }
    else if(ECALType_ == "EE" || ECALType_ == "EA") {
      cout << " End cap data : nb channels " << kEEChannels << endl;
     }
    else {
      cout << " strange ECALtype : " << ECALType_ << " abort " << endl;
      return;
    }

    const auto& ecalElectronicsMap_ = iSetup.getData(EcalElectronicsMappingToken);
    for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
      EEDetId myEEDetId = EEDetId::unhashIndex(iChannel);
      elecId[iChannel] = ecalElectronicsMap_.getElectronicsId(myEEDetId);
    }

    const auto& chStatus = iSetup.getData(ChannelStatusToken);
    EcalChannelStatusMap::const_iterator chit;
    for (int iChannel = 0; iChannel < kEBChannels; iChannel++) {
      EBDetId id = EBDetId::unhashIndex(iChannel);
      chit = chStatus.getMap().find(id.rawId());
      if( chit != chStatus.getMap().end() ) {
	EcalChannelStatusCode ch_code = (*chit);
	uint16_t statusCode = ch_code.getStatusCode() & 31;
	//      if(statusCode == 1 || (statusCode > 7 && statusCode < 13))
	if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616: status 13 or 14 means bad data!
	  maskedChannels_.push_back(iChannel);
	//      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || (statusCode > 7 && statusCode < 13))
	if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	  badChannels_.push_back(iChannel);
      }
      for (int igain = 0; igain < kGains; igain++) 
	for (int iev = 0; iev < 1000; iev++) 
	  meanEBped[igain][iChannel][iev] = -999.;
    }
    // channel 46619 is very noisy but not in the tag (??)
    badChannels_.push_back(46619);
    for (int iChannel = 0; iChannel < kEEChannels; iChannel++) {
      EEDetId id = EEDetId::unhashIndex(iChannel);
      chit = chStatus.getMap().find(id.rawId());
      if( chit != chStatus.getMap().end() ) {
	EcalChannelStatusCode ch_code = (*chit);
	uint16_t statusCode = ch_code.getStatusCode() & 31;
	//      if(statusCode == 1 || (statusCode > 7 && statusCode < 13))
	if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	  maskedEEChannels_.push_back(iChannel);
	//      if(statusCode == 1 || statusCode == 3 || statusCode == 4 || (statusCode > 7 && statusCode < 13))
	if(statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7)         //  170616
	  badEEChannels_.push_back(iChannel);
      }
      for (int igain = 0; igain < kGains; igain++) 
	for (int iev = 0; iev < 1000; iev++) 
	  meanEEped[igain][iChannel][iev] = -999.;
    }
    /*int NbOfmaskedChannels =  maskedChannels_.size();
    int i = 0;
    //cout << " Nb masked EB channels (statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7) " << NbOfmaskedChannels << endl;
    for (vector<int>::iterator iter = maskedChannels_.begin(); iter != maskedChannels_.end(); ++iter) {
      cout << *(iter) << " ";
      i++;
      if(i == 20) {
	cout << endl;
	i = 0;
      }
    }
    cout << endl;
    */
    /*NbOfmaskedChannels =  maskedEEChannels_.size();
    i = 0;
    //cout << " Nb masked EE channels (statusCode == 1 || statusCode == 3 || statusCode == 4 || statusCode > 7) " << NbOfmaskedChannels << endl;
    for (vector<int>::iterator iter = maskedEEChannels_.begin(); iter != maskedEEChannels_.end(); ++iter) {
      cout << *(iter) << " ";
      i++;
      if(i == 20) {
	cout << endl;
	i = 0;
      }
    }
    cout << endl;
    */
  }
  cnt_evt_++;
  int expectedGain = 0;
  //  int expectedGain = 0, previousExpectedGain = 0;

  if(runType_ == "Pedes") {  
    //   expectedGain = 4 - (cnt_evt_%3);
    //   if(expectedGain == 4) expectedGain = 1;
    expectedGain = 3;
    if(cnt_evt_ > 300) expectedGain = 2;
    if(cnt_evt_ > 600) expectedGain = 1;
  }
  else if(runType_ == "Ped_6") {
    expectedGain = 2;
  }
  else if(runType_ == "Other" || runType_ == "Calib") {  
    expectedGain = 1;
  }
  else {
    cout << " strange runType : " << runType_ << " abort " << endl;
    return;
  }

  //cout << "Running on event =" << cnt_evt_ << endl;
  if(cnt_evt_ < startevent_) return;
  //else cout << "starting to analyze" << endl;
  //  if(cnt_evt_%100 ==1) 
  //  unsigned int vtime = evtime.value();
  time_t unix = evtime.unixTime();
  char buf[256];
  strftime(buf, sizeof(buf), "%F %R:%S", gmtime(&unix));
  //cout << "Running on event = " << ievt << " at " << cnt_evt_ << " (UTC) " << buf << "." << evtime.microsecondOffset() << " lumi block " << lumi << endl;
  //  cout << " Unix time " << unix << " (UTC) " << buf << " micro " << evtime.microsecondOffset() << endl;
  //  if(ievt == 1) begtime[expectedGain - 1] = evtime;
  //  if(ievt == 300) endtime[expectedGain - 1] = evtime;
  time_t unixbeg = begtime[expectedGain - 1].unixTime();
  time_t unixend = endtime[expectedGain - 1].unixTime();
  time_t micbeg = begtime[expectedGain - 1].microsecondOffset();
  time_t micend = endtime[expectedGain - 1].microsecondOffset();
  if(unix < unixbeg || (unix == unixbeg && evtime.microsecondOffset() < micbeg)) begtime[expectedGain - 1] = evtime;
  if(unix > unixend || (unix == unixend && evtime.microsecondOffset() > micend)) endtime[expectedGain - 1] = evtime;
  //  if(ievt < begtime[expectedGain - 1]) begtime[expectedGain - 1] = evtime;
  //  if(ievt > endtime[expectedGain - 1]) endtime[expectedGain - 1] = evtime;

  //GET THE DIGIS
  int nebd = 0;
  int need = 0;
  const EBDigiCollection* EBdigis = 0;
  const EEDigiCollection* EEdigis = 0;
  Handle< EBDigiCollection > pEBdigis;
  try {
    //    iEvent.getByLabel( digiProducer_, EBdigiCollection_, pEBdigis);
    iEvent.getByToken(EBDigisToken, pEBdigis);
    EBdigis = pEBdigis.product(); // get a ptr to the product
    nebd = EBdigis->size();
  } catch ( std::exception& ex ) {
    cerr << "Error! can't get the product ebDigis" << endl;
    if(ECALType_ == "EB" || ECALType_ == "EA") 
      cout << "Event = " << cnt_evt_ << " EcalPedestal ebDigis not available" << endl;
    
  } //getting barrel digis

    // endcaps
  Handle< EEDigiCollection > pEEdigis;
  try {
    //    iEvent.getByLabel( digiProducer_, EEdigiCollection_, pEEdigis);
    iEvent.getByToken(EEDigisToken, pEEdigis);
    EEdigis = pEEdigis.product(); // get a ptr to the product
    need = EEdigis->size();
  } catch ( std::exception& ex ) {
    cerr << "Event = " << cnt_evt_ << "Error! can't get the product eeDigis" << endl;
    /*if(ECALType_ == "EE" || ECALType_ == "EA") 
      cout << "Event = " << cnt_evt_ << " EcalPedestal eeDigis not available" << endl;
    */
  } //getting endcap digis


  std::vector<std::vector<double>> hitTimeSamples;
  
  int xEvent = ievt;
  /*if(nebd != 60740 || need != 14523)           // number of read channels in 2018
    cout << "Number of crystals in barrel " << nebd << " endcaps " << need << endl;
  */

  if(ECALType_ == "EB" || ECALType_ == "EA") {    // barrel

    // Loop over Ecal barrel digis
    int CrystalinSM[36];
    for(int SM = 0; SM < 36; SM++) CrystalinSM[SM] = 0;
    int SM = -999;
    if(nebd != 0) {
      nevent[0]++;
      int wronggainId = 0, largediff = 0;
      int wronggainCh[10] = {10 * -999};
      for (EBDigiCollection::const_iterator digiItr = EBdigis->begin(); digiItr != EBdigis->end(); ++digiItr) {

	//std::cout<<"inside EB digis "<<std::endl;
	// Get (ieta, iphi) coordinates
	int ieta = EBDetId((*digiItr).id()).ieta();
	int iphi = EBDetId((*digiItr).id()).iphi();
	SM = EBDetId((*digiItr).id()).ism();      // Get SM number (from 1 to 36)
	int iChannel = EBDetId((*digiItr).id()).hashedIndex();      // here iChannel runs from 0 to 61200	
	CrystalinSM[SM - 1]++;

	////////////SJ
	DetId id = (*digiItr).id();
	float rhThres = -99.0;
	if (thresholds != nullptr) {
	  rhThres = (*thresholds)[id];  // access PFRechit thresholds for noise cleaning
	}


	//std::cout<<"Threhsold is "<<rhThres<<std::endl;

	double rechitEn = -99;
	EcalRecHitCollection::const_iterator it = EBRecHits->find(id);
	
	if (it != EBRecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalBarrel && it->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}

	//cout<<"rhEn "<<rechitEn<<endl;
	if(rechitEn==0) continue;
	
	///eta phi
	EBDetId det = it->id();
	//cout<<"getting Global point"<<endl;
	//cout<<"det ieta iphi "<<ieta<<" "<<iphi<<endl;
	///https://cmssdt.cern.ch/lxr/source/RecoMET/METFilters/plugins/EcalBadCalibFilter.cc

	double rheta = -99;
	double rhphi = -99;
	if (! (it == barrelRecHitsHandle->end()) ){ 
	  const GlobalPoint & rechitPoint = geo->getPosition(det);
	  //std::cout<<"Printing the values "<<endl;
	  
	  
	  double rheta = rechitPoint.eta();
	  double rhphi = rechitPoint.phi();
	  
	  //std::cout<<"Eta : Phi "<<rheta<<" "<<rhphi<<std::endl;
	}

	////////////SJ
	
	EBDataFrame df( *digiItr );
	int gainId = df.sample(0).gainId();
	if (gainId < kFirstGainId || gainId >= kFirstGainId + kGains) {
	  /*cout << "EcalPedestal::analyze: Warning: gainId = " << gainId
	       << " channel " << iChannel
	    //	     << " at eta = " << channelPos.first
	    //	     << " and phi = " << channelPos.second
	       << " eta = " << ieta
	       << " phi = " << iphi
	       << " for event " << cnt_evt_ << endl;
	  */
	  continue;
	}
	// gain check on the first read channel whatever it is
	if(digiItr == EBdigis->begin()) {
	  expectedGain = gainId;
	  /*if(ievt == 1) {
	    cout << " ********** gain at event " << cnt_evt_ << " " << gainValues[expectedGain - 1] << endl;
	    //	    previousExpectedGain = expectedGain;
	  }
	  */
	  //	  else if(expectedGain != previousExpectedGain) {
	  //	    cout << " ********** gain change at event " << cnt_evt_ << " Now : " << gainValues[expectedGain - 1] << endl;
	  //	    previousExpectedGain = expectedGain;
	  //	  }
	  if(gainId == 1) {             // gain 12
	    nevent[1]++;
	    if(xEvent < 1 || xEvent > gain12events_) cout << " **** gain 12 strange event " << ievt << " xEvent " << xEvent << endl;
	  }
	  if(gainId == 2) {             // gain 6
	    nevent[2]++;
	    if(!NormalSequence_) xEvent = ievt + gain12events_;
	    if(xEvent < gain12events_ || xEvent > gain6events_) cout << " **** gain 6 strange event " << ievt << " xEvent " << xEvent << endl;
	  }
	  if(gainId == 3) {             // gain 1
	    nevent[3]++;
	    if(!NormalSequence_) xEvent = ievt + gain6events_;
	    if(xEvent < gain6events_ || xEvent > gain1events_) cout << " **** gain 1 strange event " << ievt << " xEvent " << xEvent << endl;
	  }
	  //cout << " event " << cnt_evt_ << " gain " << expectedGain << " nevent " << nevent[expectedGain] << endl;
	  hGains->Fill(xEvent, gainValues[gainId - 1]);
	}
	if (gainId != expectedGain) {
	  if(wronggainId < 10) wronggainCh[wronggainId] = iChannel;
	  wronggainId++;
	  /*if (nevent[expectedGain] == 1)  // print only once
	    cout << "EcalPedestal::analyze: Warning: gainId = " << gainId
		 << " while " << expectedGain << " was expected " 
		 << " channel " << iChannel
		 << " eta = " << ieta
		 << " phi = " << iphi
		 << " for event " << ievt << " at " << cnt_evt_ << endl;
	  */
	  continue;
	}
	//no less than 10 samples readout
	int nSample = digiItr->size();
	if (nSample != 10)   {
	  /*cout << "EcalPedestal::analyze: Warning: N samples " << nSample
	       << " channel " << iChannel
	    //	     << " at eta = " << channelPos.first
	    //	     << " and phi = " << channelPos.second
	       << " eta = " << ieta
	       << " phi = " << iphi
	       << " for event " << cnt_evt_ << endl;
	  */
	  continue;
	}
	hEBEntry[gainId - 1]->Fill(iphi, ieta, 1.);

	int arrayId = (expectedGain - kFirstGainId) * kEBChannels + iChannel;
	Double_t adc[10];
	float pedestalevent = 0.;
	double pedestalevent2 = 0.;
	Double_t adcmin = 5000.;
	Double_t adcmax = -.1;

	///SJ
	std::vector<double> amplSamples(10);
	
	//	if(iChannel > 7814 && iChannel < 7820) cout << " Xtal " << iChannel << " ADC:";
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  int gainSample = df.sample(iSample).gainId();
	  if(gainSample != gainId) {
	    /*cout << "EcalPedestal::analyze: Warning: gain change for sample " << iSample
		 << " gain " << gainSample
		 << " channel " << iChannel
		 << " eta = " << ieta
		 << " phi = " << iphi
		 << " for event " << cnt_evt_ << endl;
	    */
	  }
	  adc[iSample] = float(df.sample(iSample).adc());
	  //	  if(iChannel > 7814 && iChannel < 7820) cout << " " << adc[iSample];
	  
	  pedestalevent += adc[iSample];
	  pedestalevent2 += adc[iSample] * adc[iSample];
	  if(SampleCorr_ && (gainId == 1)) (*vCovMatEB)[iChannel]->Fill(iSample,adc[iSample]);
	  if(adc[iSample] < adcmin) adcmin = adc[iSample];
	  if(adc[iSample] > adcmax) adcmax = adc[iSample];

	  ///SJ
	  amplSamples[iSample] = adc[iSample];
	  
	}  // sample
	float avg = pedestalevent / nSample;

	///SJ
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  amplSamples[iSample] = amplSamples[iSample] - avg;
	}
	hitTimeSamples.push_back(amplSamples);
	nCrys_++;
	  
	//	if(cnt_evt_ > 8072) cout << " event " << cnt_evt_ << " EB channel " << iChannel << " mean " << avg << endl;
	nEntriesStandard_[arrayId] += 10;
	pedestalStandard_[arrayId] += pedestalevent;
	pedestal2Standard_[arrayId] += pedestalevent2;
	//	if(iChannel > 7814 && iChannel < 7820) cout << " mean " << avg << " sum ped " << pedestalStandard_[arrayId]
	//						    << " nb entries " << nEntriesStandard_[arrayId] << endl;

	avgSamples_[arrayId] += avg;
	avgSample2_[arrayId] += avg * avg;
	nEntriesAvg_[arrayId]++;
	// vector<int>::iterator result = find(badChannels_.begin(), badChannels_.end(), iChannel);
	/*
	if (result == badChannels_.end()) {    // do not use bad channels
	  if(nevent[expectedGain] == 1) {
	    meanEBped[expectedGain][iChannel] = avg;
	    //	    if(iChannel%10000 == 0) cout << "EB event " << xEvent << " channel " << iChannel << " gain " << expectedGain 
	    //	    << " setting mean " << meanEBped[iChannel] << endl;
	  }
	  float diffped = meanEBped[expectedGain][iChannel] - avg;
	  //	  if(cnt_evt_ > 8072) cout << " event " << cnt_evt_ << " EB channel " << iChannel << " gain " << expectedGain
	  //				   << " 1st mean " << meanEBped[expectedGain][iChannel] << " diffped " << diffped << endl;
	  hGain[expectedGain - 1]->Fill(ievt, diffped);
	  hGainEB[expectedGain - 1]->Fill(ievt, diffped);
	  //	  if(iChannel%10000 == 0) cout << "EB event " << xEvent << " channel " << iChannel << " gain " << expectedGain
	  //				      << " 1st mean " << meanEBped[iChannel] << " diffped " << diffped << endl;
	  if(abs(diffped) > 10.) {
	    //	    cout << "******** EB event " << xEvent << " gain " << expectedGain << " channel " << iChannel
	    //		 << " diffped " << diffped << " phi " << iphi << " eta " << ieta << endl;
	    hEBDiff[expectedGain - 1]->Fill(iphi, ieta, 1.);
	    largediff++;
	  }
	}
	*/

	if(ievt < 1000)
	  meanEBped[expectedGain - 1][iChannel][ievt - 1] = avg;
	else if(iChannel == 1) cout << "  *** problem in event " << cnt_evt_ << " ievt = " << ievt << endl;
	//      if(cnt_evt_ > 400 && iChannel == 247) {

	// Compute HF noise
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  Double_t adcHF = adc[iSample] - avg;
	  pedavg_[arrayId] += adcHF;
	  pedavg2_[arrayId] += adcHF * adcHF;
	  hPedsample[gainId - 1]->Fill(iSample,adcHF);
	} 

	//	if(cnt_evt_ > 8072) cout << " event " << cnt_evt_ << " EB channel " << iChannel << " arrayid " << arrayId << endl;

	nEntriesSingle_[arrayId]++;
	pedestalSingle_[arrayId] += adc[6];
	pedestal2Single_[arrayId] += adc[6] * adc[6];
	//	if(cnt_evt_ > 8072) cout << " event " << cnt_evt_ << " EB channel " << iChannel << " end loop " << endl;

	///SJ
	hitsThr_.push_back(rhThres);
	hitsEnergy_.push_back(rechitEn);      
	hitsEta_.push_back(rheta);
	hitsPhi_.push_back(rhphi);
	//std::cout<<"thres and en "<<rhThres<<" "<<rechitEn<<std::endl;
      }  // loop over digis

      
    if (wronggainId > 0)  {    // print how many channels have wrong gain
      //cout << "Event = " << ievt << " at " << cnt_evt_ << " wrong gain on " << wronggainId
      //     << " channels " ;

	int totCh = wronggainId;
	if(totCh > 10) totCh = 10;
	for(int iCh = 0; iCh < totCh; iCh++){
	  //cout << wronggainCh[iCh] << " ";
	  double idiotVar = wronggainCh[iCh]; ///so that it compiles without error that wronggainCh is not used
	
	}
	// cout << endl;
    }

    /*if (largediff > 0)  // print how many channels have large difference
	cout << "Event = " << ievt << " at " << cnt_evt_ << " EB large ped difference " << largediff
	     << " channels" << endl;
    */
      for(int SM = 0; SM < 36; SM++) {
	/*if(cnt_evt_ == 1)
	  cout << " SM " << SM + 1 << " Number of read channels " << CrystalinSM[SM] << endl;
	*/
	hEntries->Fill(SM,CrystalinSM[SM]);
      }
    
      /// SJ
      hitsAmplitudes_ = hitTimeSamples;

      
      ///Fill the tree
      tree->Fill();
   
    }  //  Barrel digis present

   
    else
      cout << " No EB digi in event " << cnt_evt_ << endl;
    
  }

  // Loop over Ecal endcap digis
  if(ECALType_ == "EE" || ECALType_ == "EA") {    // endcap
    if(need != 0) {
      nevent[5]++;
      int largediff = 0;
      for (EEDigiCollection::const_iterator digiItr = EEdigis->begin(); digiItr != EEdigis->end(); ++digiItr) {
	int iChannel = EEDetId((*digiItr).id()).hashedIndex();    // here iChannel runs from 0 to 14647
	if(iChannel >= kEEChannels) {
	  cout << " ****** Error ****** EE channel " << iChannel << endl;
	  continue;
	}
	int iz = EEDetId((*digiItr).id()).zside();
	int ix = EEDetId((*digiItr).id()).ix();
	int iy = EEDetId((*digiItr).id()).iy();
	
	
	////////////SJ
	DetId id = (*digiItr).id();
	float rhThres = -99.0;
	if (thresholds != nullptr) {
	  rhThres = (*thresholds)[id];  // access PFRechit thresholds for noise cleaning
	}
	


	double rechitEn = -99;
	EcalRecHitCollection::const_iterator it = EERecHits->find(id);
	
	if (it != EERecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalBarrel && it->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}

	if(rechitEn==0) continue; 

	///eta phi
	EEDetId det = it->id();

	double rheta = -99;
	double rhphi = -99; 

	if (! (it == endcapRecHitsHandle->end()) ){  
	  const GlobalPoint & rechitPoint = geo->getPosition(det);
	  rheta = rechitPoint.eta();
	  rhphi = rechitPoint.phi();
	}

	////////////SJ

	
	EEDataFrame df( *digiItr );
	int gainId = df.sample(0).gainId();
	if (gainId < kFirstGainId || gainId >= kFirstGainId + kGains) {
	  cout << "EE crystal "<< iChannel  << " side " << iz << " x " << ix << " y " << iy
	       << " gain " << gainId  << " for event " << cnt_evt_ << endl;
	  continue;
	}
	if(digiItr == EEdigis->begin()) {      // gain check on the first read channel whatever it is
	  expectedGain = gainId;
	  if(gainId == 1) {             // gain 12
	    nevent[6]++;
	    if(xEvent < 1 || xEvent > gain12events_) cout << " **** gain 12 strange EE event " << ievt << " xEvent " << xEvent << endl;
	  }
	  if(gainId == 2) {             // gain 6
	    nevent[7]++;
	    if(!NormalSequence_) xEvent = ievt + gain12events_;
	    if(xEvent < gain12events_ || xEvent > gain6events_) cout << " **** gain 6 strange EE event " << ievt << " xEvent " << xEvent << endl;
	  }
	  if(gainId == 3) {             // gain 1
	    nevent[8]++;
	    if(!NormalSequence_) xEvent = ievt + gain6events_;
	    if(xEvent < gain6events_ || xEvent > gain1events_) cout << " **** gain 1 strange EE event " << ievt << " xEvent " << xEvent << endl;
	  }
	}
	if (gainId != expectedGain) {
	  cout << "EcalPedestal::analyze: Warning: gainId = " << gainId
	       << " while " << expectedGain << " was expected " 
	       << " channel " << iChannel
	       << " at X = " << ix << " and Y = " << iy
	       << " for event " << cnt_evt_ << endl;
	  continue;
	}
	//no less than 10 samples readout
	int nSample = digiItr->size();
	if (nSample != 10)   {
	  cout << "EcalPedestal::analyze: Warning: N samples " << nSample
	       << " channel " << iChannel
	       << " side " << iz << " ix " << ix << " iy " << iy
	       << " for event " << cnt_evt_ << endl;
	  continue;
	}

	int arrayId = (expectedGain - kFirstGainId) * kEEChannels + iChannel; 
	float adc[10];
	float pedestalevent = 0.;
	double pedestalevent2 = 0.;

	///SJ
	std::vector<double> amplSamples(10);

	
	//	if(iChannel == 10000) cout << " Xtal " << iChannel << " ADC:";
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  int gainSample = df.sample(iSample).gainId();
	  if(gainSample != gainId) {
	    cout << "EcalPedestal::analyze: Warning: gain change for sample "
		 << " EE channel " << iChannel
		 << iSample << " gain " << gainSample
		 << " side " << iz << " ix " << ix << " iy " << iy
		 << " for event " << cnt_evt_ << endl;
	  }
	  adc[iSample] = float(df.sample(iSample).adc());
	  //	  if(iChannel == 10000) cout << " " << adc[iSample];
	  pedestalevent += adc[iSample];
	  pedestalevent2 += adc[iSample] * adc[iSample];

	  ///SJ
	  amplSamples[iSample] = adc[iSample];
	  
	  if(SampleCorr_ && (gainId == 1)) (*vCovMatEE)[iChannel]->Fill(iSample,adc[iSample]);
	}  // sample
	float avg = pedestalevent / nSample;
	//	if(iChannel == 10000) cout << " mean " << avg << endl;

	///SJ
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  amplSamples[iSample] = amplSamples[iSample] - avg;
	}
	hitTimeSamples.push_back(amplSamples);
	nCrys_++;
	  
	nEntriesStandardEE_[arrayId] += 10;
	pedestalStandardEE_[arrayId] += pedestalevent;
	pedestal2StandardEE_[arrayId] += pedestalevent2;

	avgSamplesEE_[arrayId] += avg;
	avgSample2EE_[arrayId] += avg * avg;
	nEntriesAvgEE_[arrayId]++;
	// vector<int>::iterator result = find(badEEChannels_.begin(), badEEChannels_.end(), iChannel);
	/*
	if (result == badEEChannels_.end()) {    // do not use bad channels
	  if(nevent[expectedGain + 5] == 1) {
	    meanEEped[expectedGain][iChannel] = avg;   // nevent[5,6,7] for gain 1,2,3
	    //	    if(iChannel%1000 == 0) cout << "EE event " << xEvent  << " channel " << iChannel << " gain " << expectedGain 
	    //					<< " setting mean " << meanEEped[expectedGain][iChannel] << endl;
	  }
	  float diffped = meanEEped[expectedGain][iChannel] - avg;
	  hGain[expectedGain - 1]->Fill(ievt, diffped);
	  hGainEE[expectedGain - 1]->Fill(ievt, diffped);
	  //	  if(iChannel%1000 == 0) cout << "EE event " << xEvent  << " channel " << iChannel << " gain " << expectedGain
	  //				      << " 1st mean " << meanEEped[expectedGain][iChannel] << " diffped " << diffped << endl;
	  //	  if(abs(diffped) > 10.) cout << "******* EE event " << xEvent << " gain " << expectedGain << " channel " << iChannel
	  //				     << " diffped " << diffped << endl;
	  if(abs(diffped) > 10.) {
	    //	    cout << "******* EE event " << xEvent << " gain " << expectedGain << " channel " << iChannel
	    //		 << " diffped " << diffped << " iz " << iz << " ix " << ix << " iy " << iy << endl;
	    int izz = iz;
	    if(iz == -1) izz = 0;
	    hEEDiff[expectedGain - 1][izz]->Fill(ix, iy, 1.);
	    largediff++;
	  }
	}
	*/

	if(ievt < 1000)
	  meanEEped[expectedGain - 1][iChannel][ievt - 1] = avg;
	else if(iChannel == 1) cout << "  *** problem in event " << cnt_evt_ << " ievt = " << ievt << endl;

	// Compute HF noise
	for (int iSample = 0; iSample < nSample; ++iSample) {
	  float adcHF = adc[iSample] - avg;
	  pedavgEE_[arrayId] += adcHF;
	  pedavg2EE_[arrayId] += adcHF * adcHF;
	  hPedsample[gainId - 1]->Fill(iSample,adcHF);
	}

	///SJ
	hitsThr_.push_back(rhThres);
	hitsEnergy_.push_back(rechitEn);      

	hitsEta_.push_back(rheta);
	hitsPhi_.push_back(rhphi);
	
      }  // loop over digis

      /// SJ
      hitsAmplitudes_ = hitTimeSamples;

      
      ///Fill the tree
      tree->Fill();
      
      if (largediff > 0)  // print how many channels have large difference
	cout << "Event = " << ievt << " at " << cnt_evt_ << " EE large ped difference " << largediff
	     << " channels" << endl;
    }//  End cap digis present
    else
      cout << " No EE digi in event " << cnt_evt_ << endl;
  }

  //GET THE UNCALIB RECHITS
  if(ECALType_ == "EB" || ECALType_ == "EA") {
    Handle< EBUncalibratedRecHitCollection > pEBUncalibRecHits;
    const EBUncalibratedRecHitCollection*  EBuncalibRecHits = 0; 
    try {
      //      iEvent.getByLabel( EBuncalibRecHitCollection_, pEBUncalibRecHits );
      iEvent.getByToken(EBUncalibRecHitToken, pEBUncalibRecHits);
      EBuncalibRecHits = pEBUncalibRecHits.product(); // get a ptr to the product

      EBDetId EBHitId(0); 
      //      int nunhits = EBuncalibRecHits->size();
      //      if(nunhits != kEBChannels) cout << " Uncalib hits size " << nunhits << endl;
      //      if(nunhits != 61000) cout << " Uncalib hits size " << nunhits << endl;
      for(EBUncalibratedRecHitCollection::const_iterator ithit = EBuncalibRecHits->begin(); ithit != EBuncalibRecHits->end(); ++ithit) {
	EBHitId = ithit->id();
	//	int xtalNumber = EBHitId.ic() - 1; 	// here xtalNumber runs from 0 to 1699
	//	if(xtalNumber < 0 || xtalNumber > 1699) cout << " Crystal " << xtalNumber << endl;
	int iChannel = EBHitId.hashedIndex();	// here iChannel runs from 0 to 61200
	float HitAmp = ithit->amplitude();
	int arrayId = (expectedGain - kFirstGainId) * kEBChannels + iChannel; 
	//      if(xtalNumber < 10) cout  << " Crystal " << xtalNumber << " id " << arrayId 
	//				<< " amplitude hit " << HitAmp << endl;
	nEntriesHit_[arrayId]++;
	pedestalHit_[arrayId] += HitAmp;
	pedestal2Hit_[arrayId] += HitAmp * HitAmp;
      }//loop uncalib hits
    } catch ( std::exception& ex ) {
      // std::cerr << "Error! can't get the product " << EBuncalibRecHitCollection_ << endl;
      // cout << "Event = " << cnt_evt_  << "EcalPedestal error! can't get the product " << EBuncalibRecHitCollection_ << endl;
    }
    // second method MaxSample
    Handle< EBUncalibratedRecHitCollection > pEBUncalibRecHitsMax;
    const EBUncalibratedRecHitCollection*  EBuncalibRecHitsMax = 0; 
    try {
      //      iEvent.getByLabel( edm::InputTag("ecalMaxSampleUncalibRecHit:EcalUncalibRecHitsEB"), pEBUncalibRecHitsMax );
      iEvent.getByToken(EBUncalibRecHitMaxToken, pEBUncalibRecHitsMax);
      EBuncalibRecHitsMax = pEBUncalibRecHitsMax.product(); // get a ptr to the product

      EBDetId EBHitId(0); 
      for(EBUncalibratedRecHitCollection::const_iterator ithit = EBuncalibRecHitsMax->begin(); ithit != EBuncalibRecHitsMax->end(); ++ithit) {
	EBHitId = ithit->id();
	int iChannel = EBHitId.hashedIndex();  	// here iChannel runs from 0 to 61199
	float HitAmp = ithit->amplitude();
	int arrayId = (expectedGain - kFirstGainId) * kEBChannels + iChannel;
	nEntriesMinMax_[arrayId]++;
	pedestalMinMax_[arrayId] += HitAmp;
	pedestal2MinMax_[arrayId] += HitAmp * HitAmp;
      }//loop uncalib hits
    } catch ( std::exception& ex ) {
      // cout << "Event = " << cnt_evt_ << "EcalPedestal error! can't get ecalMaxSampleUncalibRecHit:EcalUncalibRecHitsEB" << endl;
    }
    // third method AlphaBetaFit
    Handle< EBUncalibratedRecHitCollection > pEBUncalibRecHitsFit;
    const EBUncalibratedRecHitCollection*  EBuncalibRecHitsFit = 0; 
    try {
      //      iEvent.getByLabel( edm::InputTag("ecalFixedAlphaBetaFitUncalibRecHit:EcalUncalibRecHitsEB"), pEBUncalibRecHitsFit );
      iEvent.getByToken(EBUncalibRecHitMaxToken, pEBUncalibRecHitsFit);
      EBuncalibRecHitsFit = pEBUncalibRecHitsFit.product(); // get a ptr to the product

      EBDetId EBHitId(0); 
      for(EBUncalibratedRecHitCollection::const_iterator ithit = EBuncalibRecHitsFit->begin(); ithit != EBuncalibRecHitsFit->end(); ++ithit) {
	EBHitId = ithit->id();
	int iChannel = EBHitId.hashedIndex();  	// here iChannel runs from 0 to 61199
	float HitAmp = ithit->amplitude();
	if(runType_ == "Ped_6") {  // gain 6
	  //	  if((iChannel < 3240 || iChannel > 3244) && HitAmp > 1000.) 
	  //	    cout << " Channel " << iChannel << " event " << ievt << " Fit amplitude " << HitAmp << endl;
	}
	else {  // gain 12
	  if(iChannel != 54145 && HitAmp > 4000.) 
	    cout << " Channel " << iChannel << " event " << ievt << " Fit amplitude " << HitAmp << endl;
	}
      }//loop uncalib hits
    } catch ( std::exception& ex ) {
      // cout << "Event = " << cnt_evt_ << "EcalPedestal error! can't get ecalFixedAlphaBetaFitUncalibRecHit:EcalUncalibRecHitsEB" << endl;
    }
  }  // barrel

  if(ECALType_ == "EE" || ECALType_ == "EA") { // endcaps
    Handle< EEUncalibratedRecHitCollection > pEEUncalibRecHits;
    const EEUncalibratedRecHitCollection*  EEuncalibRecHits = 0; 
    try {
      //      iEvent.getByLabel( EEuncalibRecHitCollection_, pEEUncalibRecHits );
      iEvent.getByToken(EEUncalibRecHitToken, pEEUncalibRecHits);
      EEuncalibRecHits = pEEUncalibRecHits.product(); // get a ptr to the product

      EEDetId EEHitId(0); 
      //      int nunhits = EEuncalibRecHits->size();
      //      if(nunhits != kEEChannels) cout << " Uncalib hits size " << nunhits << endl;
      for(EEUncalibratedRecHitCollection::const_iterator ithit = EEuncalibRecHits->begin(); ithit != EEuncalibRecHits->end(); ++ithit) {
	EEHitId = ithit->id();
	int iChannel = EEHitId.hashedIndex();    // here iChannel runs from 0 to 14647
	if(iChannel >= kEEChannels) {
	  cout << " ****** Error ****** EE channel " << iChannel << endl;
	  continue;
	}
	//	int iz = EEHitId.zside();
	//	int ix = EEHitId.ix();
	//	int iy = EEHitId.iy();
	//	if(ix%10 == 0 && iy%10 == 0) cout << " side " << iz << " ix " << ix << " iy " << iy;
	float HitAmp = ithit->amplitude();
	//	if(ix%10 == 0 && iy%10 == 0) cout << " amplitude " << HitAmp << endl;
	int arrayId = (expectedGain - kFirstGainId) * kEEChannels + iChannel; 
	nEntriesHitEE_[arrayId]++;
	pedestalHitEE_[arrayId] += HitAmp;
	pedestal2HitEE_[arrayId] += HitAmp * HitAmp;
      }//loop uncalib hits
    } catch ( std::exception& ex ) {
      // std::cerr << "Error! can't get the product " << EEuncalibRecHitCollection_ << endl;
      // cout << "Event = " << cnt_evt_ << "EcalPedestal error! can't get the product " << EEuncalibRecHitCollection_ << endl;
    }
  }  // endcap

  //uncalib rechits


  //cout << "Exiting analyze" << endl;
}//analyze

//define this as a plug-in
DEFINE_FWK_MODULE( EcalPedestal );
