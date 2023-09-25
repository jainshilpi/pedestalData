//
// File H4CovMat.cxx
//
/*!
 * \class H4CovMat
 * \brief Class used to compute the mean, variance and correlation between
 * samples. 
 *
 * $Date: 2005/05/10 15:16:23 $
 * $Author: brunelie $
 *
 * 
*/

// C & C++ standard library include files
#include <iostream>
#include <cstdlib>
#include "math.h"

// Class definition
//#include "Pit/H4CovMat/interface/H4CovMat.h"
#include "pedestalData/H4CovMat/interface/H4CovMat.h"

using namespace std;

//! Default constructor needed by root
H4CovMat::H4CovMat()
{
  ndim_ = 0;
  nevt_ = 0;
  mode_ = true;
  isStationarity_ = false;
  isLFsub_ = false;
  vmean_ = new vector<float>;
  vcov_ = new vector<float>;
  vsum_ = new vector<float>;
  vevt_ = new vector<float>;
}

//! Destructor
H4CovMat::~H4CovMat()
{
  delete vmean_;
  delete vcov_;
  delete vsum_;
  delete vevt_;
}

//! Assign
H4CovMat& H4CovMat::operator=(const H4CovMat& m)
{
  this->mode_ = m.mode_;
  this->isStationarity_ = m.isStationarity_;
  this->isLFsub_ = m.isLFsub_;
  this->ndim_ = m.ndim_;
  this->nevt_ = m.nevt_;
  this->curXtal1_ = m.curXtal1_;
  *(this->vmean_) = *(m.vmean_); // copy the vectors
  *(this->vcov_) = *(m.vcov_);
  this->lfrms_ = m.lfrms_;
  this->hfrms_ = m.hfrms_;
  this->smCrystalToIndex_ = m.smCrystalToIndex_; // copy the map
  if (this->mode_) {
    *(this->vsum_) = *(m.vsum_);
    *(this->vevt_) = *(m.vevt_);
    this->ComputeEstimators();
  }
  return *this;
}

// We have to parse the [CovMat] section to figure out whether we
// want to read or build a covariance matrix
bool H4CovMat::init(int run)
{
  ndim_ = 10;
  /*
  if (gConfigParser->isDefined("CovMat::Stationarity"))
    if (gConfigParser->readIntOption("CovMat::Stationarity") > 0)
      isStationarity_ = true;
    else
      isStationarity_ = false;

  if (gConfigParser->isDefined("CovMat::LFSubtraction"))
    if (gConfigParser->readIntOption("CovMat::LFSubtraction") > 0)
      isLFsub_ = true;
    else
      isLFsub_ = false;
  */
  isStationarity_ = false;
  isLFsub_ = false;

  mode_ = true;
  if (isStationarity_) {
    (*vmean_).resize(1);
    (*vcov_).resize(ndim_);
    (*vsum_).resize(1);
  } else {
    (*vmean_).resize(ndim_ + 1);
    (*vcov_).resize(ndim_*(ndim_ + 1)/2 + 1);
    (*vsum_).resize(ndim_);
  }
  for (unsigned int i = 0; i < vsum_->size(); i++)
    (*vsum_)[i] = 0.;

  /*
  // Look at config file to know if we read an input covariance matrix
  if (!gConfigParser->isDefined("CovMat::FileName")) return true;

  // Attempt to get the covmat file name
  string covFileName, covTreeName;
  cout << "H4CovMat::init Opening CovMat file!" << endl;
  covFileName = gConfigParser->readStringOption("CovMat::FileName");
  
  // Now, try to figure out the name of the tree
  try {
    covTreeName = gConfigParser->readStringOption("CovMat::TreeName") ;
  } catch(const char* e) {
    cout << "Warning: CovMat::TreeName not set, using the default "
	 << "(C01)" << endl ;
    covTreeName = "C01";
  }
  
  // Open input root file, tree and branch
  covmatFile_ = new TFile(covFileName.c_str(), "READ");
  if(!covmatFile_) 
    cout << "H4CovMat::init: ERROR, unable to open covmat file "
	 << covFileName << " for reading" << endl;

  covmatTree_ = (TTree*)covmatFile_->Get(covTreeName.c_str());
  if(!covmatTree_) 
    cout << "H4CovMat::init: ERROR, unable to load covmat Tree:"
	 << covTreeName << " from file " << covFileName << endl;

  covmatHist_ = (TH1F*)covmatFile_->Get("hRunNumber");
  string branchName;
  if (covmatHist_->GetNbinsX() == 1)
    branchName = "BH4CovMat";
  else {
    int runNumber = H4CovMat::FindNearestRun();
    branchName = (string)Form("BH4CovMat_%d",runNumber);
  }
  covmatBranch_ = (TBranch*)covmatTree_->GetBranch(branchName.c_str());

  // Look for map smCrystal number -> tree entry
  int firstCrystal = -1;
  TH1I* hSMCrystalToIndex = (TH1I*)covmatFile_->Get("hSMCrystalToIndex");
  for (int smCrystal = 1; smCrystal <= hSMCrystalToIndex->GetNbinsX(); 
       smCrystal++) {
    int index = (int)hSMCrystalToIndex->GetBinContent(smCrystal);
    if (index >= 0) {
      smCrystalToIndex_[smCrystal] = index;
      if (firstCrystal == -1) firstCrystal = smCrystal;
    }
  }
  delete hSMCrystalToIndex;

  // One entry / crystal
  geom_.getTowerCrystalNumber((int)curTower1_, (int)curXtal1_, firstCrystal);
  curTower2_ = curTower1_;
  curXtal2_ = curXtal1_;
  H4CovMat::ReadBranch(curTower1_, curXtal1_, curTower2_, curXtal2_);
  */
  return true;
}

//! Switch between mode_ = true (fill the input array) and mode_ = false
void H4CovMat::SetMode(bool mode)
{
  if (mode == mode_) return;
  mode_ = mode;
  if (mode_) {
    delete vsum_;
    vsum_ = new vector<float>(ndim_, 0.);
  } else {
    delete vevt_;
    delete vsum_;
    vevt_ = new vector<float>;
    vsum_ = new vector<float>;
  }
  return;
}

//! Fill vector storing events
void H4CovMat::Fill(size_t sample, float val)
{
  if (mode_ == false) {
    cout << "H4CovMat::Fill: arrays you want to fill do not exist"
	 << " in this mode" << endl;
    return;
  }
  H4CovMat::CheckSize(sample,"H4CovMat::Fill");
  if (sample == 0) nevt_++;
  vevt_->push_back(val);
  if (isStationarity_)
    (*vsum_)[0] += val;
  else
    (*vsum_)[sample] += val;
  return;
}

//! Compute sample mean and covariance matrices inside a given crystal
bool H4CovMat::ComputeEstimators()
{
  // Compute the estimators (mean, covariance) from vevt
  // Only in mode_ = true 
  //
  if (mode_ == false) {
    cout << "H4CovMat::ComputeEstimators: You cannot compute the estimators"
	 << " in this mode" << endl;
    return false;
  }
  if (nevt_ < 2) return false;
  lfrms_ = 0.;
  hfrms_ = 0.;
  if (isStationarity_) {
    // Compute the mean values from vsum
    (*vmean_)[0] = (*vsum_)[0]/float(ndim_*nevt_);
    // Loop on all the events to compute the covariance matrix
    for (unsigned int i = 0; i < ndim_; i++)
      (*vcov_)[i] = 0.;
    for (unsigned int k = 0; k < nevt_; k++) {
      float lfmean = 0.;
      for (unsigned int i = 0; i < ndim_; i++)
	lfmean += (*vevt_)[k*ndim_+i];
      lfmean /= (float)ndim_;
      float mean = (*vmean_)[0];
      lfrms_ += (lfmean - mean)*(lfmean - mean);
      if (isLFsub_) mean = lfmean;
      for (unsigned int i = 0; i < ndim_; i++) {
	for (unsigned int j = i; j < ndim_; j++)
	  (*vcov_)[j - i] += ((*vevt_)[k*ndim_+i] - mean)* 
	                     ((*vevt_)[k*ndim_+j] - mean);
	hfrms_ += ((*vevt_)[k*ndim_+i] - lfmean)*
	          ((*vevt_)[k*ndim_+i] - lfmean);
      }
    }
    for (unsigned int i = 0; i < ndim_; i++)
      (*vcov_)[i] = (*vcov_)[i]/((ndim_ - i)*nevt_ - 1);
  } else {
    // Compute the mean values from vsum
    double sum = 0.;
    for (unsigned int i = 0; i < ndim_; i++) {
      (*vmean_)[i] = (*vsum_)[i]/float(nevt_);
      sum += (*vsum_)[i];
    }
    (*vmean_)[ndim_] = sum/float(ndim_*nevt_); 
    // Loop on all the events to compute the covariance matrix
    sum = 0.;
    for (unsigned int i = 0; i < ndim_; i++)
      for (unsigned int j = 0; j <= i; j++)
	(*vcov_)[i*(i + 1)/2 + j] = 0.;
    for (unsigned int k = 0; k < nevt_; k++) {
      float mean = (*vmean_)[ndim_];
      float lfmean = 0.;
      for (unsigned int i = 0; i < ndim_; i++)
	lfmean += (*vevt_)[k*ndim_+i];
      lfmean /= (float)ndim_;
      lfrms_ += (lfmean - mean)*(lfmean - mean);
      for (unsigned int i = 0; i < ndim_; i++) {
	float meani = (*vmean_)[i];
	for (unsigned int j = 0; j <= i; j++) {
	  float meanj = (*vmean_)[j];
	  if (isLFsub_) // we overwrite mean values by lfmean
	    mean = meani = meanj = lfmean; 
	  (*vcov_)[i*(i + 1)/2 + j] += ((*vevt_)[k*ndim_+i] - meani)* 
	                               ((*vevt_)[k*ndim_+j] - meanj);
	  if (i == j) {
	    sum += ((*vevt_)[k*ndim_+i] - mean)* 
	           ((*vevt_)[k*ndim_+j] - mean);
	    hfrms_ += ((*vevt_)[k*ndim_+i] - lfmean)* 
	              ((*vevt_)[k*ndim_+j] - lfmean);
	  }
	}
      }
    }
    for (unsigned int i = 0; i < ndim_; i++)
      for (unsigned int j = 0; j <= i; j++)
	(*vcov_)[i*(i + 1)/2 + j] /= (nevt_ - 1.);
    (*vcov_)[ndim_*(ndim_ + 1)/2] = sum/(ndim_*nevt_ - 1.);
  }
  lfrms_ /= float(nevt_ - 1.);
  hfrms_ /= float(ndim_*nevt_ - 1.);
  return true;
}

//! Get pedestal mean value
float H4CovMat::getPedestal()
{

  if (isStationarity_) return (*vmean_)[0];
  return (*vmean_)[ndim_];
}

//! Get pedestal sigma
float H4CovMat::getPedestalSigma()
{
  if (isStationarity_) return sqrt((*vcov_)[0]);
  return sqrt((*vcov_)[ndim_*(ndim_ + 1)/2]);
}

//! Get correlation
float H4CovMat::getCorrelation(int sample1, int sample2)
{
  H4CovMat::CheckSize(sample1, "H4CovMat::getCorrelation");
  H4CovMat::CheckSize(sample2, "H4CovMat::getCorrelation");
    if (sample2 > sample1) H4CovMat::Swap(&sample1, &sample2);  
    if (sample1 == sample2) return 1.;
    if (isStationarity_) {
      float rms2 = (*vcov_)[0];
      if (rms2 <= 0.) return 0.;
      return (*vcov_)[sample1 - sample2]/rms2;
    } else {
      float rms2 = (*vcov_)[ndim_*(ndim_ + 1)/2]; // take the average value
      if (rms2 <= 0.) return 0.;
      return (*vcov_)[sample1*(sample1 + 1)/2 + sample2]/rms2;
    }
  return (*vcov_)[sample1*ndim_ + sample2];
}

//! Check if i < ndim_
void H4CovMat::CheckSize(size_t i, string sroutine) const
{
  if (i >= ndim_) {
    cerr << sroutine << ": FATAL: element number i = " << i
	 << " is greater than the array size dim = " << ndim_ << endl;
    abort();
  }
  return;
}

//! Swap i & j
void H4CovMat::Swap(int *i, int *j) const
{
  int k = *i;
  *i = *j;
  *j = k;
  return;
}
