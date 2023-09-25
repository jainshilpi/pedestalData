// H4CovMat.h
//
// Class used to compute & get the mean, variance and correlation between 
// samples and xtals. 
//
// last change : $Date: 2004/11/03 09:57:35 $
// by          : $Author: brunelie $
//

#ifndef H4COVMAT_H
#define H4COVMAT_H

// C++ standard library include files
#include <string>
#include <vector>
#include <map>

class H4CovMat {
 public:

  // Default Constructor, mainly for Root
  H4CovMat();
  ~H4CovMat();
  H4CovMat& operator=(const H4CovMat&); // Assign

  bool init(int run);

  float getPedestal();
  float getPedestalSigma();
  float getCorrelation(int sample1, int sample2);

  void SetMode(bool mode);
  void Fill(size_t sample, float val); // only if mode_ == true

  bool ComputeEstimators(); // only if mode_ == true
  bool GetMode() const { return mode_; }  // Which mode_ are we using 
  bool IsStationarity() const { return isStationarity_; }
  bool IsLFSubtraction() const { return isLFsub_; }

  size_t GetDim() const { return ndim_; } // Return covariance matrix dimension
  size_t GetNevt() const { return nevt_; } // Return the number of events

  float GetLFRMS() const { return sqrt(lfrms_); } // low frequency rms
  float GetHFRMS() const { return sqrt(hfrms_); } // high frequency rms

 private:
  // Private methods
  void CheckSize(size_t i, std::string sroutine) const; // Check if i < ndim_
  void Swap(int *i, int *j) const; // Swap i & j

 private:
  // Class members
  bool           mode_;       // do we store ndim_ values/event into an array ?
  bool           isStationarity_; // do we assume noise stationarity ?
  bool           isLFsub_;    // do we subtract lfmean for each event ?
  size_t         ndim_;       // square matrice & vector dimension
  size_t         nevt_;       // number of events used to compute estimators
  size_t         curTower1_;  // current tower number
  size_t         curXtal1_;   // current crystal number
  size_t         curTower2_;  // current tower number
  size_t         curXtal2_;   // current crystal number
  std::vector<float>* vmean_; // estimated mean value
  std::vector<float>* vcov_;  // estimated sample by sample covariance matrix
  std::vector<float>* vsum_;  // == sum_{k = 0}^{N} x_k
  std::vector<float>* vevt_;  // ndim values / event
  float          lfrms_;      // low frequency rms
  float          hfrms_;      // high frequency rms
  float          lfcorr_;     // low frequency correlation
  std::map<size_t,size_t> smCrystalToIndex_; // SM crystal number -> tree entry
};

#endif
