// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
#define vec std::vector


// constexpr const bool useBoost = false;
#include <boost/math/special_functions/beta.hpp>
// [[Rcpp::depends(BH)]]
#include "h/textbookRegularizedIncompleteBeta.h"
#include "h/regularizedIncompleteBetaFun.h"


#include "h/processPmfs.h"
#include "h/trbCharlie.hpp"


#include "h/distance.hpp"
#include "h/charlieThreadPool.hpp"
#include "solve-d.hpp"
#include "longestSubseq.hpp"
#include "movingAverage.hpp"
#include "autoCorr.hpp"
#include "discretize.hpp"
#include "unbiasedDiscretize.hpp"


#include "h/regrid.h"
#include "h/mergeRegrid.h"
#include "h/mixRegrid.h"
#include "makeEmpDistr.hpp"
#include "mixDistr.hpp"
#include "interpolateEmpDistr.hpp"
#include "assignP0.hpp"
#include "extractMain.hpp"
#include "h/upscaleBoundedMean.h"
#include "correctMeanBias.hpp"
#include "exponentialScaling.hpp"
#include "upscaleAndBoundPMF.hpp"
// #include "fineTuneDiscretizationFulfillQArequirements.hpp"
#include "fineTuneDiscretizationFulfillQArequirements-002.hpp"
#undef vec


#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include "lbfgs/LBFGSBcharlie.hpp"
#include "fitObj.hpp"
#include "distanceAPI.hpp"
#include "ignoranceScore.hpp"









