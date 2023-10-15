// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;
// constexpr const bool useBoost = false;
#include <boost/math/special_functions/beta.hpp>
// [[Rcpp::depends(BH)]]
#include "../h/textbookRegularizedIncompleteBeta.h"
#include "../h/regularizedIncompleteBetaFun.h"


#include "../h/processPmfs.h"
#include "../h/trbCharlie.hpp"


#include "../h/charlieThreadPool.hpp"
#include "../h/regrid.h"
#include "../h/upscaleBoundedMean.h"
#include "../fineTuneDiscretizationFulfillQArequirements.hpp"



