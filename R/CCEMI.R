#' Constrained Chained Equations Multiple Imputation (CCEMI)
#'
#' Multiple Imputation of More than One Exposure with Non-differential Measurement Error
#'
#' Constrained Chained Equations Multiple Imputation (CCEMI) treats the true concentrations in the new samples as missing and utilize multiple imputation to fill in these values based on the non-differential measurement error (NDME) assumption.
#'
#' @param data A data frame or a matrix containing the incomplete data. Missing values are coded as NA.
#' @param predM A numeric matrix of length(blocks) rows and ncol(data) columns, containing 0/1 data specifying the set of predictors to be used for each target column. Each row corresponds to a variable block, i.e., a set of variables to be imputed. A value of 1 means that the column variable is used as a predictor for the target block (in the rows). The matrix should base on the calibration design you want.
#' @param nCalib Row number of calibration data
#' @param nMain Row number of Main data
#' @param model The model you want to analysis, example: "y ~ x1 + x2 + x3 + z"
#' @param nImp Imputation number per boot, normally set as 2.
#' @param nBoot Boot number, normally set >=200
#' @param seed An integer that is used as argument by the set.seed() for offsetting the random number generator. Default is to leave the random number generator alone.
#' @return impst: A list of imputed datasets. imp_complete: Complete imputed data. result: Point estimate, ci, variance
#' @author Yuanzhi Yu, Qixuan Chen
#' @references \url{https://github.com/yy3019/CCEMI}
#' @importFrom bootImpute mice AER MCMCglmm
#' @export
#' @examples see \url{https://github.com/yy3019/CCEMI}
#' @name CCEMI
#'
#'




library(bootImpute)
library(mice)

CCEMI = function(data, predM, nCalib, nMain, model, nImp = 2, nBoot = 20, ...){

imp1 = data %>% as.tibble()

impOnce <- function(inputData,M) {
  miceImps <- mice::mice(inputData, m=M, maxit = 20, print = F, method = "norm", predictorMatrix = predM) #, remove.collinear = FALSE
  imps <- vector("list", M)
  for (i in 1:M) {
    imps[[i]] <- mice::complete(miceImps,i)
  }
  imps
}
impst = bootImpute(imp1, impOnce, nBoot=nBoot, nImp = nImp, M = nImp)

imp_complete = complete(mice(imp1, m=20, maxit = 20, print = F, method = "norm", predictorMatrix = predM))

t1 = nCalib + 1
t2 = nCalib + nMain
analyseImp <- function(inputData) {
  inputData2 = inputData[t1:t2,]
  coef(lm(model,data=inputData2))
}
resultt = bootImputeAnalyse(impst, analyseImp, quiet = TRUE)


return(list(impst = impst, imp_complete = imp_complete, result = resultt))

}

