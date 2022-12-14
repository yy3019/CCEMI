#' Constrained Iterative Conditional Imputation (CICI)
#'
#' Multiple Imputation of More than One Exposure with Non-differential Measurement Error
#'
#' Constrained Iterative Conditional Imputation (CICI) treats the true concentrations in the new samples as missing and utilize multiple imputation to fill in these values based on the non-differential measurement error (NDME) assumption.
#'
#' @param data A data frame or a matrix containing the incomplete data. Missing values are coded as NA.
#' @param predM A numeric matrix of length(blocks) rows and ncol(data) columns, containing 0/1 data specifying the set of predictors to be used for each target column. Each row corresponds to a variable block, i.e., a set of variables to be imputed. A value of 1 means that the column variable is used as a predictor for the target block (in the rows). The matrix should base on the calibration design you want.
#' @param nCalib Row number of calibration data
#' @param nMain Row number of Main data
#' @param method Can be either a single string, or a vector of strings with length of ncol(data). Example: c("RUTR", "RUTR", "RUTR", "norm", "norm", "norm", "norm", "norm")
#' @param model The model you want to analysis, example: "y ~ x1 + x2 + x3 + z"
#' @param upper_bound The upper bound value of LOD in main (observed) data. example: c(0.05, 0.1, 0.3)
#' @param nImp Imputation number per boot, normally set as 2.
#' @param nBoot Boot number, normally set >=200
#' @param seed An integer that is used as argument by the set.seed() for offsetting the random number generator. Default is to leave the random number generator alone.
#' @return impst: A list of imputed datasets. imp_complete: Complete imputed data. result: Point estimate, ci, variance
#' @author Yuanzhi Yu, Qixuan Chen
#' @references \url{https://github.com/yy3019/CICI}
#' @importFrom bootImpute mice AER MCMCglmm
#' @export
#' @examples see \url{https://github.com/yy3019/CICI}
#' @name CICI_LOD
#'
#'


library(bootImpute)
library(mice)
library(AER)
library(MCMCglmm)

CICI_LOD = function(data, predM, nCalib, nMain, model, method, upper_bound, nImp = 2, nBoot = 20, seed = NA){
imp2 = data %>% as.tibble()
impOnce <- function(inputData,M) {
  miceImps <- mice::mice(inputData, m=M, maxit = 20, print = F, method = method, predictorMatrix = predM, upper_bound = upper_bound, type_c = colnames(imp2)) #, remove.collinear = FALSE
  imps <- vector("list", M)
  for (i in 1:M) {
    imps[[i]] <- mice::complete(miceImps,i)
  }
  imps
}

impst = bootImpute(imp2, impOnce, nBoot=nBoot, nImp = nImp, M = nImp)

imp_complete = complete(mice(imp2, m=20, maxit = 20, print = F, method = method, predictorMatrix = predM, upper_bound = upper_bound, type_c = colnames(imp2)))

t1 = nCalib + 1
t2 = nCalib + nMain

analyseImp <- function(inputData) {
  inputData2 = inputData[t1:t2,]
  coef(lm(model,data=inputData2))
}
resultt = bootImputeAnalyse(impst, analyseImp, quiet = TRUE)


return(list(impst = impst, imp_complete = imp_complete, result = resultt))

}


#' RUTR
#'
#' Calculates imputations for univariate missing data by Tobit regression
#'
#' Although the mice package does not include an option for fitting the tobit regression model, we have derived a custom tobit function that can be integrated with the mice package to allow the imputation of W for those below LODs using tobit imputation models.
#'
#'
#' @author Yuanzhi Yu, Qixuan Chen
#' @seealso \code{\link[tools]{file_ext}}, \code{\link[tools]{file_path_sans_ext}}
#' @references \url{https://github.com/yihui/rmini}
#' @export
#' @examples mice.impute.RUTR(df2)
#' @name mice.impute.RUTR
#'
#'



mice.impute.RUTR = function(y, ry, type, x, wy = NULL, upper_bound, type_c,...){

  ### Choose upper bound for y
  '%ni%' <- Negate('%in%')
  upper_bound_y <- upper_bound[which(type_c %ni% attributes(type)$name)[1]]

  ### Data for tobit
  y1 = y
  y1[wy] = upper_bound_y - 1
  dt = data.frame(y1, x)

  ### Tobit regression, Mu & Sigma
  utobit <- tobit(y1 ~ ., data = dt, left = upper_bound_y, dist = "gaussian")
  mu <- predict(utobit)
  mu_u <- mu[wy]
  sigma_u <- utobit$scale

  ### Impute LOD
  lod = length(which(wy*1 == 1))
  imp_u = c()
  for(i in 1:lod){
    imp_u[i] = rtnorm(1, mean = mu_u[i],sd = sigma_u, upper = upper_bound_y)}

  return(imp_u)
}


