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
