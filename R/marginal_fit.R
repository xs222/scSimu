
#' Estimate the marginal parameters of genes
#'
#' This function estimates the marginal parameters of genes by using the glm_gp()
#' function in glmGamPoi package.
#'
#' @param ori_ct Gene-by-cell count matrix with genes in rows and cells in columns.
#' @return A data.frame containing:
#' \itemize{
#'   \item \code{gene}: Gene name
#'   \item \code{mu}: The mean expression level of gene on log10 scale
#'   \item \code{alpha}: The shape parameter of gene for the Gamma distribution on log10 scale
#'   \item \code{deviances}: Deviance of the fit
#' }
#' @export


marginal_fit <- function(ori_ct){
  library(glmGamPoi)
  size_factors_sel <- colSums(ori_ct)
  gp_ex <- glm_gp(ori_ct, size_factors = size_factors_sel, verbose = T,
                  overdispersion_shrinkage = T, do_cox_reid_adjustment = T)
  marginal_para <- data.frame(gene=rownames(ori_ct),
                              mu = log10(exp(gp_ex$Beta[,1])),
                              alpha = -log10(gp_ex$overdispersions),
                              deviances = gp_ex$deviances)
  return(marginal_para)
}


#' Fit a smooth curve between mu and alpha
#'
#' This function fits a kernel regression smooth curve between mu and alpha.
#'
#' @param marginal_para A data frame where
#' @return A data.frame containing:
#' \itemize{
#'   \item \code{gene}: Gene name
#'   \item \code{mu}: The mean expression level of gene on log10 scale
#'   \item \code{alpha}: The shape parameter of gene for the Gamma distribution on log10 scale
#'   \item \code{deviances}: Deviance of the fit
#' }
#' @export

fit_smooth <- function(marginal_para){
  # fit the kernel regression
  marginal_para$up <- ifelse(marginal_para$alpha>2.5, "upper","lower")
  ## fit line (only use the lower cluster)
  marginal_fit_PNAS_sel <- marginal_para[marginal_para$up=="lower",]
  ## kernel smooth
  km5 <- ksmooth(marginal_fit_PNAS_sel$mu, marginal_fit_PNAS_sel$alpha,
                 kernel="normal", bandwidth = bw.SJ(marginal_fit_PNAS_sel$mu)*5)
}


# library(roxygen2)
# roxygen2::roxygenize()
