#' Simulate drop-seq scRNA count data
#'
#' This function simulates Drop-seq scRNA count data using negative binomial distribution.
#' The input can be a count matrix with gene mean expression level and dispersion parameter or user defined parameters.
#'
#' @param mu A vector contains gene mean expression level.
#' @param alpha A vector contains gene dispersion parameter (Shape parameter in the gamma distribution).
#' @param count_dat A gene x cell count matrix that user wants to mimic. If a count matrix is provided, then gene name, cell name and sequencing depth will be calculated automatically.
#' @param gene_name A vector contains gene names.
#' @param cell_name A vector contains cell name.
#' @param seq_depth A vector contains cell sequencing depth.
#' @param IND Simulate independent genes or not.
#' @param cor_mat A correlation matrix for simulating correlated genes. If count matrix is not provided, this is required for simulating correlated genes.
#' @param cor_gene User can also provides a vector of gene names if the count matrix is provided. scSimu will automatically estimate the correlation matrix for these genes.
#' @param sig_level The p-value threshold if a vector of gene names is provided to estimate the correlation matrix.
#' @param strength_level The correlation strength threshold if a vector of gene names is provided to estimate the correlation matrix.
#' @param seed Random seed
#' @return A simulated count data
#' @export
#'


scSimu <- function(mu, alpha, count_dat=NULL, gene_name=NULL, cell_name=NULL, seq_depth=NULL, IND=T,
                   cor_mat=NULL, cor_gene=NULL, sig_level=0.05, strength_level=0, seed=7232024){
  if (is.null(count_dat)){
    if (is.null(gene_name) | is.null(cell_name) | is.null(seq_depth)){
      stop("You need to provide a gene x cell count matrix or a set of parameters including gene name, cell name, sequencing depth.")
    } else{
      simu_count <- NB_copula(gene_name, cell_name, seq_depth, mu, alpha,IND, cor_mat, seed)
    }
  } else{
    gene_name <- rownames(count_dat)
    cell_name <- colnames(count_dat)
    seq_depth <- colSums(count_dat)

    if (!IND){
      if (is.null(cor_mat)){
        library(CSCORE)
        library(Seurat)
        sc_obj <- CreateSeuratObject(count_dat)
        if (length(cor_gene)>1){
          cscore_est <- CSCORE(sc_obj, genes = cor_gene, seq_depth = seq_depth)
        } else{
          stop("You need to provide a correlation matrix or specify the correlated gene: Highly expressed, Highly variable, or a vector of gene names.")
        }

        filtered_ori_ests <- cscore_est$est
        filtered_ori_ests[MatrixBH(cscore_est$p_value) >= sig_level] <- 0
        filtered_ori_ests[abs(filtered_ori_ests)<strength_level] <- 0
        cor_mat <- (filtered_ori_ests+t(filtered_ori_ests))/2
      }
    }

    simu_count <- NB_copula(gene_name, cell_name, seq_depth, mu, alpha,IND, cor_mat, seed)


  }
}

MatrixBH <- function(p_matrix){
  p_Matrix_BH = p_matrix - p_matrix
  p_Matrix_BH[upper.tri(p_Matrix_BH)] = p.adjust(p_matrix[upper.tri(p_matrix)], method = "BH")
  return(p_Matrix_BH + t(p_Matrix_BH))
}


NB_copula <- function(gene_name, cell_name, seq_depth, mu, alpha,
                               IND=F, cor_mat=NULL, seed=7232024){

  set.seed(seed)
  ngene <- length(mu)
  ncell <- length(seq_depth)

  # get beta based on mu=alpha*beta
  beta <- mu/alpha

  ###### generate independent expression matrix based on gamma(alpha, beta) ######
  exp_mat <- matrix(rgamma(ngene*ncell, shape = rep(alpha,ncell),
                           scale = rep(beta,ncell)), nrow = ngene)
  colnames(exp_mat) <- cell_name
  rownames(exp_mat) <- gene_name

  ##################### generate correlated gene expression ######################
  if (!IND){
    if (is.null(cor_mat)){
      stop("You need to provide a correlation matrix with gene name for simulating correlated data")
    }
    ncor_gene <- nrow(cor_mat)
    cor_gene_name <- colnames(cor_mat)
    cor_alpha <- alpha[cor_gene_name]
    cor_beta <- beta[cor_gene_name]

    # copula
    copula <- copula_fun(cor_mat, ncor_gene, ncell)
    cor_exp_matrix <- matrix(NA, nrow=ncor_gene, ncol=ncell)
    rownames(cor_exp_matrix) <- cor_gene_name
    for (i in 1:ncor_gene){
      cor_exp_matrix[i,] <- qgamma(copula[i,], shape=cor_alpha[i], scale=cor_beta[i])
    }

    # replace the value for correlated gene in the independent matrix
    exp_mat[cor_gene_name,] <- cor_exp_matrix
  }

  ################ generate count matrix by draw from poisson ####################
  seq_depth_matrix <- matrix(seq_depth, nrow = ngene, ncol = ncell, byrow = T)
  pois_para <- exp_mat*seq_depth_matrix
  count_mat <- matrix(rpois(ngene*ncell, lambda = c(pois_para)), nrow=ngene)
  colnames(count_mat) <- cell_name
  rownames(count_mat) <- gene_name
  return(count_mat)
}



copula_fun <- function(cor_mat, ncor_gene, ncell){
  tryCatch(
    {
      cor_mat_up <- chol(cor_mat)
      copula <- matrix(rnorm(ncor_gene*ncell),nrow = ncor_gene)
      copula <- pnorm(t(cor_mat_up)%*%copula)
      return(copula)

    }, error = function(e) {
      print(e)
      print("The correlation matrix is not positive semi-definite. Use eigen decomposition.")
      mvnrv <- mvtnorm::rmvnorm(ncell, mean = rep(0, dim(cor_mat)[1]), sigma = cor_mat,
                                checkSymmetry = FALSE, method = "eigen")
      copula <- pnorm(t(mvnrv))
      return(copula)
    }
  )
}




