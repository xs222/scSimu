#' Simulate scRNA count data using zero inflated model
#'
#' This function generates zero-inflated scRNA-seq count data. It first generates a baseline count matrix from a negative binomial distribution. Then, it introduces additional zeros into the matrix by assuming the zero-inflation is a technique artifact distinct from biological zeros.

#' @param mu A vector contains gene mean expression level.
#' @param alpha A vector contains gene dispersion parameter (Shape parameter in the gamma distribution).
#' @param zero_inflated_gene A vector contains the names of zero-inflated genes
#' @param zero_inflated_p A vector contains the zero-inflation parameters. The order of elements in the vector must match the order in the zero_inflated_gene.
#' @param count_dat A gene x cell count matrix that user wants to mimic. If a count matrix is provided, then gene name, cell name and sequencing depth will be calculated automatically.
#' @param gene_name A vector contains gene names.
#' @param cell_name A vector contains cell name.
#' @param seq_depth A vector contains cell sequencing depth.
#' @param IND Simulate independent genes or not.
#' @param cor_mat A correlation matrix for simulating correlated genes. If count matrix is not provided, this is required for simulating correlated genes. For ZINB model, we recommend that users provide their own pre-estimated correlation structure.
#' @param cor_gene User can also provides a vector of gene names if the count matrix is provided. scSimu will automatically estimate the correlation matrix for these genes.
#' @param sig_level The p-value threshold if a vector of gene names is provided to estimate the correlation matrix.
#' @param strength_level The correlation strength threshold if a vector of gene names is provided to estimate the correlation matrix.
#' @param seed Random seed
#' @return A simulated count data
#' @export
#'


scSimu_ZINB <- function(mu, alpha, zero_inflated_gene, zero_inflated_p,
                        count_dat=NULL, gene_name=NULL, cell_name=NULL,
                        seq_depth=NULL, IND=T, cor_mat=NULL, cor_gene=NULL,
                        sig_level=0.05, strength_level=0, seed=7232024){
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
  # introduce zero inflation
  set.seed(seed)
  zero_inflated_matrix <- simu_count
  for (g in rownames(simu_count)) {
    num_cells <- ncol(simu_count)
    if (g %in% zero_inflated_gene){
      rate = zero_inflated_p[zero_inflated_gene==g]
      zero_mask = rbinom(n = num_cells, size = 1, prob = rate)
    } else{
      zero_mask = rep(0, num_cells)
    }

    zero_inflated_matrix[g, ] <- zero_inflated_matrix[g, ] * (1 - zero_mask)
  }
  return(zero_inflated_matrix)
}




