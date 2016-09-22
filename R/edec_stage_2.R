#' Run EDec stage 2 algorithm
#'
#' This function implements the second stage of the EDec method. It takes as
#' input the gene expression profiles of complex tissue samples, and the
#' proportions of constituent cell types in each sample. It then estimates
#' average and standard errors of cell type specific gene expression profiles.
#'
#' EDec assumes that the gene expression profiles of complex tissue samples
#' correspond to the linear combination of cell type proportions and gene
#' expression profiles of each cell type. Given the gene expression profiles of
#' a set of complex tissue samples and the proportions of constituent cell types
#' in each sample, this function estimates average gene expression profiles of
#' constituent cell types by solving constrained least squares problems through
#' quadratic programming. The constraint is that the gene expression profiles of
#' constituent cell types are numbers greater than or equal to zero.
#'
#' @param gene_exp_bulk_samples Matrix of methylation profiles of bulk complex
#'   tissue samples. Columns correspond to different samples and rows correspond
#'   to different loci/probes.
#' @param cell_type_props Matrix of proportions of constituent cell types.
#'   Columns correspond to different cell types and rows correspond to different
#'   bulk tissue samples.
#'
#' @return A list with the following components:
#' @return \describe{
#'    \item{\code{means}}{A matrix with the estimated average gene expression
#'    profiles of constituent cell types. Rows correspond to different genes.
#'    Columns correspond to different cell types.} \item{\code{std.errors}}{A
#'    matrix with estimated standard errors for each cell type specific gene
#'    expression estimate. Rows correspond to different genes. Columns
#'    correspond to different cell types.}
#'    \item{\code{degrees.of.freedom}}{Number of degrees of freedom for
#'    estimates of cell type specific gene expression}
#'    \item{\code{explained.variances}}{Vector with the proportion of variance
#'    in input expression of each gene across samples explained by the final
#'    model.} \item{\code{residuals}}{Matrix with the difference between the
#'    original gene expression values and the linear combination between
#'    proportions of constituent cell types and gene expression profiles of
#'    constituent cell types.}
#'  }
#'
#' @export
run_edec_stage_2 <- function(gene_exp_bulk_samples, cell_type_props) {

  # ---------------------------------------------------------------------------
  # Check for NA values in input methylation profiles
  # ---------------------------------------------------------------------------

  if (sum(is.na(gene_exp_bulk_samples)) > 0) {
    warning("Your input expression profiles contain NA values.
            Loci with NA values in any samples will not be included
            in the analysis, and will not be present in cell type
            specific methylation profiles.")
  }
  gene_exp_bulk_samples <- as.matrix(stats::na.omit(gene_exp_bulk_samples))

  num_cell_types <- ncol(cell_type_props)
  num_genes <- nrow(gene_exp_bulk_samples)
  num_samples <- nrow(cell_type_props)

  # ---------------------------------------------------------------------------
  # Specify the constraints for the least squares solution in the format
  # appropriate for quadratic programming (see quadprog::solve.QP)
  # ---------------------------------------------------------------------------

  a_matrix <- diag(rep(1, num_cell_types))
  b_vector <- rep(0, num_cell_types)

  # ---------------------------------------------------------------------------
  # Estimate average expression profiles of constituent cell types by solving
  # constrained least squares problems through quadratic programming
  # ---------------------------------------------------------------------------

  d_matrix <- t(cell_type_props) %*% cell_type_props

    estimate_cell_type_exp_single_gene <- function(x) {
      d_vector <- x %*% cell_type_props
      result <- quadprog::solve.QP(Dmat = d_matrix,
                                   dvec = d_vector,
                                   Amat = a_matrix,
                                   bvec = b_vector,
                                   meq = 0)
      result$solution
    }

  estimated_cell_type_gene_exp <- t(apply(gene_exp_bulk_samples, 1,
                                      estimate_cell_type_exp_single_gene))

  rownames(estimated_cell_type_gene_exp) <- rownames(gene_exp_bulk_samples)
  colnames(estimated_cell_type_gene_exp) <- colnames(cell_type_props)

  # ---------------------------------------------------------------------------
  # Compute goodness of fit metrics
  # ---------------------------------------------------------------------------

  residuals <- gene_exp_bulk_samples -
    estimated_cell_type_gene_exp %*% t(cell_type_props)

  mean_squared_residuals <- apply(residuals^2, 1, sum)/
    (num_samples - num_cell_types)

  explained_variances <- 1 - (apply(residuals^2, 1, sum)/
                                apply((gene_exp_bulk_samples^2), 1, sum))

  # ---------------------------------------------------------------------------
  # Estimate standard errors for each estimate
  # ---------------------------------------------------------------------------

  m <- solve(t(cell_type_props) %*% cell_type_props)
  diag_m = diag(m)

  compute_variances <- function(x){
    variance = x * diag_m
    return(variance)
  }
  variances = t(sapply(mean_squared_residuals,compute_variances))

  std_devs <- sqrt(variances)
  rownames(std_devs) <- rownames(estimated_cell_type_gene_exp)
  colnames(std_devs) <- colnames(estimated_cell_type_gene_exp)

  # ---------------------------------------------------------------------------
  # Add all output variables to a list and return
  # ---------------------------------------------------------------------------

  result <- list(means = estimated_cell_type_gene_exp,
                 std.errors = std_devs,
                 degrees.of.freedom = num_samples - num_cell_types,
                 explained.variances = explained_variances,
                 residuals = residuals)
  return(result)
}
