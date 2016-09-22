#' Run EDec stage 1 algorithm.
#'
#' \code{run_edec_stage_1} takes as input the methylation profiles of complex tissue
#' samples, the set of loci with high variability in methylation across cell
#' types, and the number of constituent cell types. It then estimates the
#' average methylation profiles of constituent cell types, and the proportions
#' of constituent cell types in each input sample.
#'
#' The first stage of EDec performs constrained matrix factorization to find
#' cell type specific methylation profiles and constituent cell type proportions
#' that minimize the Euclidian distance between their linear combination and the
#' original matrix of tissue methylation profiles. The minimization algorithm
#' involves an iterative procedure that, in each round, alternates between
#' estimating constituent cell type proportions (using \code{\link{estimate_props_qp}}
#' function) and methylation profiles (using \code{\link{estimate_meth_qp}} function)
#' by solving constrained least squares problems through quadratic programming.
#' The minimization problem is made tractable by the constraints that
#' methylation measurements (beta values) and cell type proportions are numbers
#' in the [0,1] interval, and that cell type proportions within a sample add up
#' to one. These constraints restrict the space of possible solutions, thus
#' making it possible for the local iterative search to reproducibly find a
#' global minimum and an accurate solution. One key requirement for EDec is that
#' cell type proportions vary across samples. A second requirement is that there
#' must be significant differences across constituent cell type methylation
#' profiles. The latter requirement can be met by providing EDec with loci
#' expected to vary in methylation levels across constituent cell types.
#'
#' @param meth_bulk_samples Matrix with methylation profiles of bulk tissue
#'   samples. Rows correspond to loci/probes and columns correspond to different
#'   samples.
#' @param informative_loci A vector containing names (strings) of rows
#'   corresponding to loci/probes that are informative for distinguishing cell
#'   types.
#' @param num_cell_types Number of cell types to use in deconvolution.
#' @param max_its Maximum number of iterations after which the algorithm will
#'   stop.
#' @param rss_diff_stop Maximum difference between the residual sum of squares
#' of the model in two consecutive iterations for the algorithm to converge.
#' @return A list with the following components:
#' @return \describe{
#'  \item{\code{methylation}}{A matrix with average methylation profiles of
#'   constituent cell types. Rows represent different loci/probes and columns
#'   represent different cell types.}
#'  \item{\code{proportions}}{A matrix with proportions of constituent cell
#'   types in each input sample. Rows represent different samples. Columns
#'   represent different cell types.}
#'  \item{\code{iterations}}{Number of iterations the method went through before
#' reaching convergence or maximum number of iterations.}
#'  \item{\code{explained.variance}}{Proportion of variance in input methylation
#'   profiles over informative loci explained by the final model.}
#'  \item{\code{res.sum.squares}}{Residual sum of squares for the final model
#'   over the set of informative loci.}
#'  \item{\code{aic}}{Akaike Information Criterion for the final model over the
#'   set of informative loci.}
#'  \item{\code{rss.per.iteration}}{Vector of residual sum of squares for the
#'   models generated in each iteration of the algorithm.}
#'  }
#'
#' @export
run_edec_stage_1 <- function(meth_bulk_samples,
                             informative_loci,
                             num_cell_types,
                             max_its = 2000,
                             rss_diff_stop = 1e-10) {

  # ---------------------------------------------------------------------------
  # Check for NA values in input methylation profiles
  # ---------------------------------------------------------------------------

  if (sum(is.na(meth_bulk_samples)) > 0) {
    warning("Your input methylation profiles contain NA values.
            Loci with NA values in any samples will not be included
            in the analysis, and will not be present in cell type
            specific methylation profiles.")
  }

  # ---------------------------------------------------------------------------
  # Prepare data for analysis
  # ---------------------------------------------------------------------------

  meth_bulk_samples <- as.matrix(stats::na.omit(meth_bulk_samples))
  informative_loci <- intersect(informative_loci, rownames(meth_bulk_samples))
  meth_bulk_over_inf_loci <- meth_bulk_samples[informative_loci, ]
  num_samples <- ncol(meth_bulk_samples)
  num_loci <- nrow(meth_bulk_over_inf_loci)

  # ---------------------------------------------------------------------------
  # Generate initial estimates of proportions and methylation profiles of
  # constituent cell types over informative loci
  # ---------------------------------------------------------------------------

  estimated_proportions <- gtools::rdirichlet(num_samples,
                                              1 * rep(1/num_cell_types, num_cell_types))

  estimated_cell_type_meth <- estimate_meth_qp(meth_bulk_over_inf_loci,
                                               estimated_proportions)

  # ---------------------------------------------------------------------------
  # Compute residual sum of squares for initial estimates
  # ---------------------------------------------------------------------------

  rss <- norm(meth_bulk_over_inf_loci -
                estimated_cell_type_meth %*% t(estimated_proportions),
              "F")^2

  # ---------------------------------------------------------------------------
  # Go through iterative rounds of estimating proportions and methylation
  # profiles of constituent cell types over informative loci until either
  # convergence or maximum number of iterations are reached
  # ---------------------------------------------------------------------------
  it <- 0
  rss_diff <- Inf
  rss_per_iteration <- c()

  while (it < max_its & rss_diff > rss_diff_stop) {
    rss_per_iteration <- c(rss_per_iteration, rss)
    estimated_proportions <- estimate_props_qp(meth_bulk_over_inf_loci,
                                               estimated_cell_type_meth)
    estimated_cell_type_meth <- estimate_meth_qp(meth_bulk_over_inf_loci,
                                                 estimated_proportions)
    new_rss <- norm(meth_bulk_over_inf_loci -
                      estimated_cell_type_meth %*% t(estimated_proportions),
                    "F")^2
    rss_diff <- rss - new_rss
    rss <- new_rss
    it <- it + 1
  }

  # ---------------------------------------------------------------------------
  # Compute goodness of fit metrics
  # ---------------------------------------------------------------------------

  num_samples <- num_samples * num_loci
  num_estimated_variables <- num_cell_types * (num_samples + num_loci)

  aic <- num_samples * log(rss/num_samples) +
    2 * num_estimated_variables +
    (2 * num_estimated_variables * (num_estimated_variables + 1))/
    (num_samples - num_estimated_variables - 1)

  explained_variance <- 1 - (rss/(norm(meth_bulk_over_inf_loci, "F")^2))

  # ---------------------------------------------------------------------------
  # Estimate cell type specific methylation profiles for all loci given the
  # final estimates of proportions of constituent cell types
  # ---------------------------------------------------------------------------

  full_estimated_cell_type_meth <- estimate_meth_qp(meth_bulk_samples,
                                                    estimated_proportions)

  # ---------------------------------------------------------------------------
  # Add all output variables to a list and return
  # ---------------------------------------------------------------------------

  result <- list(methylation = full_estimated_cell_type_meth,
                 proportions = estimated_proportions,
                 iterations = it,
                 explained.variance = explained_variance,
                 res.sum.squares = rss,
                 aic = aic,
                 rss.per.iteration = rss_per_iteration)
  return(result)
  }

#' Estimate cell type proportions
#'
#' This function will estimate the proportions of constituent cell types in
#' each input sample, given methylation profiles of complex tissue samples and
#' methylation profiles of constituent cell types.
#'
#' EDec assumes that the methylation profiles of complex tissue samples
#' correspond to the linear combination of cell type proportions and
#' methylation profiles of each cell type. Given the methylation profiles of a
#' set of complex tissue samples and the methylation profiles of constituent
#' cell types this function estimates cell type proportions in each sample by
#' solving constrained least squares problems through quadratic programming. The
#' constraints are that the proportions of constituent cell types are numbers in
#' the [0,1] interval and that the proportions of all cell types in each sample
#' sum up to one.
#'
#' @param meth_bulk_samples Matrix of methylation profiles of bulk complex
#'   tissue samples. Columns correspond to different samples and rows correspond
#' to different loci/probes.
#' @param cell_type_specific_meth Matrix of methylation profiles of constituent
#'   cell types. Columns correspond to different cell types and rows correspond
#'   to different loci/probes.
#'
#' @return Matrix with estimated proportions of constituent cell types in each
#'   sample.
estimate_props_qp <- function(meth_bulk_samples, cell_type_specific_meth) {

  num_cell_types <- ncol(cell_type_specific_meth)
  num_samples <- ncol(meth_bulk_samples)

  # ---------------------------------------------------------------------------
  # Specify the constraints for the least squares solution in the format
  # appropriate for quadratic programming (see quadprog::solve.QP)
  # ---------------------------------------------------------------------------

  a_matrix <- cbind(rep(-1, num_cell_types), diag(num_cell_types))
  b_vector <- c(-1, rep(0, num_cell_types))

  # ---------------------------------------------------------------------------
  # Estimate proportions of constituent cell types in each sample by solving
  # constrained least squares problems through quadratic programming
  # ---------------------------------------------------------------------------

  d_matrix <- t(cell_type_specific_meth) %*% cell_type_specific_meth

  estimate_proportions_single_sample <- function(x) {
    d_vector <- t(cell_type_specific_meth) %*% x
    result <- quadprog::solve.QP(Dmat = d_matrix,
                                 dvec = d_vector,
                                 Amat = a_matrix,
                                 bvec = b_vector,
                                 meq = 1)
    result$solution
  }
  estimated_proportions <- t(apply(meth_bulk_samples, 2,
                                   estimate_proportions_single_sample))
  rownames(estimated_proportions) <- colnames(meth_bulk_samples)
  colnames(estimated_proportions) <- colnames(cell_type_specific_meth)

  return(estimated_proportions)
}

#' Estimate cell type specific methylation
#'
#' Given methylation profiles of complex tissue samples and propotions of
#' contituent cell types in each sample, \code{estimate_meth_qp} will estimate
#' average methylation profiles of constituent cell types across the set of
#' input samples.
#'
#' EDec assumes that the methylation profiles of complex tissue samples
#' correspond to the linear combination of cell type proportions and methylation
#' profiles of each cell type. Given the methylation profiles of a set of
#' complex tissue samples and the proportions of constituent cell types in each
#' sample, this function estimates average methylation profiles of constituent
#' cell types by solving constrained least squares problems through quadratic
#' programming. The constraint is that the methylation profiles of constituent
#' cell types are numbers in the [0,1] interval.
#'
#' @param meth_bulk_samples Matrix of methylation profiles of bulk complex
#'   tissue samples. Columns correspond to different samples and rows correspond
#'   to different loci/probes.
#' @param cell_type_props Matrix of proportions of constituent cell types.
#' Columns correspond to different cell types and rows correspond to different
#' bulk tissue samples.
#'
#' @return Matrix with estimated average methylation profiles of constituent
#'   cell types.
estimate_meth_qp <- function(meth_bulk_samples, cell_type_props) {

  num_cell_types <- ncol(cell_type_props)
  num_loci <- nrow(meth_bulk_samples)

  # ---------------------------------------------------------------------------
  # Specify the constraints for the least squares solution in the format
  # appropriate for quadratic programming (see quadprog::solve.QP)
  # ---------------------------------------------------------------------------

  a_matrix <- diag(rep(1, num_cell_types))
  a_matrix <- cbind(a_matrix, diag(rep(-1, num_cell_types)))
  b_vector <- c(rep(0, num_cell_types), rep(-1, num_cell_types))

  # ---------------------------------------------------------------------------
  # Estimate average methylation profiles of constituent cell types by solving
  # constrained least squares problems through quadratic programming
  # ---------------------------------------------------------------------------

  d_matrix <- t(cell_type_props) %*% cell_type_props

  estimate_cell_type_meth_single_locus <- function(x) {
    d_vector <- x %*% cell_type_props
    result <- quadprog::solve.QP(Dmat = d_matrix,
                                 dvec = d_vector,
                                 Amat = a_matrix,
                                 bvec = b_vector,
                                 meq = 0)
    result$solution
  }

  estimated_cell_type_meth <- t(apply(meth_bulk_samples, 1,
                                      estimate_cell_type_meth_single_locus))
  rownames(estimated_cell_type_meth) <- rownames(meth_bulk_samples)
  colnames(estimated_cell_type_meth) <- colnames(cell_type_props)

  return(estimated_cell_type_meth)
}
