#' Select informative loci (EDec stage 0)
#'
#' \code{run_edec_stage_0} selects loci/probes that display cell type specific
#' patterns of methylation, and are likely to be informative for deconvolution.
#' This selection is based on comparison of reference methylation profiles
#' representing different cell types.
#'
#' There are two versions of this method. The user can choose which one to use
#' by setting the \code{version} argument to either "one.vs.rest" (default) or
#' "each.pair". In the first version, the method will perform two-sample T-tests
#' to compare the levels of methylation on each locus/probe between references
#' of one class against all other references. These comparisons will be done for
#' each reference class in turn. The method will then select the most hyper- and
#' hypo-methylated loci/probes that showed p-value less than or equal to
#' \code{max_p_value}. The same number of markers will be selected from each
#' comparison, adding up to \code{num_markers}. The difference between the
#' "one.vs.rest" version and the "each.pair" version of the method is that,
#' instead of comparing one class of references against all others as is done in
#' "one.vs.rest" version, the "each.pair" version will perform comparisons
#' between each pair of reference classes.
#'
#' @param reference_meth A matrix of reference methylation profiles. Rows
#'   correspond to different loci/probes. Columns correspond to different
#'   reference samples.
#' @param reference_classes A vector of numbers or strings representing the classes
#'   associated with each reference methylation profile in
#' \code{reference_meth}.
#' @param max_p_value The maximum p-value, from T-tests comparing reference
#' classes, for a locus to be considered as part of the chosen set of marker
#' loci.
#' @param num_markers Number of marker loci to be extracted.
#' @param version Either "one.vs.rest" or "each.pair". Specifies whether the marker
#'   selection procedure should compare each reference class against all others,
#'   or if independent comparisons should be made for each pair of reference
#'   classes.
#'
#' @return A vector with the names of the loci/probes selected as markers.
#'
#' @examples
#' if (requireNamespace("EDecExampleData",quietly=TRUE)) {
#'   informative_loci <-
#'     run_edec_stage_0(reference_meth = EDecExampleData::reference_meth,
#'                      reference_classes = EDecExampleData::reference_meth_class,
#'                      max_p_value = 1e-5,
#'                      num_markers = 500,
#'                      version = "one.vs.rest")
#' } else {
#'   print("To run this example, please install EDecExampleData package from
#'     github by running devtools::install_github('BRL-BCM/EDecExampleData')")
#' }
#'
#'
#' @export
run_edec_stage_0 <- function(reference_meth,
                             reference_classes,
                             max_p_value,
                             num_markers,
                             version = "one.vs.rest"){

  num_comparisons = NULL
  num_hyper_markers_per_comp = NULL
  t_test_results = NULL

  if (version == "one.vs.rest") {

    num_comps <- length(unique(reference_classes))

    t_test_results <- perform_t_tests_all_classes_one_vs_rest(reference_meth,
                                                              reference_classes)

    num_hyper_markers_per_comp <- ceiling((num_markers / num_comps) / 2)

  } else if (version == "each.pair") {

    t_test_results <- perform_t_tests_all_classes_each_pair(reference_meth,
                                                            reference_classes)

    num_comps <- ncol(t_test_results$p.values)

    num_hyper_markers_per_comp <- ceiling((num_markers / num_comps) / 2)

  } else{
    stop("Error: type has to be either 'one.vs.rest' or 'each.pair'")
  }

  markers <- NULL
  for (i in 1:num_comps) {
    significant_loci <- names(which(t_test_results$p.values[ ,i] <=
                                    max_p_value))
    significant_loci <- setdiff(significant_loci, markers)
    hyper_loci <-
      names(utils::tail(sort(t_test_results$diff.means[significant_loci,i]),
                 num_hyper_markers_per_comp))

    hypo_loci <-
      names(utils::head(sort(t_test_results$diff.means[significant_loci,i]),
                 num_hyper_markers_per_comp))

    markers <- c(markers, hyper_loci, hypo_loci)
  }

  num_extra_markers <- length(markers) - num_markers

  if (num_extra_markers > 0) {
    min_p_values <- apply(t_test_results$p.values[markers, ], 1, min)
    markers_to_drop <- utils::tail(sort(min_p_values, index.return=TRUE)$ix,
                            num_extra_markers)
    markers <- markers[-markers_to_drop]
  }

  return(markers)
}


#' T-tests for comparing all rows of two matrices
#'
#' \code{perform_t_tests_all_rows} performs two sample T-tests comparing the
#' values in each row of one matrix against the values of the matching row in a
#' second matrix.
#'
#' @param data_group_1 A matrix of numbers.
#' @param data_group_2 A matrix of numbers.
#'
#' @return A matrix where each row corresponds to comparisons between the
#'   corresponding rows in \code{data_group_1} and \code{data_group_2}. The
#'   first column contains the p-values, the second column contains the means of
#'   data_group_1 for each row, and the third column contains the means of
#'   data_group_2 for each row.
perform_t_tests_all_rows = function(data_group_1, data_group_2){
  nGroup1 = ncol(data_group_1)
  nGroup2 = ncol(data_group_2)
  dataAll = cbind(data_group_1, data_group_2)
  tTestWithErrorHandling = function(x){
    testResult = try(stats::t.test(x[1:nGroup1],
                            x[(nGroup1 + 1):(nGroup1 + nGroup2)]),
                     silent=TRUE);
    if (is.character(testResult)) {
      warning(testResult)
      c(NA, NA, NA)
    } else {
      c(testResult$p.value, testResult$estimate)
    }
  }
  results = matrix(unlist(apply(dataAll, 1, tTestWithErrorHandling)), ncol=3,
                   byrow=TRUE)
  colnames(results) = c("P.value", "Mean.group.1", "Mean.group.2")
  rownames(results) = rownames(data_group_1)
  return(results)
}

#' Comparison of one class against all others
#'
#' For each row of a matrix, this function performs T-tests comparing the values
#' of that row for samples in one class against the values in that row for all
#' other samples. Such comparisons are done for each reference class in turn.
#'
#' @param data_matrix Matrix containing the data for all samples. Columns
#'   correspond to different samples.
#' @param class_vector A vector of numbers or strings representing the classes
#'   associated with each column in \code{data_matrix}.
#'
#' @return A list containing:
#' @return \describe{
#'   \item{\code{p.values}}{A matrix containing the p-values
#'     for each comparison. Each column corresponds to the comparison between a
#'     class and all others.
#'   }
#'   \item{\code{diff.means}}{A matrix containing the
#'     difference between group means for each comparison. Each column corresponds
#'     to the comparison between a class and all others.
#'   }
#' }
perform_t_tests_all_classes_one_vs_rest = function(data_matrix, class_vector){

  if (ncol(data_matrix) != length(class_vector)) {
    stop("Number of columns of data matrix must be equal to the length of the
         class vector")
  }
  possibleClasses = unique(class_vector)
  nClasses = length(possibleClasses)

  allPvalues = matrix(NA, nrow = nrow(data_matrix), ncol = nClasses)
  allDiffMeans = matrix(NA, nrow = nrow(data_matrix), ncol = nClasses)
  colnames(allPvalues) = possibleClasses
  rownames(allPvalues) = rownames(data_matrix)
  colnames(allDiffMeans) = possibleClasses
  rownames(allDiffMeans) = rownames(data_matrix)

  for (i in 1:nClasses) {
    class = possibleClasses[i]
    resultTest = perform_t_tests_all_rows(data_matrix[ ,class_vector == class],
                                          data_matrix[ ,class_vector != class])
    allPvalues[ ,i] = resultTest[ ,1]
    allDiffMeans[ ,i] = resultTest[ ,2] - resultTest[ ,3]
  }
  result = list(allPvalues, allDiffMeans)
  names(result) = c("p.values", "diff.means")
  return(result)
}

#' Comparison of each pair of classes
#'
#' For each row of a matrix, this function performs T-tests comparing the values
#' of that row for samples in one class against the values in that row samples
#' in another class. Such comparisons are done for each pair of reference
#' classes.
#'
#' @param data_matrix Matrix containing the data for all samples. Columns
#'   correspond to different samples.
#' @param class_vector A vector of numbers or strings representing the classes
#'   associated with each column in \code{data_matrix}.
#'
#' @return A list containing:
#' @return \describe{
#'   \item{\code{p.values}}{A matrix containing the p-values
#'     for each comparison. Each column corresponds to the comparison between a
#'     pair of classes.
#'   }
#'   \item{\code{diff.means}}{A matrix containing the
#'     difference between group means for each comparison. Each column
#'     corresponds to the comparison between a pair of classes.
#'   }
#' }
perform_t_tests_all_classes_each_pair = function(data_matrix, class_vector){
  if (ncol(data_matrix)!=length(class_vector)) {
    stop("Number of columns of data matrix must be equal to the length
         of the class vector")
  }
  possibleClasses = unique(class_vector)
  nClasses = length(possibleClasses)

  allPValues = NULL
  allDiffMeans = NULL
  names = NULL
  for (i in 1:(nClasses - 1)) {
    for (j in (i + 1):nClasses) {
      class1 = possibleClasses[i]
      class2 = possibleClasses[j]
      names = c(names, paste(class1, class2, sep="."))
      result = perform_t_tests_all_rows(data_matrix[ ,class_vector == class1],
                                        data_matrix[ ,class_vector == class2])
      allPValues = cbind(allPValues, result[ ,1])
      allDiffMeans = cbind(allDiffMeans, result[ ,2] / result[ ,3])
    }
  }
  colnames(allPValues) = names
  colnames(allDiffMeans) = names
  result = list(allPValues, allDiffMeans)
  names(result) = c("p.values", "diff.means")
  return(result)
}
