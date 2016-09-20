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
#' @return \code{methylation} A matrix with average methylation profiles of
#'   constituent cell types. Rows represent different loci/probes and columns
#'   represent different cell types.
#' @return \code{proportions} A matrix with proportions of constituent cell
#'   types in each input sample. Rows represent different samples. Columns
#'   represent different cell types.
run_edec_stage_1 <- function(meth_bulk_samples,
                      informative_loci,
                      num_cell_types,
                      max_its = 2000,
                      rss_diff_stop = 1e-10) {

    if (sum(is.na(meth_bulk_samples)) > 0) {
        warning("Your input methylation profiles contain NA values.
            Loci with NA values in any samples will not be included in the analysis,
            and will not be present in cell type specific methylation profiles.")
    }
    meth_bulk_samples <- as.matrix(na.omit(meth_bulk_samples))
    informative_loci <- intersect(informative_loci, rownames(meth_bulk_samples))
    mixturesOverCTSLoci <- meth_bulk_samples[informative_loci, ]
    nMixtures <- ncol(meth_bulk_samples)
    nLoci <- nrow(mixturesOverCTSLoci)

    estProps <- rdirichlet(nMixtures, 1 * rep(1/num_cell_types, num_cell_types))
    estMeth <- estimate_meth_qp(mixturesOverCTSLoci, estProps)

    rss <- norm(mixturesOverCTSLoci - estMeth %*% t(estProps), "F")^2
    it <- 0
    rssDiff <- Inf
    rssPerIteration <- c()

    while (it < max_its & rssDiff > rss_diff_stop) {
        rssPerIteration <- c(rssPerIteration, rss)
        estProps <- estimate_props_qp(mixturesOverCTSLoci, estMeth)
        estMeth <- estimate_meth_qp(mixturesOverCTSLoci, estProps)
        newRss <- norm(mixturesOverCTSLoci - estMeth %*% t(estProps), "F")^2
        rssDiff <- rss - newRss
        rss <- newRss
        it <- it + 1
    }

    nSamples <- nMixtures * nLoci
    nParam <- num_cell_types * (nMixtures + nLoci)
    aic <- nSamples * log(rss/nSamples) + 2 * nParam + (2 * nParam * (nParam + 1))/(nSamples - nParam -
        1)
    expVar <- 1 - (rss/(norm(mixturesOverCTSLoci, "F")^2))

    fullEstMeth <- estimate_meth_qp(meth_bulk_samples, estProps)

    result <- list(fullEstMeth, estProps, it, expVar, rss, aic, rssPerIteration)
    names(result) <- c("methylation", "proportions", "iterations", "explained.variance", "rss", "aic",
        "rss.per.iteration")
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
#'   to different loci/probes.
#' @param cell_type_specific_meth Matrix of methylation profiles of constituent
#'   cell types. Columns correspond to different cell types and rows correspond
#'   to different loci/probes.
estimate_props_qp <- function(meth_bulk_samples, cell_type_specific_meth) {
    num_cell_types <- ncol(cell_type_specific_meth)
    nMixtures <- ncol(meth_bulk_samples)

    estimate <- matrix(0, nMixtures, num_cell_types)
    rownames(estimate) <- colnames(meth_bulk_samples)
    colnames(estimate) <- colnames(cell_type_specific_meth)

    Amat <- cbind(rep(-1, num_cell_types), diag(num_cell_types))
    b0vec <- c(-1, rep(0, num_cell_types))
    Dmat <- t(cell_type_specific_meth) %*% cell_type_specific_meth

    deconv <- function(x) {
        dvec <- t(cell_type_specific_meth) %*% x
        result <- solve.QP(Dmat, dvec, Amat, b0vec, meq = 1)
        result$solution
    }
    estimate <- t(apply(meth_bulk_samples, 2, deconv))
    return(estimate)
}

#' Estimate methylation profiles of constituent cell types
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
#'   Columns correspond to different cell types and rows correspond to different
#'   bulk tissue samples.
estimate_meth_qp <- function(meth_bulk_samples, cell_type_props) {
    num_cell_types <- ncol(cell_type_props)
    nLoci <- nrow(meth_bulk_samples)

    estimate <- matrix(0, nrow = nrow(meth_bulk_samples), ncol = ncol(cell_type_props))
    rownames(estimate) <- rownames(meth_bulk_samples)
    colnames(estimate) <- colnames(cell_type_props)

    Amat <- diag(rep(1, num_cell_types))
    Amat <- cbind(Amat, diag(rep(-1, num_cell_types)))
    bvec <- c(rep(0, num_cell_types), rep(-1, num_cell_types))
    Dmat <- t(cell_type_props) %*% cell_type_props
    deconv <- function(x) {
        dvec <- x %*% cell_type_props
        result <- solve.QP(Dmat, dvec, Amat, bvec = bvec, meq = 0)
        result$solution
    }
    estimate <- t(apply(meth_bulk_samples, 1, deconv))
    return(estimate)
}


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
run_edec_stage_2 <- function(gene_exp_bulk_samples, cell_type_props) {
    gene_exp_bulk_samples <- as.matrix(gene_exp_bulk_samples)

    num_cell_types <- ncol(cell_type_props)
    nGenes <- nrow(gene_exp_bulk_samples)
    nSamples <- nrow(cell_type_props)

    estimate <- matrix(0, nrow = nrow(gene_exp_bulk_samples), ncol = ncol(cell_type_props))
    rownames(estimate) <- rownames(gene_exp_bulk_samples)
    colnames(estimate) <- colnames(cell_type_props)

    Amat <- diag(rep(1, num_cell_types))
    bvec <- rep(0, num_cell_types)

    for (gene in 1:nGenes) {
        Dmat <- t(cell_type_props) %*% cell_type_props
        dvec <- gene_exp_bulk_samples[gene, ] %*% cell_type_props
        result <- solve.QP(Dmat, dvec, Amat, bvec = bvec, meq = 0)
        estimate[gene, ] <- result$solution
    }
    residuals <- gene_exp_bulk_samples - estimate %*% t(cell_type_props)
    meanSquaredResiduals <- apply(residuals^2, 1, sum)/(nSamples - num_cell_types)
    explainedVariances <- 1 - (apply(residuals^2, 1, sum)/apply((gene_exp_bulk_samples^2), 1, sum))
    m <- solve(t(cell_type_props) %*% cell_type_props)
    vars <- matrix(0, nGenes, num_cell_types)
    for (i in 1:nGenes) {
        vars[i, ] <- meanSquaredResiduals[i] * diag(m)
    }
    sds <- sqrt(vars)
    rownames(sds) <- rownames(estimate)
    colnames(sds) <- colnames(estimate)
    result <- list(estimate, sds, nSamples - num_cell_types, explainedVariances, residuals)
    names(result) <- c("Means", "Std.Errors", "Degrees.of.freedom", "Explained.Variances", "Residuals")
    return(result)
}
