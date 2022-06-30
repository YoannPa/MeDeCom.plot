
#' Plots MeDeComSet CV and RMS error following kappa and lambda values tested
#'
#' @param MDCset     A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.
#' @param err.type   A \code{character} specifying the type of error you want to
#'                   explore:
#'                   \itemize{
#'                    \item{'cve' - to explore cross-validation error.}
#'                    \item{'rmse' - to explore root mean square error.}
#'                   }
#' @param plot.title A \code{character} to specify another title than the
#'                   default one for your plot.
#' @return A \code{trellis} 3D lattice graph of the error following kappa and
#'         lambda parameters values.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Load the example data sets
#' data(example.dataset, package="MeDeCom")
#' # Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
#' medecom.result <- runMeDeCom(
#'   D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
#' # Plot default cross-validation error
#' errorplot3d(MDCset = medecom.result)
#' # Plot RMS error
#' errorplot3d(MDCset = medecom.result, err.type = "rmse")
#' @references Lutsik P. et al., MeDeCom: discovery and quantification of latent
#'             components of heterogeneous methylomes. Genome Biol. 2017 Mar
#'             24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
#'             PMCID: PMC5366155.
#' @references [Package ‘lattice’](
#'             https://cran.r-project.org/web/packages/lattice/lattice.pdf)
#' @references [3D Graph - Lattice Package](
#'             http://myrcodes.blogspot.com/2015/10/3d-graph-lattice-package.html)
#' @keywords internal

errorplot3d <- function(
    MDCset, err.type = "cve", plot.title = NULL){
  # Extract and format results
  dt_err <- data.table::melt.data.table(data = data.table::as.data.table(
    MDCset@outputs$`1`[[err.type]], keep.rownames = "kappa"), id.vars = "kappa")
  dt_err[, kappa := gsub(pattern = "K_", replacement = "", x = kappa)]
  dt_err[, variable := gsub(
    pattern = "lambda_", replacement = "", x = variable)]
  dt_err[, kappa := as.factor(kappa)]
  dt_err[, kappa := factor(x = kappa, levels = sort(as.integer(levels(kappa))))]
  dt_err[, variable := as.factor(variable)]
  dt_err[, variable := factor(
    x = variable, levels = sort(as.numeric(levels(variable))))]
  dt_err[, position := seq(.N), by = variable]

  #Get valleys for each lambda distributions
  dt_valleys <- dt_err[
    order(variable)][, quantmod::findValleys(value), by = variable]
  if(nrow(dt_valleys) > 0){
    #Get error values for each valleys
    minimas <- unlist(lapply(X = seq(nrow(dt_valleys)), FUN = function(i){
      dt_err[
        variable == dt_valleys[i]$variable &
          position == dt_valleys[i]$V1-1]$value
    }))
    #Get kappa values for each valleys
    kappa.min <- unlist(lapply(X = seq(nrow(dt_valleys)), FUN = function(i){
      dt_err[
        variable == dt_valleys[i]$variable &
          position == dt_valleys[i]$V1-1]$kappa
    }))
    #Add kappas & error values to dt_valleys
    dt_valleys[, c("Kappa", "valleys") := .(kappa.min, minimas)]
    dt_valleys <- dt_valleys[, -c("V1"), ]
    #Get the minimum among minimas for each kappa value
    dt_minimas <- dt_valleys[, min(valleys), by = Kappa]
    colnames(dt_minimas) <- c("Kappa", "Error")
    #Get final minimas table with lambda and kappa for all minimum errors
    dt_minimas <- data.table::rbindlist(lapply(X = seq(nrow(dt_minimas)), function(i){
      dt_valleys[Kappa == dt_minimas[i]$Kappa & valleys == dt_minimas[i]$Error]
    }))
    colnames(dt_minimas)[c(1,3)] <- c("Lambda", "Error")
    #Create legend keys for local minimas
    dt_minimas[, legend_keys := paste0(
      "L=", Lambda, " K=", Kappa, " Err=", round(Error, digits = 2))]
  } else { dt_minimas <- data.table() }

  if(err.type %in% c("cve", "rmse")){
    colnames(dt_err)[1:3] <- c("Kappa", "Lambda", "Error")
  } else { stop("Unknown error type. Supported: c('cve', 'rmse')") }

  if(err.type == "cve"){
    z_label <- "C.V.\nError"
    if(is.null(plot.title)){
      plot.title <- paste(
        "Cross-validation error evolution following MeDeCom",
        "Kappa and Lambda parameters values", sep = "\n")
    }
  } else {
    z_label <- "RMSE"
    if(is.null(plot.title)){
      plot.title <- paste(
        "Root mean square error evolution following MeDeCom",
        "Kappa and Lambda parameters values", sep = "\n")
    }
  }
  #Lattice 3D plot template for mixing surfaces and dots
  mypanel <- function(x, y, z, x2, y2, z2, ...){
    lattice::panel.wireframe(x, y, z, ...)
    lattice::panel.cloud(x2, y2, z2, ...)
  }
  #Prepare legend for local minimas
  if(nrow(dt_valleys) > 0){
    minimas_legend <- list(
      x = 0.03, y = 0.25, title = "Local minimas", cex.title = 1.2,
      points = list(col = rep("black", nrow(dt_minimas)), cex = 1.4, pch = 19),
      text = list(lab = dt_minimas$legend_keys, cex = 0.95))
  } else { minimas_legend <- NULL }
  # Error plot
  cat("TEST")
  err.plt <- lattice::wireframe(
    x = Error ~ Lambda * Kappa, data = dt_err, drape = TRUE, colorkey = TRUE,
    screen = list(z = -135, x = -80, y = 0), at = sort(dt_err$Error),
    col.regions = colorRampPalette(c("green", "red"))(256),
    zlab = z_label,
    scales = list(
      arrows = FALSE, x = list(labels = levels(dt_err$Lambda)),
      y = list(labels = levels(dt_err$Kappa)),
      z = list (labels = round(seq(
        from = min(dt_err$Error), to = max(dt_err$Error), length.out = 10)))
    ),
    main = plot.title, panel = mypanel,
    x2 = dt_minimas$Lambda, y2 = dt_minimas$Kappa, z2 = dt_minimas$Error,
    pch = 19, alpha = 1, cex = 1.4, col = "black", key = minimas_legend)
  return(err.plt)
}

#' Plots MeDeComSet C.V. error following kappa and lambda values tested
#'
#' @param MDCset     A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.
#' @param plot.title A \code{character} to specify another title than the
#'                   default one for your plot.
#' @return A \code{trellis} 3D lattice graph of the cross-validation error
#'         following kappa and lambda parameters values.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Load the example data sets
#' data(example.dataset, package="MeDeCom")
#' # Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
#' medecom.result <- runMeDeCom(
#'   D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
#' # Plot cross-validation error
#' cve_plot3d(MDCset = medecom.result)
#' @references Lutsik P. et al., MeDeCom: discovery and quantification of latent
#'             components of heterogeneous methylomes. Genome Biol. 2017 Mar
#'             24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
#'             PMCID: PMC5366155.
#' @references [Package ‘lattice’](
#'             https://cran.r-project.org/web/packages/lattice/lattice.pdf)
#' @references [3D Graph - Lattice Package](
#'             http://myrcodes.blogspot.com/2015/10/3d-graph-lattice-package.html)

cve_plot3d <- function(MDCset, plot.title = NULL){
  errplt <- errorplot3d(MDCset = MDCset, plot.title = plot.title)
  return(errplt)
}

#' Plots MeDeComSet R.M.S. error following kappa and lambda values tested
#'
#' @param MDCset     A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.
#' @param plot.title A \code{character} to specify another title than the
#'                   default one for your plot.
#' @return A \code{trellis} 3D lattice graph of the cross-validation error
#'         following kappa and lambda parameters values.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Load the example data sets
#' data(example.dataset, package="MeDeCom")
#' # Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
#' medecom.result <- runMeDeCom(
#'   D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
#' # Plot root mean square error
#' rmse_plot3d(MDCset = medecom.result)
#' @references Lutsik P. et al., MeDeCom: discovery and quantification of latent
#'             components of heterogeneous methylomes. Genome Biol. 2017 Mar
#'             24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
#'             PMCID: PMC5366155.
#' @references [Package ‘lattice’](
#'             https://cran.r-project.org/web/packages/lattice/lattice.pdf)
#' @references [3D Graph - Lattice Package](
#'             http://myrcodes.blogspot.com/2015/10/3d-graph-lattice-package.html)

rmse_plot3d <- function(MDCset, plot.title = NULL){
  errplt <- errorplot3d(
    MDCset = MDCset, err.type = "rmse", plot.title = plot.title)
  return(errplt)
}

#' Plots a MDS of LMCs and samples.
#'
#' @param MDCset     A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.
#' @param D          A methylation \code{matrix} you ran MeDeCom on.
#' @param k          Am \code{integer} to specify the value of the kappa
#'                   parameter.
#' @param lambda     A \code{numeric} to specify the value of the lambda
#'                   parameter.
#' @param cor.method A \code{character} to specify the method to use for
#'                   computing correlation between LMCs and samples
#'                   (Default: cor.method = 'pearson'; For a list of all
#'                   supported correlation methods see \link[stats]{cor}).
#' @param plot.title A \code{character} to specify another title than the
#'                   default one for your plot.
#' @return A \code{gg} plot of the MDS on LMCs and samples.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Load the example data sets
#' data(example.dataset, package="MeDeCom")
#' # Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
#' medecom.result <- runMeDeCom(
#'   D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
#' # Plot Pearson-based MDS
#' mds_plot(MDCset = medecom.result, D = D, k = 5, lambda = 10^-2)
#' # Plot Spearman-based MDS
#' mds_plot(
#'   MDCset = medecom.result, D = D, k = 5, lambda = 10^-2,
#'   cor.method = "spearman")
#' @references Lutsik P. et al., MeDeCom: discovery and quantification of latent
#'             components of heterogeneous methylomes. Genome Biol. 2017 Mar
#'             24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
#'             PMCID: PMC5366155.

mds_plot <- function(
    MDCset, D, k, lambda, cor.method = "pearson", plot.title = NULL){
  #Get LMCs
  That <- MeDeCom::getLMCs(MDCset, k, lambda, 1)
  if(is.null(colnames(That))) { colnames(That) <- paste("LMC", 1:ncol(That)) }
  mdd <- cbind(That, D)
  #Compute correlation and distance
  d <- as.dist(1 - cor(mdd, method = cor.method))
  fit <- cmdscale(d, k = 2, eig = FALSE)
  #Format
  dt_mds <- as.data.table(fit, keep.rownames = "IDs")
  dt_mds[grepl(pattern = "^LMC\\s\\d+$", x = dt_mds$IDs), data.type := "LMCs"]
  dt_mds[is.na(data.type), data.type := "Samples"]
  dt_mds[data.type == "LMCs", label := IDs]
  dt_mds[data.type == "LMCs", text := ""]
  dt_mds[data.type == "Samples", text := IDs]
  dt_mds[data.type == "Samples", label := ""]
  #Make title
  if(is.null(plot.title)){
    substr(x = cor.method, 1, 1) <- toupper(substr(x = cor.method, 1, 1))
    plot.title <- paste0(
      cor.method, "-based MDS of Latent Methylation Components (LMCs) and",
      " samples\n(kappa = ", k, " ; lambda = ", lambda, ")")
  }
  #Plot MDS
  mds.plot <- ggplot2::ggplot(data = dt_mds, mapping = ggplot2::aes(
    x = V1, y = V2, color = data.type, size = data.type)) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(
      data = dt_mds, mapping = ggplot2::aes(label = text), size = 3.2) +
    ggrepel::geom_label_repel(
      data = dt_mds, mapping = ggplot2::aes(label = label), size = 4) +
    ggplot2::scale_color_manual(values = c("red", "black")) +
    ggplot2::scale_size_manual(values = c(2.5, 1.5)) +
    ggplot2::labs(x = "Dimension 1", y = "Dimension 2", color = "Data type") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 13), axis.text = ggplot2::element_text(size = 12),
      legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      panel.grid = ggplot2::element_line(color = NULL),
      plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggtitle(plot.title)
  return(mds.plot)
}

#' Plots a heatmap of LMCs proportions in all methylation samples.
#'
#' @param MDCset      A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.
#' @param D           A methylation \code{matrix} you ran MeDeCom on.
#' @param k           Am \code{integer} to specify the value of the kappa
#'                    parameter.
#' @param lambda      A \code{numeric} to specify the value of the lambda
#'                    parameter.
#' @param dist.method A \code{character} to specify the method to use for
#'                    computing distances between LMCs proportions in
#'                    methylation samples(Default: dist.method = 'manhattan';
#'                    For all supported distance methods see
#'                    \link[parallelDist]{parDist}).
#' @param margins     A ggplot2 \code{margin} object specifying widths of plot
#'                    margins
#'                    (Default: margins = ggplot2::margin(0, 1.2, 0, 0, "cm")).
#' @param plot.title  A \code{character} to specify another title than the
#'                    default one for your plot.
#' @return A \code{type} object returned description.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Load the example data sets
#' data(example.dataset, package="MeDeCom")
#' # Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
#' medecom.result <- runMeDeCom(
#'   D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
#' #Plot LMC proportions heatmap
#' LMCprop_heatmap(MDCset = medecom.result, D = D, k = 5, lambda = 10^-2)
#' @references Pageaud Y. et al., BiocompR - Advanced visualizations for data
#'             comparison.
#' @references Lutsik P. et al., MeDeCom: discovery and quantification of latent
#'             components of heterogeneous methylomes. Genome Biol. 2017 Mar
#'             24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
#'             PMCID: PMC5366155.

LMCprop_heatmap <- function(
    MDCset, D, k, lambda, dist.method = "manhattan",
    margins = ggplot2::margin(0, 1.2, 0, 0, "cm"), plot.title = NULL){
  #Get LMCs proportions
  proportions_mat <- MeDeCom::getProportions(MDCset, K = k, lambda = lambda)
  #Get sample IDs
  colnames(proportions_mat) <- colnames(D)
  #Set title
  if(is.null(plot.title)){
    plot.title <- "Heatmap of LMC proportions in methylation samples"
  }
  #Plot heatmap
  prop_heatmap <- BiocompR::gg2heatmap(
    m = proportions_mat, show.annot = FALSE, y.axis.right = TRUE,
    dend.size = c(0, 4), dist.method = dist.method,
    theme_heatmap = ggplot2::theme(
      axis.text.y.right = ggplot2::element_text(size = 12, color = "black"),
      axis.title.x = ggplot2::element_text(size = 14),
      panel.border = ggplot2::element_rect(
        color = "black", fill = "transparent"),
      plot.margin = margins),
    scale_fill_grad = ggplot2::scale_fill_gradient(
      low = "white", high = "black", limits = c(0, 1)),
    guide_custom_bar = ggplot2::guide_colorbar(
      title = "Proportions", title.vjust = 0.86, ticks.linewidth = 3,
      barwidth = 15, ticks.colour = "red", frame.colour = "black"),
    plot.labs = ggplot2::labs(title = plot.title))
  return(prop_heatmap)
}
