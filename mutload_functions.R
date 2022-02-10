#' Make matrices for depth, variant bases, and reference bases for depth of
#' coverage files
#'
#' doc_matrices creates the matrices that are consecutively used by the mutload
#' function. It will parse depth-of-coverage files to recover total read depth,
#' variant bases, and the sample reference base for all samples at all
#' positions. The functions are made for mutation load estimation, which in this
#' form can only be used on duplex consensus reads.
#'
#' @param files Character vector with file paths to depth of coverage files
#' @param pos Character vector with positions of bases that need to be
#'   extracted. They are generally of the form "chr1:16451706"
#' @param bed Character string with file path to a bed-file with the regions to
#'   be analyzed. Using this argument will override pos
#' @param offset Numeric. Specify how many bases before or after the region
#'   specified by bed should be used. If this argument has length one, it will
#'   only apply to start base. If it has length two, it will apply to the start
#'   and end base respectively.
#' @param names Character. Replaces column names with this.
#'
#' @details Depth of coverage files should be tab-delimited files and are
#'   expected to have a header row and columns for position, A, C, G, T (N is
#'   not used). Positions given in pos should match the position description in
#'   the depth of coverage files; most likely a character vector, each string
#'   containing the chromosome and position separated by a colon. If a bed-file
#'   is used, this is expected to be a tab-delimited file without header and
#'   containing at least three columns with the chromosome, start position and
#'   end position respectively. Offset can be used to pad (or shrink) regions.
#'   For instance, to pad a region with 5 bases, use offset = c(-5, 5).
#'
#' @return Returns a list with three matrices: depth, vb, and base. All matrices
#'   have the base positions in rows and samples in columns.
#'
#' @author Jos B. Poell
#'
#' @export

doc_matrices <- function(files, pos, bed, offset = c(1,0), names) {
  
  if (missing(pos) && missing(bed)) stop("This function requires position information")
  if (!missing(pos) && !missing(bed)) warning("Only using bed-file for position information")
  if (!missing(bed)) {
    if (length(offset) == 1) {offset <- c(offset, 0)}
    requireNamespace("foreach")
    bed <- read.delim(bed, header = F, stringsAsFactors = F)
    pos <- foreach(i=seq_len(nrow(bed)), .combine = c) %do% {
      paste0(bed[i,1], ":", seq(bed[i,2] + offset[1], bed[i,3] + offset[2]))
    }
  }
  if (!missing(names)) {
    if (length(names) != length(files)) {
      warning("names and files of unequal length: using file names")
      names <- files
    }
  } else {names <- files}
 
  depth <- matrix(data = 0, nrow = length(pos),
                  ncol = length(files),
                  dimnames = list(pos, names))
  vb <- depth
  base <- depth
  for (s in seq_along(files)) {
    f <- files[s]
    message(paste0("processing ", f))
    dat <- read.delim(f, stringsAsFactors = FALSE)
    dat <- dat[dat$Locus %in% pos,]
    tempdat <- matrix(data = 0, nrow = length(pos), ncol = 4, dimnames = list(pos, c("A","C","G","T")))
    tempdat[match(dat$Locus, row.names(tempdat)), 1] <- dat[,2]
    tempdat[match(dat$Locus, row.names(tempdat)), 2] <- dat[,3]
    tempdat[match(dat$Locus, row.names(tempdat)), 3] <- dat[,4]
    tempdat[match(dat$Locus, row.names(tempdat)), 4] <- dat[,5]
    depth[,s] <- rowSums(tempdat)
    vb[,s] <- apply(tempdat, 1, function(x) {
      ref <- which.max(x)
      sum(x[-ref])
    })
    base[,s] <- apply(tempdat, 1, which.max)
  }
  # Reference base defaults to 1 (A) when there is no data. To fix:
  base[depth == 0] <- 0
  return(list(depth=depth, vb=vb, base=base))
}


#' Calculate mutation load on a group of samples
#'
#' The premise for mutation load in this function are rare somatic mutations.
#' These variants only have to occur with a single variant-supporting base.
#' However, this base must have been determined with extremely high specificity.
#' Such methods have been described by Schmitt et al. (PNAS 2012) and Hoang et
#' al. (PNAS 2016). Germline variants and clonal variants are largely excluded
#' by having a variant allele frequency significantly higher than the VAF cutoff
#' (default 0.05). Additionally, positions are filtered out that contain variant
#' bases too many samples. The fraction of true signal lost by this filtering
#' step can be set manually (default 0.01). The matrices depth, vb, and base are
#' sample-by-position matrices as provided by doc_matrices.
#'
#' @param depth Integer matrix with total base coverage
#' @param vb Integer matrix with number of variant bases
#' @param base Matrix with integers specifying the most-sequence base for that
#'   particular sample. 1 = A, 2 = C, 3 = G, 4 = T, 0 = no base coverage.
#' @param VAF_cutoff Numeric. Default = 0.05
#' @param signal_cutoff Numeric. Default = 0.01
#' @param noisy Numeric. Default = 0.0008
#' @param FDR Numeric. Default = 0.1
#' @param pos_info Logical. If set to TRUE, a data frame with information on the
#'   base positions is returned. Default = FALSE
#' @param conf.level Logical or numeric. If set, the respective confidence
#'   interval for the mutation load will be reported. If TRUE, the default of
#'   poisson.test is used. Default = FALSE
#'
#' @details In principle, this function counts positions with variant bases, and
#'   divides that number with the total number of bases covered (or sequenced).
#'   Importantly, specific positions are excluded on either a sample or cohort
#'   basis when they are highly unlikely to represent rare somatic mutations.
#'   First of all variants are filtered out when the variant allele frequency is
#'   higher than the VAF_cutoff (default 0.05) with a false discovery rate lower
#'   than FDR (default 0.1). The probability is calculated with the pbinom
#'   function and corrected for multiple testing with p.adjust(method = "BH).
#'   Further power of this function comes from having information from multiple
#'   samples. It is recommended to include all data in a single analysis. This
#'   will give the best informed selection of excluded positions. Positions will
#'   be excluded if the average VAF is higher than noisy (default 0.0008) or if
#'   a higher number of samples have variants in the position than expected by
#'   chance. The cutoff for this is based on what proportion of true signal is
#'   explained by positions with x or more samples with variant bases. The
#'   maximum fraction of signal loss is specified with signal_cutoff (default
#'   0.01). The resulting cut-off in terms of samples is reported in a message.
#'
#' @return Returns a data frame with sample name, (rare somatic) mutations,
#'   bases covered, and mutation load (rare somatic mutations per million bases
#'   covered). If pos_info is set to TRUE, it will return a list with two data
#'   frames: the mutation load data and a data frame with position info:
#'   position, average VAF, samples with variant bases, and whether the row was
#'   included in the analysis.
#'
#' @note Bases covered is the cumulative coverage of all positions that are
#'   included in the analysis (coverage of positions with common SNPs or too
#'   high background error rate is excluded).
#'
#' @author Jos B. Poell
#'
#' @export

mutload <- function(depth, vb, base, VAF_cutoff = 0.05, signal_cutoff = 0.01,
                    noisy = 0.0008, FDR = 0.1, pos_info = FALSE, conf.level = FALSE) {
  pmat <- pbinom(vb-1, depth, VAF_cutoff, lower.tail = F)
  qmat <- apply(pmat, 2, function(x) p.adjust(x, "BH"))
  fvb <- vb
  fvb[qmat < FDR] <- NaN
  fdepth <- depth
  fdepth[qmat < FDR] <- NaN
  fvafrow <- rowSums(fvb, na.rm = T)/rowSums(fdepth, na.rm = T)
  omr <- sum(fvb>0, na.rm = T)/sum(!is.na(fvb))
  n <- ncol(fvb)
  x <- max(which(sapply(seq_len(ceiling(sqrt(n))), function(x) {
    (x+1)*pbinom(x, n, omr, lower.tail = F)/(n*omr)
  }) > signal_cutoff)) + 1
  message(paste0("Positions with ", x+1, " or more samples carrying a variant are excluded"))
  rex <- fvafrow > noisy | is.na(fvafrow) | rowSums(fvb == 0, na.rm = TRUE) < ncol(fvb) - x
  mutations <- apply(fvb[!rex,], 2, function(x) {sum(x > 0, na.rm = T)})
  bases_covered <- apply(fdepth[!rex,], 2, function(x) {sum(x, na.rm = T)})
  sample_info <- data.frame(sample = names(mutations), 
                            mutations, bases_covered, 
                            mutload = 10^6*mutations/bases_covered)
  if (conf.level != FALSE) {
    sample_info$lower = sapply(seq_along(mutations), function(i) {
      poisson.test(mutations[i], bases_covered[i]/1000000)$conf.int[1]})
    sample_info$upper = sapply(seq_along(mutations), function(i) {
      poisson.test(mutations[i], bases_covered[i]/1000000)$conf.int[2]})
  }
  if (pos_info == FALSE) {
    return(sample_info)
  } else {
    pos_info <- data.frame(position = rownames(depth),
                           VAF = fvafrow,
                           nvb = rowSums(fvb > 0, na.rm = TRUE),
                           included = !rex)
    return(list(sample_info=sample_info, pos_info=pos_info))
  }
}
