###############################################################################################
## TRENDS IN FUNCTIONAL COMPOSITION OF SMALL MAMMAL COMMUNITIES ACROSS MILLENIAL TIME SCALES ##
## Collin S. VanBuren and Marta A. Jarzyna ##
## Ecography 2022 ##
## code by Collin S. VanBuren ##
###############################################################################################

#### Read in packages ####

require(FD) #v. 1.0-12
require(stringr) #v. 1.4.0
require(inla) #v. 20.03.17
require(tidyverse) #v. 1.2.1
require(Bchron) #v. 4.3.0
require(neotoma) #v. 1.7.4
require(rgeos) #v. 0.5-2

#### Modified dbFD and fdisp functions to allow for species not present in the community to be included in the trait matrix  ####
dbFD2 <- function (x, a, w, w.abun = TRUE, stand.x = TRUE, ord = c("podani", 
                                                                   "metric"), asym.bin = NULL, corr = c("sqrt", "cailliez", 
                                                                                                        "lingoes", "none"), calc.FRic = TRUE, m = "max", stand.FRic = FALSE, 
                   scale.RaoQ = FALSE, calc.FGR = FALSE, clust.type = "ward", 
                   km.inf.gr = 2, km.sup.gr = nrow(x) - 1, km.iter = 100, km.crit = c("calinski", 
                                                                                      "ssi"), calc.CWM = TRUE, CWM.type = c("dom", "all"), 
                   calc.FDiv = TRUE, dist.bin = 2, print.pco = FALSE, messages = TRUE) 
{
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  CWM.type <- match.arg(CWM.type)
  km.crit <- match.arg(km.crit)
  if (!is.logical(messages)) 
    stop("'messages' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.FRic)) 
    stop("'stand.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(stand.x)) 
    stop("'stand.x' must be TRUE or FALSE.", "\n")
  if (!is.logical(w.abun)) 
    stop("'w.abun' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FRic)) 
    stop("'calc.FRic' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FDiv)) 
    stop("'calc.FDiv' must be TRUE or FALSE.", "\n")
  if (!is.logical(calc.FGR)) 
    stop("'calc.FGR' musts be TRUE or FALSE.", "\n")
  if (!is.logical(calc.CWM)) 
    stop("'calc.CWM' must be TRUE or FALSE.", "\n")
  if (!is.logical(scale.RaoQ)) 
    stop("'scale.RaoQ' must be TRUE or FALSE.", "\n")
  if (!is.logical(print.pco)) 
    stop("'print.pco' must be TRUE or FALSE.", "\n")
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- dim(x)[1]
    t.x <- dim(x)[2]
    if (is.null(row.names(x))) 
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x))) 
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels"))) 
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (missing(a)) {
    ab.names <- list("Community1", x.rn)
    a <- matrix(1, 1, s.x, dimnames = ab.names)
  }
  else {
    if (is.matrix(a) | is.data.frame(a)) {
      s.a <- dim(a)[2]
      ab.t <- t(a)
      if (is.null(row.names(ab.t))) 
        stop("'a' must have column names.", "\n")
      else ab.t.row <- row.names(ab.t)
      a <- as.matrix(a)
    }
    if (is.vector(a)) {
      s.a <- length(a)
      if (is.null(names(a))) 
        stop("'a' must have names.", "\n")
      else ab.t.row <- names(a)
      ab.names <- list("Community1", ab.t.row)
      a <- matrix(a, 1, s.a, dimnames = ab.names)
    }
    if (s.x != s.a) 
      stop("Different number of species in 'x' and 'a'.", 
           "\n")
    if (any(ab.t.row != x.rn)) 
      stop("Species labels in 'x' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
           "\n")
  }
  a <- as.matrix(a)
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  
  abun.sum2 <- apply(a, 2, sum)
  
  if (!missing(w) & is.dist.x) 
    stop("When 'x' is a distance matrix, 'w' should be left missing.", 
         "\n")
  if (!missing(w) & !is.dist.x) {
    if (!is.numeric(w) | length(w) != t.x) 
      stop("'w' should be a numeric vector of length = number of traits.", 
           "\n")
    else w <- w/sum(w)
  }
  if (missing(w)) 
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) 
        x[, x.class == "character"] <- as.factor(x[, 
                                                   x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1])) 
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          x.rn <- x.rn[-pos.NA]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)])) 
            stop("'dist.bin' must be an integer between 1 and 10.", 
                 "\n")
          x.dist <- dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = T, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.", 
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x))) 
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).", 
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist))) 
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0)) 
      stop("At least one species has no trait data.", "\n")
  }
  c <- dim(a)[1]
  if (!w.abun) 
    for (h in 1:c) {
      abpos <- which(a[h, ] > 0)
      a[h, abpos] <- 1
    }
  attr(x.dist, "Labels") <- x.rn
  if (is.euclid(x.dist)) 
    x.dist2 <- x.dist
  if (!is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- lingoes(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", 
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", 
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2)) 
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", 
             "\n")
      if (is.euclid(x.dist2)) 
        if (messages) 
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", 
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages) 
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", 
            "\n")
    }
  }
  x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  traits <- round(x.pco$li, .Machine$double.exponent)
  nb.sp <- numeric(c)
  for (i in 1:c) {
    sp.pres <- which(a[i, ] > 0)
    traits.sp.pres <- traits[sp.pres, , drop = F]
    traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < 
                     tol] <- 0
    nb.sp[i] <- nrow(unique(traits.sp.pres))
  }
  names(nb.sp) <- row.names(a)
  min.nb.sp <- min(nb.sp)
  if (min.nb.sp < 3) 
    if (messages) 
      cat("FEVe: Could not be calculated for communities with <3 functionally singular species.", 
          "\n")
  if (min.nb.sp < 2) 
    if (messages) 
      cat("FDis: Equals 0 in communities with only one functionally singular species.", 
          "\n")
  if (calc.FRic) {
    x.class2 <- sapply(x, data.class)
    if (all(x.class2 == "factor" | x.class2 == "ordered")) {
      if (length(x.class2) == 1 & x.class2[1] == "ordered") {
        traits.FRic1 <- rank(x[, 1])
        names(traits.FRic1) <- x.rn
        traits.FRic <- data.frame(traits.FRic1)
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one ordinal trait present in 'x'. FRic was measured as the range of the ranks, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when 'x' is a single ordinal trait.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      else {
        traits.FRic <- x
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only categorical and/or ordinal trait(s) present in 'x'. FRic was measured as the number of unique trait combinations, NOT as the convex hull volume.", 
              "\n")
        if (stand.FRic) 
          FRic.all <- nrow((unique(traits.FRic)))
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot be computed when only categorical and/or ordinal trait(s) present in 'x'.", 
                "\n")
        }
      }
    }
    else {
      if (x.pco$nf == 1) {
        traits.FRic <- x.pco$li
        qual.FRic = 1
        if (messages) 
          cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.", 
              "\n")
        if (calc.FDiv) {
          calc.FDiv <- FALSE
          if (messages) 
            cat("FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.", 
                "\n")
        }
        if (stand.FRic) {
          traits.range <- range(traits.FRic[, 1])
          FRic.all <- traits.range[2] - traits.range[1]
        }
      }
      if (x.pco$nf > 1) {
        warning <- FALSE
        m.max <- min.nb.sp - 1
        if (m == "min") {
          warning <- TRUE
          if (min.nb.sp < 4) {
            nb.sp2 <- nb.sp[nb.sp > 3]
            m.min <- floor(log2(min(nb.sp2)))
            if (messages) 
              cat("FRic: To respect s >= 2^t, FRic could not be calculated for communities with <4 functionally singular species.", 
                  "\n")
          }
          else m.min <- floor(log2(min.nb.sp))
        }
        else {
          if (min.nb.sp < 3) {
            nb.sp2 <- nb.sp[nb.sp > 2]
            m.max <- min(nb.sp2) - 1
            if (messages) 
              cat("FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species.", 
                  "\n")
          }
          else m.max <- m.max
        }
        if (is.numeric(m) & m <= 1) 
          stop("When 'm' is an integer, it must be >1.", 
               "\n")
        if (is.numeric(m) & m > m.max) 
          m <- m.max
        if (m == "min") 
          m <- m.min
        if (m == "max") 
          m <- m.max
        if (!is.numeric(m) & m != "min" & m != "max") 
          stop("'m' must be an integer >1, 'min', or 'max'.", 
               "\n")
        if (m < x.pco$nf) {
          traits.FRic <- x.pco$li[, 1:m]
          if (x.pco$nf - m == 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last PCoA axis (out of", 
                  x.pco$nf, "in total) was removed.", "\n")
          if (x.pco$nf - m > 1) 
            if (messages) 
              cat("FRic: Dimensionality reduction was required. The last", 
                  x.pco$nf - m, "PCoA axes (out of", x.pco$nf, 
                  "in total) were removed.", "\n")
          if (is.euclid(x.dist)) {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr != "none") {
            qual.FRic <- sum(x.pco$eig[1:m])/sum(x.pco$eig)
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (based on corrected distance matrix) =", 
                  qual.FRic, "\n")
          }
          if (!is.euclid(x.dist) & corr == "none") {
            delta <- -0.5 * bicenter.wt(x.dist * x.dist)
            lambda <- eigen(delta, symmetric = TRUE, 
                            only.values = TRUE)$values
            sum.m <- sum(lambda[1:m])
            sum.n <- sum(lambda)
            lambda.neg <- c(lambda[lambda < 0])
            max.neg <- abs(min(lambda.neg))
            qual.FRic <- (sum.m + (length(lambda[1:m]) * 
                                     max.neg))/(sum.n + ((length(lambda) - 1) * 
                                                           max.neg))
            if (messages) 
              cat("FRic: Quality of the reduced-space representation (taking into account the negative eigenvalues) =", 
                  qual.FRic, "\n")
          }
        }
        if (m >= x.pco$nf) {
          qual.FRic = 1
          traits.FRic <- x.pco$li
          if (x.pco$nf == 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. The 2 PCoA axes were kept as 'traits'.", 
                  "\n")
          if (x.pco$nf > 2) 
            if (messages) 
              cat("FRic: No dimensionality reduction was required. All", 
                  x.pco$nf, "PCoA axes were kept as 'traits'.", 
                  "\n")
        }
        if (stand.FRic) {
          hull.all <- convhulln(traits.FRic, "FA")
          FRic.all <- hull.all$vol
        }
      }
    }
  }
  if (!calc.FRic & calc.FDiv) 
    cat("FDiv: Cannot be computed when 'calc.FRic' is FALSE.", 
        "\n")
  if (calc.FRic & calc.FDiv) 
    if (min.nb.sp < 3) 
      if (messages) 
        cat("FDiv: Could not be calculated for communities with <3 functionally singular species.", 
            "\n")
  if (calc.FGR) {
    if (clust.type == "kmeans") {
      tr.clust <- cascadeKM(traits, km.inf.gr, km.sup.gr, 
                            km.iter, km.crit)
      cat("FGR: Summary of kmeans clustering\n")
      cat("\nPartition\n")
      print(tr.clust$partition)
      cat("\nResults\n")
      print(tr.clust$results)
      cat("\nSize\n")
      print(tr.clust$size)
      plot(tr.clust)
      part.names <- colnames(tr.clust$partition)
      part.names <- as.numeric(substr(part.names, 1, 1))
      cat("\nFGR: How many groups?", "\n")
      cut.g <- toupper(scan(file = "", what = "character", 
                            nlines = 1, quiet = T))
      cut.gr <- as.integer(cut.g)
      if (cut.gr < km.inf.gr | cut.gr > km.sup.gr) 
        stop("You must type an integer between 'km.ing.gr' and 'km.sup.gr'.", 
             "\n")
      spfgr.all <- tr.clust$partition[, part.names == cut.gr]
      names(spfgr.all) <- x.rn
    }
    else {
      tr.clust <- hclust(x.dist, method = clust.type)
      plot(tr.clust, main = "Cluster dengrogram of species based on functional traits")
      cat("FGR: Do you want to cut the dendrogram from height or from the number of groups? Type 'h' for height, 'g' for groups.", 
          "\n")
      cut <- toupper(scan(file = "", what = "character", 
                          nlines = 1, quiet = T))
      if (cut == "H") {
        cat("FGR: At what height do you want the dendrogram to be cut?", 
            "\n")
        cut.d <- toupper(scan(file = "", what = "character", 
                              nlines = 1, quiet = T))
        cut.dist <- as.numeric(cut.d)
        spfgr.all <- cutree(tr.clust, h = cut.dist)
      }
      if (cut == "G") {
        cat("FGR: How many groups?", "\n")
        cut.g <- toupper(scan(file = "", what = "character", 
                              nlines = 1, quiet = T))
        cut.gr <- as.integer(cut.g)
        spfgr.all <- cutree(tr.clust, k = cut.gr)
      }
      if (cut != "H" & cut != "G") 
        stop("You must type 'h' or 'g'", "\n")
    }
    a.t <- t(a)
    by.gr <- list(spfgr.all)
    gr.abun <- aggregate(a.t, by.gr, sum)
    lab <- paste("group", gr.abun[, 1], sep = "")
    gr.abun <- data.frame(t(gr.abun[, -1]))
    colnames(gr.abun) <- lab
    rownames(gr.abun) <- rownames(a)
  }
  if (is.matrix(x) | is.data.frame(x) & calc.CWM) {
    CWM <- functcomp(x, a, CWM.type = CWM.type)
  }
  if (calc.CWM & class(x)[1] == "dist" | class(x)[1] == "dissimilarity") 
    if (messages) 
      cat("CWM: When 'x' is a distance matrix, CWM cannot be calculated.", 
          "\n")
  divc <- function(df, dis = NULL, scale = FALSE) {
    if (!inherits(df, "data.frame")) 
      stop("Non convenient df")
    if (any(df < 0)) 
      stop("Negative value in df")
    if (!is.null(dis)) {
      if (!inherits(dis, "dist")) 
        stop("Object of class 'dist' expected for distance")
      dis <- as.matrix(dis)
      if (nrow(df) != nrow(dis)) 
        stop("Non convenient df")
      dis <- as.dist(dis)
    }
    if (is.null(dis)) 
      dis <- as.dist((matrix(1, nrow(df), nrow(df)) - diag(rep(1, 
                                                               nrow(df)))) * sqrt(2))
    div <- as.data.frame(rep(0, ncol(df)))
    names(div) <- "diversity"
    rownames(div) <- names(df)
    for (i in 1:ncol(df)) {
      if (sum(df[, i]) < 1e-16) 
        div[i, ] <- 0
      else div[i, ] <- (t(df[, i]) %*% (as.matrix(dis)^2) %*% 
                          df[, i])/2/(sum(df[, i])^2)
    }
    if (scale == TRUE) {
      divmax <- divcmax(dis)$value
      div <- div/divmax
    }
    return(div)
  }
  RaoQ <- divc(data.frame(t(a)), x.dist, scale = scale.RaoQ)
  RaoQ <- RaoQ[, 1]
  names(RaoQ) <- rownames(a)
  disp <- fdisp2(x.dist, a)
  FDis <- disp$FDis
  nbsp <- rep(NA, c)
  names(nbsp) <- row.names(a)
  FRic <- rep(NA, c)
  names(FRic) <- row.names(a)
  FEve <- rep(NA, c)
  names(FEve) <- row.names(a)
  FGR <- rep(NA, c)
  names(FGR) <- row.names(a)
  FDiv <- rep(NA, c)
  names(FDiv) <- row.names(a)
  for (i in 1:c) {
    sppres <- which(a[i, ] > 0)
    S <- length(sppres)
    nbsp[i] <- S
    tr <- data.frame(traits[sppres, ])
    if (calc.FRic) 
      tr.FRic <- data.frame(traits.FRic[sppres, ])
    ab <- as.matrix(a[i, sppres])
    abundrel <- ab/sum(ab)
    if (calc.FRic) {
      if (all(x.class2 == "factor" | x.class2 == "ordered")) {
        if (length(x.class2) == 1 & x.class2[1] == "ordered") {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
        else {
          if (!stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))
          if (stand.FRic) 
            FRic[i] <- nrow((unique(tr.FRic)))/FRic.all
        }
      }
      else {
        if (dim(tr.FRic)[2] > 1 & nb.sp[i] >= 3) {
          if (warning) 
            thresh <- 4
          if (!warning) 
            thresh <- 3
          if (nb.sp[i] >= thresh) {
            convhull <- convhulln(tr.FRic, "FA")
            if (!stand.FRic) 
              FRic[i] <- convhull$vol
            if (stand.FRic) 
              FRic[i] <- convhull$vol/FRic.all
          }
          else {
          }
        }
        if (dim(tr.FRic)[2] == 1) {
          tr.range <- range(tr.FRic[, 1])
          t.range <- tr.range[2] - tr.range[1]
          if (!stand.FRic) 
            FRic[i] <- t.range
          if (stand.FRic) 
            FRic[i] <- t.range/FRic.all
        }
      }
    }
    if (nb.sp[i] >= 3) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] + 
        abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (m in 1:((S - 1) * S/2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m]/(abund2vect[m])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    }
    if (calc.FDiv & calc.FRic) {
      if (any(x.class2 == "numeric") & dim(tr.FRic)[2] > 
          1 & nb.sp[i] >= 3) {
        vert0 <- convhulln(tr.FRic, "Fx TO 'vert.txt'")
        vert1 <- scan("vert.txt", quiet = T)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        trvertices <- tr.FRic[vertices, ]
        baryv <- apply(trvertices, 2, mean)
        distbaryv <- rep(0, S)
        for (j in 1:S) distbaryv[j] <- (sum((tr.FRic[j, 
                                                     ] - baryv)^2))^0.5
        meandB <- mean(distbaryv)
        devdB <- distbaryv - meandB
        abdev2 <- abundrel * devdB
        ababsdev2 <- abundrel * abs(devdB)
        FDiv[i] <- (sum(abdev2) + meandB)/(sum(ababsdev2) + 
                                             meandB)
      }
    }
    if (calc.FGR) 
      FGR[i] <- length(unique(spfgr.all[sppres]))
  }
  res <- list()
  res$nbsp <- nbsp
  res$sing.sp <- nb.sp
  if (calc.FRic) 
    res$FRic <- FRic
  if (calc.FRic) 
    res$qual.FRic <- qual.FRic
  res$FEve <- FEve
  if (calc.FDiv) 
    res$FDiv <- FDiv
  res$FDis <- FDis
  res$RaoQ <- RaoQ
  if (calc.FGR) {
    res$FGR <- FGR
    res$spfgr <- spfgr.all
    res$gr.abun <- gr.abun
  }
  if (is.matrix(x) | is.data.frame(x) & calc.CWM) 
    res$CWM <- CWM
  if (print.pco) {
    res$x.values <- x.pco$eig
    res$x.axes <- x.pco$li
  }
  invisible(res)
}



fdisp2 <-   function (d, a, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (!is.matrix(a)) 
    stop("'a' must be a matrix.")
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
         "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  #if (any(abun.sum == 0)) 
  #  stop("At least one community has zero-sum abundances (no species).", 
  #       "\n")
  abun.sum2 <- apply(a, 2, sum)
  #if (any(abun.sum2 == 0)) 
  #  stop("At least one species does not occur in any community (zero total abundance across all communities).", 
  #       "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}


#### Read in Data (need to clean up)####

#### R Script, with text from the Methods ####


######## READING IN OCCURRENCE DATA ########

#### "SC data were downloaded from the Neotoma Paleoecology Database #####
SC.data <- get_site(sitename = "Samwell Cave") 
SC.data <- get_dataset(SC.data)
SC.data <- get_download(SC.data)
SC.data <- SC.data[[4]] 

SC.time <- neo.SC$chronologies$`Blois et al. 2010`$age
SC.time <- SC.time/(-10000)

#### "EMC data...from Cuenca-Bescós et al. (Cuenca-Bescós et al., 2009)...Straus and González Morales (2003, 2007, 2016, 2018) and Straus et al. (2015)" ####
EMC.main <- read.csv("data/VanBuren_Jarzyna_EMC_DataFile.csv", header = T, stringsAsFactors = F, fileEncoding = "latin1")

EMC.data <- EMC.main[,4:ncol(EMC.main)]
EMC.data <- EMC.data[which(!is.na(EMC.data$Arvicola.terrestris)),] # Removes placeholder NA rows from main data file

#### "Each radiocarbon date at EMC was calibrated using the...package ‘Bchron’" ####
chron.dat <- EMC.main[,1:3]
chron.dat$median <- NA
mean <- chron.dat[,2]
sd <- chron.dat[,3]
for(i in 1:(length(mean) - 1)){
  if(!is.na(sd[i])){
    newage <- BchronCalibrate(ages = mean[i], ageSds = sd[i], calCurves = "intcal13")
    age_samples <- sampleAges(newage, n_samp = 20000) # "A distribution of 20,000 ages for each layer..."
    med.age <- median(age_samples) # "...the median value was extracted from this distribution..."
    chron.dat$median[i] <- med.age
  }
}

ords <- unique(chron.dat$Level)
chron.dat2 <- chron.dat %>%
  group_by(Level) %>%
  summarise(median.age = mean(median, na.rm = T))
chron.dat2 <- as.data.frame(chron.dat2)
chron.dat2 <- chron.dat2[match(ords, chron.dat2$Level),]

#### "We removed the two oldest dated layers...the most recent layer...[and] layer 115" ####
EMC.time1 <- chron.dat2$median.age 
EMC.time1 <- EMC.time1[-c(1, 32, 41:43)] 
EMC.time1 <- EMC.time1/(-10000) 
EMC.time <- EMC.time1[which(!is.na(EMC.time1))] # keep EMC.time1 to use below

#### Cleaning up Samwell Cave download data ####
SC.spp <- SC.data$taxon.list
SC.counts <- SC.data$counts

dels <- c("Carnivora", "Chiroptera", "Leporidae", "Artiodactyla") # "we...removed all occurrences not resolved to at least the genus level

SC.counts <- SC.counts[,str_detect(SC.spp$taxon.name, paste(dels, collapse = '|')) == F]
SC.spp <- SC.spp[str_detect(SC.spp$taxon.name, paste(dels, collapse = '|')) == F,]

#### Cleaning up El Mirón Cave data ####
exclude <- c("Miniopterus", # "we excluded bats from the analysis..."
             "Myotis",
             "Chiroptera",
             "Lagomorpha",  # "Lagamorpha...were removed"
             "Pliomys") # "Pliomys lenki...is extinct"

EMC.data <- EMC.data[,str_detect(colnames(EMC.data), paste(exclude, collapse = '|')) == F]

EMC.data <- EMC.data[-c(1, 32, 41:43),] # Remove layers with time issues, as above
EMC.data <- EMC.data[-which(is.na(EMC.time1)),] # Remove layers without time, as above

#### Replacing colnames in SC.counts with names that match Elton Traits, including removing 'cf' from some occurrences ####
dels <- c("Scapanus latimanus", "Arborimus albipes", "megalotis",
         "sabrinus", "californicus", "Thomomys sp.", "Thomomys cf. T. mazama", "Tamias ")
repl <- c("Scapanus latimanus", "Arborimus albipes", "Reithrodontomys megalotis",
          "Glaucomys sabrinus", "Myodes californicus", "Thomomys", "Thomomys mazama", "Tamias")

levels(SC.spp$taxon.name) <- c(levels(SC.spp$taxon.name), repl)

for(i in 1:length(dels)){
  SC.spp$taxon.name[str_detect(SC.spp$taxon.name, dels[i])] <- repl[i]
}

SC.spp$taxon.name <- droplevels(SC.spp$taxon.name)

## Simplifying the names of genus-level occurrences; Thomomys is simplified above because some occurrences are resolved to the species level; Tamias is used above because of how commonly "Tamias" is used
genera <- c("Microtus", "Neotoma", "Peromyscus", "Sorex")

levels(SC.spp$taxon.name) <- c(levels(SC.spp$taxon.name), genera)
for(i in 1:length(genera)){
  SC.spp$taxon.name[str_detect(SC.spp$taxon.name, genera[i])] <- genera[i]
}
SC.spp$taxon.name <- droplevels(SC.spp$taxon.name)

colnames(SC.counts) <- SC.spp$taxon.name

#### Replacing colnames in EMC.data with names that match Elton Traits ####
EMC.spp <- c("Arvicola amphibius",
             "Arvicola sapidus",
             "Microtus lusitanicus",
             "Microtus oeconomus",
             "Microtus agrestis",
             "Microtus arvalis",
             "Chionomys nivalis",
             "Myodes glareolus",
             "Microtus gregalis",
             "Apodemus sylvaticus",
             "Eliomys quercinus",
             "Glis glis",
             "Sciurus vulgaris",
             "Sorex minutus",
             "Sorex coronatus",
             "Neomys fodiens",
             "Crocidura russula",
             "Galemys pyrenaicus",
             "Talpa europaea",
             "Erinaceus europaeus",
             "Mustela nivalis")

colnames(EMC.data) <- EMC.spp

EMC.data <- EMC.data[,which(colSums(EMC.data) > 0)] # Remove any additional columns

EMC.data <- EMC.data[,(order(colnames(EMC.data)))]
EMC.rel <- t(apply(EMC.data, 1, function(x) x/sum(x, na.rm = T))) #relative abundance
EMC.rel <- EMC.rel[-c(1, 32, 41:43),] # "We removed the two oldest dated layers...the most recent layer...[and] layer 115"
EMC.rel <- EMC.rel[which(!is.na(EMC.time)),]


#### ####


######## REGIONAL SPECIES POOL AND TRAITS ########

#### "We downloaded range maps...from the IUCN Red List" ####

SC_IUCN <- readOGR(dsn = "data/SC_MC_IUCN/")

EMC_IUCN <- readOGR("data/ElMironCave_pool/")

#### "...that extended +/-1° Latitude and +/-1° Longitude from each site" ####
#SC_coords <- SpatialPoints(coords = data.frame(-122.2378, 40.91691)) # coordinates of Samwell Cave
#EMC_coords <- SpatialPoints(coords = data.frame(-3.45111, 43.246389)) #coordinates of El Mirón Cave

SC_polygon_coords <- data.frame(lon = c(-123.2378, -123.2378, -121.2378, -121.2378, -123.2378),
                               lat = c(39.91691, 41.91691, 41.91691, 39.91691, 39.91691))

SC_polygon <- Polygons(list(Polygon(SC_polygon_coords)), ID = 1)
SC_polygon <- SpatialPolygons(list(SC_polygon))
proj4string(SC_polygon) <- proj4string(SC_IUCN)

SC_regional <- gIntersects(SC_IUCN, SC_polygon, byid = T) #stored as sc_reg
SC_regional <- unique(SC_IUCN@data$BINOMIAL[colSums(SC_regional) > 0]) # drop species that do not occur in the polygon
SC_regional <- droplevels(SC_regional)


EMC_polygon_coords <- data.frame(lon = c(-4.451111, -4.451111, -2.451111, -2.451111, -4.451111),
                                lat = c(42.246389, 44.246389, 44.246389, 42.246389, 42.246389))

EMC_polygon <- Polygons(list(Polygon(EMC_polygon_coords)), ID = 1)
EMC_polygon <- SpatialPolygons(list(EMC_polygon))
proj4string(EMC_polygon) <- proj4string(EMC_IUCN)

EMC_regional <- gIntersects(EMC_IUCN, EMC_polygon, byid = T) #stored as emc_reg
EMC_regional <- unique(EMC_IUCN@data$BINOMIAL[colSums(EMC_regional) > 0])
EMC_regional <- droplevels(EMC_regional)





#### "We based estimates of all...diversity metrics on...traits in Wilman et al. 2014" ####
Elton <- read.delim("data/EltonTraits.txt")
Elton <- Elton[order(Elton$Scientific),]

#### Aligning SC_regional and Elton Traits ####
SC_regional <- SC_regional[SC_regional %in% "Cervus canadensis" == F] # Needs to be removed and is not a match with Elton Traits
errors <- SC_regional[SC_regional %in% Elton$Scientific == F]
repls <- c("Spermophilus beecheyi", # Replacing species names with those used by Elton Traits
           "Spermophilus beldingi",
           "Spermophilus lateralis",
           "Tamias quadrimaculatus",
           "Tamias siskiyou",
           "Tamias senex",
           "Tamias amoenus",
           "Tamias minimus",
           "Tamias sonomae",
           "Tamias speciosus")
levels(SC_regional) <- c(levels(SC_regional), repls)

for(i in 1:length(repls)){
  SC_regional[SC_regional %in% errors[i]] <- repls[i]
}

SC_regional <- droplevels(SC_regional)

SC_regional <- as.character(SC_regional)
SC_regional <- SC_regional[order(SC_regional)] 

#### Aligning EMC_regional and Elton Traits #### 
errors <- EMC_regional[EMC_regional %in% Elton$Scientific == F] #a bat and a whale, so they are removed
EMC_regional <- EMC_regional[EMC_regional %in% errors == F]
EMC_regional <- droplevels(EMC_regional)

EMC_regional <- as.character(EMC_regional)
EMC_regional <- EMC_regional[order(EMC_regional)] 

#### "...all bats and mammals > 1,100g were removed from the dataset." ####

## Samwell Cave
SC_traits <- Elton[Elton$Scientific %in% SC_regional,]

SC_regional <- SC_regional[-which(SC_traits$BodyMass.Value > 1100)]
SC_traits <- SC_traits[-which(SC_traits$BodyMass.Value > 1100),]

bats <- c("Vespertilionidae", "Molossidae") # "we excluded bats from the analysis..."
SC_regional <- SC_regional[SC_traits$MSWFamilyLatin %in% bats == F]
SC_traits <- SC_traits[SC_traits$MSWFamilyLatin %in% bats == F,]
rownames(SC_traits) <- SC_traits$Scientific

SC_traits <- SC_traits[,c(4:13,16,19:21, 24)] # Only the columns with trait data
SC_traits[,11] <- as.character(SC_traits[,11]) # This is required by the FD function later


## El Miron Cave
EMC_traits <- Elton[Elton$Scientific %in% EMC_regional,]


EMC_regional <- EMC_regional[-which(EMC_traits$BodyMass.Value > 1100)]
EMC_traits <- EMC_traits[-which(EMC_traits$BodyMass.Value > 1100),]

bats <- c("Vespertilionidae", "Molossidae", "Rhinolophidae")
EMC_regional <- EMC_regional[EMC_traits$MSWFamilyLatin %in% bats == F] # "we excluded bats from the analysis..."
EMC_traits <- EMC_traits[EMC_traits$MSWFamilyLatin %in% bats == F,]
rownames(EMC_traits) <- EMC_traits$Scientific
EMC_traits <- EMC_traits[,c(4:13,16,19:21, 24)] # Only the columns with trait data
EMC_traits[,11] <- as.character(EMC_traits[,11]) # This is required by the FD function later

#### "For all genus-level occurrences, genus-level mean trait values were used." ####

SC_genera <- c("Tamias ",
               "Thomomys",
               "Sorex",
               "Peromyscus",
               "Neotoma",
               "Microtus")

Thomomys <- c("Thomomys bottae", "Thomomys mazama")
Thoms <- suppressWarnings(SC_traits[str_detect(rownames(SC_traits), Thomomys),]) # Keep individual trait data for these two species of Thomomys

tr.hold <- NULL
for(i in 1:length(SC_genera)){ 
  hold <- SC_traits[str_detect(rownames(SC_traits), SC_genera[i]),]
  SC_traits <- SC_traits[str_detect(rownames(SC_traits), SC_genera[i]) == F,]
  x <- data.frame(c(apply(hold[,1:10], 2, mean),
                    unique(hold[,11:14]),
                    mean(hold[,15])))
  colnames(x) <- colnames(SC_traits)
  rownames(x) <- SC_genera[i]
  tr.hold <- rbind(tr.hold, x)
}

SC_traits <- rbind(SC_traits, tr.hold, Thoms) # Add trait data back for Thomomys bottae and T. mazama
SC_traits <- SC_traits[order(rownames(SC_traits)),]
rownames(SC_traits)[rownames(SC_traits) %in% "Tamias " == T] <- "Tamias" # correct this from above
#saveRDS(SC_traits, file = "SC_grid_rsp.rds")

#### "The taxonomic list [at EMC] includes occurrences for Apodemus gr. sylvaticus-flavicollis and Sorex gr. coronatus-araneus" ####

sorex <- c("Sorex coronatus", "Sorex araneus")
apedomus <- c("Apodemus sylvaticus", "Apodemus flavicollis")

sor.tr <- Elton[str_detect(Elton$Scientific, paste(sorex, collapse = '|')),]
apo.tr <- Elton[str_detect(Elton$Scientific, paste(apedomus, collapse = '|')),]

### "...we used mean trait values of the two species" ####
EMC_traits$BodyMass.Value[rownames(EMC_traits = "Apodemis sylvaticus")] <- mean(apo.tr$BodyMass.Value) # These species only differ in body mass, so their average trait values for other traits are the same
EMC_traits$BodyMass.Value[rownames(EMC_traits) == "Sorex coronatus"] <- mean(sor.tr$BodyMass.Value)
#EMC_traits <- MC_trs[-1,]
#saveRDS(EMC_traitss, "EMC_traits.rds")


#### "We then expanded the regional species pool to include species...present in the fossil record...but have since been...extirpated..." ####

hold.names <- colnames(SC.counts)
SC.counts <- data.frame(SC.counts)
colnames(SC.counts) <- hold.names

SC.counts[,(ncol(SC.counts) + 1):nrow(SC_traits)] <- NA #ncol(SC.counts) + 1 = 19, used below
colnames(SC.counts)[19:ncol(SC.counts)] <- rownames(SC_traits)[rownames(SC_traits) %in% colnames(SC.counts) == F]
SC.counts[is.na(SC.counts)] <- 0
SC.counts <- SC.counts[,order(colnames(SC.counts))]
#saveRDS(counts, "SC_rawcounts_grid.rds")



EMC.data[,colnames(EMC.data) %in% rownames(EMC_traits) == F] #Arvicola amphibius and Microtus oeconomus not in Regional Pool

names <- c("Arvicola amphibius", "Microtus oeconomus") #adding them to the EMC_trait dataset
add <- Elton[Elton$Scientific %in% names == T,]
add <- add[,c(4:13,16,19:21, 24)]
add[,11] <- as.character(add[,11])
rownames(add) <- names
EMC_traits <- rbind(EMC_traits, add)
EMC_traits <- EMC_traits[order(rownames(EMC_traits)),]

EMC.data[,(ncol(EMC.data) + 1):nrow(EMC_traits)] <- NA
colnames(EMC.data)[20:ncol(EMC.data)] <- rownames(EMC_traits)[rownames(EMC_traits) %in% colnames(EMC.data) == F] #Expanding the regional pool
EMC.data[is.na(EMC.data)] <- 0
EMC.data <- EMC.data[,order(colnames(EMC.data))]
#saveRDS(EMC.data, "EMC_counts.rds")






#### ####


######## CALCULATING FUNCTIONAL DIVERSITY - WARNING: COMPUTATIONALLY INTENSIVE ########

#### "At SC...we subsampled to...130...a layer" ####
data <- SC_counts
time <- SC.time
traits <- SC_traits
Weights <- c(rep(.25/10, 10), .25, rep(.25/3, 3), .25) # "Equal weights were given to each of the trait categories and to each axis within the trait categories..."


data.sub <- data0[which(rowSums(data0) > 129),]
time.sub <- time[which(rowSums(data0) > 129)]


sub.rich <- NULL
sub.FD <- list()
null.FD <- list()
abun <- list()
MIN <- 130
for(t in 1:500){ # "Subsampling was repeated 500 times."
  sample.list <- list()
  
  for(i in 1:nrow(data.sub50)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = MIN, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub)[str_detect(colnames(data.sub), 
                                                                                  paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3) # (m = 3) = "...we based FR estimation on the first three principal coordinate axes..."
    for(i in 1:100){ # "...we randomized (100 times) the abundances across species in the regional species pool
      x2 <- apply(abun.mat, 1, sample)
      x2 <- t(x2)
      colnames(x2) <- colnames(abun.mat)
      x2 <- as.data.frame(x2)
      null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
      null.FD[[(t - 1)*100 + i]] <- null
    }
  print(t)
}

SC_FD <- sub.FD
SC_nullFD <- null.FD
SC_communities <- abun

#write_rds(sub.FD, "SC_grid_FD.rds")
#write_rds(null.FD, "SC_grid_FD_null.rds")
#write_rds(abun, "SC_grid_communities.rds")



#### "At EMC...we used...45 individuals...as the sampling target" ####

data <- EMC.data
time <- EMC.time
traits <- EMC_traits
Weights <- c(rep(.25/10, 10), .25, rep(.25/3, 3), .25) # "Equal weights were given to each of the trait categories and to each axis within the trait categories..."

data.sub <- data[which(rowSums(data) > 44),] # "All layers with fewer than 45 individuals at EMC were excluded..." 
time.sub <- time[which(rowSums(data) > 44)]

EMC.timesub <- time.sub # to use later since not all layers at EMC were included in the analysis

sub.rich <- NULL
sub.FD <- list()
null.FD <- list()
abun <- list()
MIN <- 45
for(t in 1:500){ # "Subsampling was repeated 500 times."
  sample.list <- list()
  
  for(i in 1:nrow(data.sub)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = MIN, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub)[str_detect(colnames(data.sub), 
                                                                                  paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3) # (m = 3) = "...we based FR estimation on the first three principal coordinate axes..."
    for(i in 1:100){ # "...we randomized (100 times) the abundances across species in the regional species pool
      x2 <- apply(abun.mat, 1, sample)
      x2 <- t(x2)
      colnames(x2) <- colnames(abun.mat)
      x2 <- as.data.frame(x2)
      null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
      null.FD[[(t - 1)*100 + i]] <- null
    }
  print(t)
}


EMC_FD <- sub.FD
EMC_nullFD <- null.FD
EMC_communities <- abun

#write_rds(EMC_FD, "EMC_FD.rds")
#write_rds(EMC_nullFD, "EMC_FD_null.rds")
#write_rds(EMC_communities, "EMC_communities.rds")



#### ####


######## MODEL TESTING (with time) ########

#### "...we tested...a linear relationship...and...a quadratic relationship..." ####
#### Samwell Cave ####

FD.dat <- SC_FD
time <- SC.time

waic <- data.frame(matrix(NA, 100, 2))

## FRic
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_FR_waic.rds")

## FDis

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_FDS_waic.rds")

## FDiv

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_FDV_waic.rds")

## FEve
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_FE_waic.rds")

#### Traits ####
### Diet ####

waic <- data.frame(matrix(NA, 100, 2))

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Fruit_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Inv_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_PlantO_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Scav_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Seed_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Vend_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vunk/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_Vunk_waic.rds")

### Foraging Stratum ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_SA_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_SG_waic.rds")

### Activity Pattern ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_AC_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_AD_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_AN_waic.rds")

### Body Mass ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
#saveRDS(waic, "SC_BM_waic.rds")

#### El Miron Cave #### 

FD.dat <- EMC_FD
time <- EMC.timesub

waic <- data.frame(matrix(NA, 100, 2))

## FRic
waic <- data.frame(matrix(NA, 100, 2))

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_FR_waic.rds")

## FDis
waic <- data.frame(matrix(NA, 100, 2))

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_FDS_waic.rds")

## FDiv
waic <- data.frame(matrix(NA, 100, 2))

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_FDV_waic.rds")

## FEve 
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_FE_waic.rds")

#### Traits ####
### Diet ####

waic <- data.frame(matrix(NA, 100, 2))

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Fruit_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Inv_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_PlantO_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Scav_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Seed_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vect/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Vect_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Vend_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vfish/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_Vfish_waic.rds")

### Foraging Stratum ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_SA_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_SG_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_S,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_SS_waic.rds")

### Activity Pattern ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_AC_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_AD_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_AN_waic.rds")

### Body Mass ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  form2 <- FD ~ time + I(time^2)
  mod.i2 <- inla(form2, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                 control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.i$waic$waic
  waic[i,2] <- mod.i2$waic$waic
}

colnames(waic) <- c("lin", "quad")
saveRDS(waic, "EMC_BM_waic.rds")

#### WAIC comparing ####
waic <- readRDS("SC_FR_waic.rds") # use whichever input you're interested in analyzing
length(which(waic[,1] <= (waic[,2] + 7)))

#### ####


######## REGRESSIONS AGAINST TIME ########

#### "The linear model was preferred..." #### 
func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope <- data.frame(matrix(NA,500,7))

#### Samwell Cave ####
FD.dat <- SC_FD
time <- SC.time

## FRic
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FR_fit.rds")

## FDis
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FDS_fit.rds")

## FDiv
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FDV_fit.rds")

## FEve
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FE_fit.rds")

#### Traits ####
### Diet ####
func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope <- data.frame(matrix(NA,500,7))

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Fruit_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Inv_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_PlantO_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Scav_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Seed_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Vend_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vunk/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_Vunk_fit.rds")


### Foraging Stratum ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_SA_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_SG_fit.rds")


### Activity pattern ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_AC_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_AD_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_AN_fit.rds")


### Body Mass ####
for (i in 1:500){
  data <- data.frame("FD" = log(FD.dat[[i]]$CWM$BodyMass.Value),
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_BM_fit.rds")


#### El Mirón Cave ####
FD.dat <- EMC_FD
time <- EMC.timesub

## FRic
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FR_fit.rds")


## FDis
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FDS_fit.rds")

## FDiv
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FDV_fit.rds")

## FEve 
func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope <- data.frame(matrix(NA,500,7))

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  hold <- mod.i$summary.fitted.values
  colnames(hold) <- c("mean","sd","lci","med","uci","mode")
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FE_fit.rds")

#### Traits ####
### Diet ####
func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope <- data.frame(matrix(NA,500,7))

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_feve_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #feve[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Fruit_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_feve_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #feve[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Inv_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_PlantO_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Scav_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Seed_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vect/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Vect_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Vend_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vfish/100,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_Vfish_fit.rds")



### Foraging Stratum ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_SA_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_SG_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_S,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_SS_fit.rds")



### Activity pattern ####
func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope <- data.frame(matrix(NA,500,7))

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_AC_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_AD_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_AN_fit.rds")



### Body Mass ####
for (i in 1:500){
  data <- data.frame("FD" = log(FD.dat[[i]]$CWM$BodyMass.Value),
                     "time" = time)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_BM_fit.rds")



#### ####


######## NULL MODELS (time) - WARNING: COMPUTATIONALLY INTENSIVE ########

func <- list()
func$intercept <- data.frame(matrix(NA,50000,7))
func$slope <- data.frame(matrix(NA,50000,7))

#### Samwell Cave ####
FD.dat <- SC_nullFD
time <- SC.time

for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FR_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FDS_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FDV_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="SC_FE_null.rds")

#### El Mirón Cave ####
FD.dat <- EMC_nullFD
time <- EMC.timesub

func <- list()
func$intercept <- data.frame(matrix(NA,50000,7))
func$slope <- data.frame(matrix(NA,50000,7))

for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FR_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FDS_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FDV_null.rds")


for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "time" = time)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ time
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  #saveRDS(mod.i, file=paste0("SC_func_",i,".rds") #if we want to save all the models (they're small, but will add up)
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  #func[[i]] <- mod.i
}

colnames(func[[1]]) <- colnames(mod.i$summary.fixed)
colnames(func[[2]]) <- colnames(mod.i$summary.fixed)
#saveRDS(func, file="EMC_FE_null.rds")

#### ####


######## ENVIRONMENTAL DATA (from Beyer et al. 2020) ########
env <- readRDS("data/ENV4REGS_SC.rds")
x <- env[,c(2,3,5,6)] # "...mean annual temperature, annual precipitation, leaf area index (LAI), net primary productivity (NPP)..."
x <- apply(x, 2, scale) # "...transformed...into z-scores..."
SC_env <- princomp(x) # "...summarized using a principal components analysis..."

env <- readRDS("data/ENV4REGS_EMC.rds")
x <- env[,c(2,3,5,6)] # "...mean annual temperature, annual precipitation, leaf area index (LAI), net primary productivity (NPP)..."
x <- apply(x, 2, scale) # "...transformed...into z-scores..."
EMC_env <- princomp(x) # "...summarized using a principal components analysis..."

#### ####


######## MODEL TESTING (with environment) ########

#### "For each metric, we fitted 500 models with additive and 500 models with interaction effects..." ####
#### Samwell Cave ####
FD.dat <- SC_FD
env <- SC_env

axis1 <- env$scores[,1]
axis2 <- env$scores[,2]

waic <- data.frame(matrix(NA, 100, 2))

#### FRic
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  
  mod.1 <- inla(form1, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_FR_ENV_waic.rds")

#### FDis
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_FDS_ENV_waic.rds")

#### FDiv 

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_FDV_ENV_waic.rds")

#### FEve
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_FE_ENV_waic.rds")


#### Traits ####
### Diet ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Fruit_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Inv_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_PlantO_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Scav_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Seed_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Vend_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vunk/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_Vunk_ENV_waic.rds")


### Foraging Stratum ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_SA_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_SG_ENV_waic.rds")


### Activity Pattern ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_AC_ENV_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_AD_ENV_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_AN_ENV_waic.rds")


### Body Mass ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "SC_BM_ENV_waic.rds")



#### El Mirón Cave ####
FD.dat <- EMC_FD
env <- EMC_env

axis1 <- env$scores[,1]
axis2 <- env$scores[,2]

waic <- data.frame(matrix(NA, 100, 2))

#### FRic
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  
  mod.1 <- inla(form1, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_FR_ENV_waic.rds")


#### FDis
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_FDS_ENV_waic.rds")

#### FDiv 

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2

  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_FDV_ENV_waic.rds")

#### FEve
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_FE_ENV_waic.rds")


#### Traits ####
#### Diet ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Fruit_ENV_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Inv_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_PlantO_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Scav_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Seed_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vect/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Vect_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Vend_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vfish/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_Vfish_ENV_waic.rds")


#### Foraging Stratum ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_SA_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_SG_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_S,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_SS_ENV_waic.rds")


#### Activity Pattern ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_AC_ENV_waic.rds")


for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_AD_ENV_waic.rds")

for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_AN_ENV_waic.rds")

#### Body Mass ####
for (i in 1:100){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form1 <- FD ~ axis1 + axis2
  form2 <- FD ~ axis1 + axis2 + axis1*axis2
  mod.1 <- inla(form1, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  mod.2 <- inla(form2, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  waic[i,1] <- mod.1$waic$waic
  waic[i,2] <- mod.2$waic$waic
}
colnames(waic) <- c("model1", "model2")

#saveRDS(waic, "EMC_BM_ENV_waic.rds")



#### ####


######## REGRESSIONS AGAINST ENVIRONMENT ########

#### "The model that captured an additive effect of PC1 and PC2 was superior..." ####

func <- list()
func$intercept <- data.frame(matrix(NA,500,7))
func$slope1 <- data.frame(matrix(NA,500,7))
func$slope2 <- data.frame(matrix,500,7)

#### Samwell Cave ####
FD.dat <- SC_FD
env <- SC_env

axis1 <- env$scores[,1]
axis2 <- env$scores[,2]

## FRic
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func,"SC_FR_ENV_fit.rfs")

## FDis
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func,"SC_FDS_ENV_fit.rfs")

## FDiv
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func,"SC_FDV_ENV_fit.rfs")

## FEve
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func,"SC_FE_ENV_fit.rfs")

#### Traits ####
### Diet ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Fruit_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Inv_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_PlantO_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Scav_ENV_fit.rds")

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Seed_ENV_fit.rds")

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Vend_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vunk/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_Vunk_ENV_fit.rds")

### Foraging Stratum ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_SA_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_SG_ENV_fit.rds")



### Activity Pattern ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_AC_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_AD_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_AN_ENV_fit.rds")

### Body Mass ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "SC_BM_ENV_fit.rds")

#### El Mirón Cave ####
FD.dat <- EMC_FD
env <- EMC_env

axis1 <- env$scores[,1]
axis2 <- env$scores[,2]

## FRic
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_FR_ENV_fit.rfs")

## FDis
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_FDS_ENV_fit.rfs")

## FDiv
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_FDV_ENV_fit.rfs")

## FEve
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_FE_ENV_fit.rfs")

#### Traits ####
### Diet ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Fruit/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Fruit_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Inv/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Inv_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.PlantO/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_PlantO_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Scav/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Scav_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Seed/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Seed_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vect/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Vect_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vend/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Vend_ENV_fit.rds")

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Diet.Vfish/100,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_Vfish_ENV_fit.rds")

### Foraging Stratum ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_Ar,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_SA_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_G,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_SG_ENV_fit.rds")

for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$ForStrat.Value_S,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_SS_ENV_fit.rds")

### Activity Pattern #### 
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Crepuscular_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_AC_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Diurnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_AD_ENV_fit.rds")


for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$Activity.Nocturnal_1,
                     "axis1" = axis1,
                     "axis2" = axis2)
  data$FD[data$FD == 0] <- 0.00001
  data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_AN_ENV_fit.rds")


### Body Mass ####
for (i in 1:500){
  data <- data.frame("FD" = FD.dat[[i]]$CWM$BodyMass.Value,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="lognormal", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  func[[1]][i,] <- mod.i$summary.fixed[1,]
  func[[2]][i,] <- mod.i$summary.fixed[2,]
  func[[3]][i,] <- mod.i$summary.fixed[3,]
}

#saveRDS(func, "EMC_BM_ENV_fit.rds")


#### ####


######## NULL MODELS (environment) - WARNING: COMPUTATIONALLY INTENSIVE ########

res <- list()
res$intercept <- data.frame(matrix(NA,50000,7))
res$slope1 <- data.frame(matrix(NA,50000,7))
res$slope2 <- data.frame(matrix(NA,50000,7))

#### Samwell Cave ####
FD.dat <- SC_nullFD
env <- SC_env

### FRic
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="SC_FD_ENV_null.rds")

### FDis
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="SC_FR_ENV_null.rds")

### FDiv
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="SC_FDV_ENV_null.rds")

### FEve
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="SC_FE_ENV_null.rds")

#### El Mirón Cave ####
FD.dat <- EMC_nullFD
time <- EMC.timesub

### FRic
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FRic,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="EMC_FR_ENV_null.rds")

### FDis
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDis,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="gamma", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="EMC_FDS_ENV_null.rds")


### FDiv
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FDiv,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="EMC_FDV_ENV_null.rds")

### FEve
for (i in 1:50000){
  data <- data.frame("FD" = FD.dat[[i]]$FEve,
                     "axis1" = axis1,
                     "axis2" = axis2)
  #data$FD[data$FD == 0] <- 0.00001
  #data$FD[data$FD == 1] <- 0.99999
  form <- FD ~ axis1 + axis2
  mod.i <- inla(form, data=data, family="beta", control.compute = list(dic = TRUE, waic = TRUE), 
                control.predictor=list(compute=TRUE, link=1))
  res[[1]][i,] <- mod.i$summary.fixed[1,]
  res[[2]][i,] <- mod.i$summary.fixed[2,]
  res[[3]][i,] <- mod.i$summary.fixed[3,]
}

colnames(res[[1]]) <- colnames(mod.i$summary.fixed)
colnames(res[[2]]) <- colnames(mod.i$summary.fixed)
colnames(res[[3]]) <- colnames(mod.i$summary.fixed)
#saveRDS(res, file="EMC_FE_ENV_null.rds")


#### ####


######## CODE TO CALCULATE SIGNIFICANCE ########
# Regressions against time (1 slope)
obs <- readRDS("SC_FR_fit.rds") # Use whichever output combination you are interested in looking at
null <- readRDS("SC_FR_null.rds")

slope1 <- obs[[2]]
null1 <- null[[2]]

# Regressions against environment (2 slopes)
obs <- readRDS("SC_FR_ENV_fit.rds") # Use whichever output combination you are interested in looking at
null <- readRDS("SC_FR_ENV_null.rds")

slope1 <- obs[[2]]
null1 <- null[[2]]

slope2 <- obs[[3]]
null2 <- null[[3]]


#### Simple Majority ####
pos <- length(which(slope1[,3]) > 0) # number with lower CI above zero, i.e. positive slope
neg <- length(which(slope1[,5]) < 0) # number with upper CI below zero, i.e. negative slope


#### Standardized Effect Size and Quantile Score ####

slope <- slope1
null <- null1

SES <- NULL
QuantScore <- NULL
for(i in 1:100){
  obs.sub <- slope[i,4]
  null.sub <- NULL
  for(t in 1:100){
    null.sub[t] <- null[(i-1)*100 + t, 4]
  }
  QuantScore[i] <- rank(c(obs.sub, null.sub))[1]/101
  SES[i] <- (obs.sub + mean(null.sub))/sd(null.sub)
}

mean(QuantScore)
mean(SES)



#### ####


######## ANALYSES: SUPPORTING INFORMATION ########
#### "...we reduced the taxonomic resolution at EMC to the genus level..." ####

EMC.genus <- EMC.data
time3 #from earlier code
MC.rsp.tr #from earlier code

genera <- c("Microtus", "Sorex", "Arvicola") # "The genera at EMC for which multiple species were identified—Arvicola, Microtus, and Sorex..."

for(i in 1:3){
  colnames(EMC.genus)[str_detect(colnames(EMC.genus), genera[i])] <- genera[i]
}

t_dat <- t(EMC.genus)

EMC.genus <- aggregate(t_dat, by = list(rownames(t_dat)), sum)
colnames(EMC.genus) <- NULL
rownames(EMC.genus) <- EMC.genus[,1]
EMC.genus <- EMC.genus[,-1]
EMC.genus <- t(EMC.genus)
median(rowSums(EMC.genus)) #still 45

#saveRDS(dat2, "EMC_Genus_counts.rds")

EMC_traitgen <- EMC_traits
tr.hold <- NULL

for(i in 1:3){
  hold <- EMC_traitgen[str_detect(rownames(EMC_traitgen), genera[i]),]
  EMC_traitgen <- EMC_traitgen[str_detect(rownames(EMC_traitgen), genera[i]) == F,]
  x <- data.frame(c(apply(hold[,1:10], 2, mean),
                    unique(hold[,11:14]),
                    mean(hold[,15])))
  colnames(x) <- colnames(EMC_traitgen)
  rownames(x) <- genera[i]
  tr.hold <- rbind(tr.hold, x)
}

EMC_traitgen <- rbind(EMC_traitgen, tr.hold)
EMC_traitgen <- EMC_traitgen[order(rownames(EMC_traitgen)),]

#saveRDS(EMC_traitgen, "EMC_Genus_Traits.rds")


data.sub <- EMC_genus
traits <- EMC_traitgen

data.sub <- data.sub[which(rowSums(data.sub) > 44),]


sub.rich <- NULL
sub.FD <- list()
#null.FD <- list()
abun <- list()
for(t in 1:500){
  sample.list <- list()
  
  for(i in 1:nrow(data.sub)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = 45, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub50)[str_detect(colnames(data.sub), 
                                                                                  paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3)
  #  for(i in 1:100){                                  # If you are interested in null models
  #    x2 <- apply(abun.mat, 1, sample)
  #    x2 <- t(x2)
  #    colnames(x2) <- colnames(abun.mat)
  #    x2 <- as.data.frame(x2)
  #    null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
  #    null.FD[[(t - 1)*100 + i]] <- null
  #  }
  print(t)
}

EMC_genFD <- sub.FD

#saveRDS(EMC_genFD, "EMC_Genus_grid_FD.rds")



#### "We quantifyied NISP and MNI-based metrics for SC..." ####
# "...after removing all Thomomys sp., Sorex sp., and Scapanus latimanus observations..."
#### MNI ####
SC_MNI <- read.csv("data/SC_MNI_counts.csv") # From supplemental info in Blois et al. (2010)
traits <- SC_traits


colnames(SC_MNI) <- rownames(SC_traits)
SC_MNI <- SC_MNI[,-which(colnames(SC_MNI) %in% c("Thomomys", "Sorex", "Scapanus latimanus") == T)]
traits <- traits[-which(rownames(traits) %in% c("Thomomys", "Sorex", "Scapanus latimanus") == T),]
rowSums(SC_MNI, na.rm = T) #minimum of 27 individuals

data.sub <- SC_MNI

data.sub[is.na(data.sub)] <- 0

sub.rich <- NULL
sub.FD <- list()
#null.FD <- list()
abun <- list()
MIN <- 27
for(t in 1:500){
  sample.list <- list()
  
  for(i in 1:nrow(data.sub)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = MIN, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub)[str_detect(colnames(data.sub), 
                                                                                  paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3)
  #  for(i in 1:100){
  #    x2 <- apply(abun.mat, 1, sample)
  #    x2 <- t(x2)
  #    colnames(x2) <- colnames(abun.mat)
  #    x2 <- as.data.frame(x2)
  #    null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
  #    null.FD[[(t - 1)*100 + i]] <- null
  #  }
  print(t)
}

SC_MNIFD <- sub.FD

#saveRDS(SC_MNIFD, file = "SC_grid_MNI_FD.rds")

#### NISP ####
SC.NISP <- SC.counts
traits <- SC_traits

SC.NISP <- SC.NISP[,-which(colnames(SC.NISP) %in% c("Thomomys", "Sorex", "Scapanus latimanus") == T)]
traits <- traits[-which(rownames(traits) %in% c("Thomomys", "Sorex", "Scapanus latimanus") == T),]

rowSums(SC.NISP, na.rm = T) #minimum of 112 individuals

data.sub <- SC.NISP

sub.rich <- NULL
sub.FD <- list()
#null.FD <- list()
abun <- list()
MIN <- 112
for(t in 1:500){
  sample.list <- list()
  
  for(i in 1:nrow(data.sub)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = MIN, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub)[str_detect(colnames(data.sub), 
                                                                                  paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3)
  #  for(i in 1:100){
  #    x2 <- apply(abun.mat, 1, sample)
  #    x2 <- t(x2)
  #    colnames(x2) <- colnames(abun.mat)
  #    x2 <- as.data.frame(x2)
  #    null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
  #    null.FD[[(t - 1)*100 + i]] <- null
  #  }
  print(t)
}
SC_NISPFD <- sub.FD

#saveRDS(SC_NISPFD, file = "SC_NISP_FD.rds")

#### "...we subsampled SC to 45 individuals..." #### 
data <- SC_counts
time <- SC.time
traits <- SC_traits


data.sub <- data[which(rowSums(data) > 129),]
time.sub <- time[which(rowSums(data) > 129)]


sub.rich <- NULL
sub.FD <- list()
#null.FD <- list()
abun <- list()
MIN <- 45 # "...45 individuals..."
for(t in 1:500){ # "Subsampling was repeated 500 times."
  sample.list <- list()
  
  for(i in 1:nrow(data.sub50)){
    comm <- NULL
    for(j in 1:length(colnames(data.sub))){
      names <- rep(colnames(data.sub)[j], data.sub[i,j])
      comm <- c(names, comm)
    }
    sample.list[[i]] <- comm
  }
  
  abun.mat <- NULL
  
  for(i in 1:length(sample.list)){
    sub <- sample(sample.list[[i]], size = MIN, replace = F)
    sub <- as.data.frame(table(sub))
    sub$bin <- i
    abun.mat <- rbind(abun.mat, sub)
  }
  
  abun.mat <- abun.mat %>% 
    group_by(sub, bin) %>%
    spread(key = sub, value = Freq)
  abun.mat <- as.data.frame(abun.mat)
  abun.mat <- abun.mat[,-1]
  if(ncol(abun.mat != ncol(data.sub))){
    start <- ncol(abun.mat) + 1
    
    abun.mat[,start:ncol(data.sub)] <- NA
    colnames(abun.mat)[start:ncol(data.sub)] <- colnames(data.sub)[str_detect(colnames(data.sub), 
                                                                              paste(colnames(abun.mat), collapse = '|')) == F]
  }
  
  abun.mat <- abun.mat[,order(colnames(abun.mat))]
  abun.mat[is.na(abun.mat)] <- 0
  abun[[t]] <- abun.mat
  sub.rich <- cbind(sub.rich, specnumber(abun.mat))
  sub.FD[[t]] <- dbFD2(x = traits, a = abun.mat, w = Weights, CWM.type = "all", m = 3)
  #for(i in 1:100){ 
  #  x2 <- apply(abun.mat, 1, sample)
  #  x2 <- t(x2)
  #  colnames(x2) <- colnames(abun.mat)
  #  x2 <- as.data.frame(x2)
  #  null <- dbFD2(x = traits, w = Weights, a = x2, m = 3, CWM.type = "all")
  #  null.FD[[(t - 1)*100 + i]] <- null
  #}
  print(t)
}

SC_45FD <- sub.FD


#write_rds(SC_45FD, "SC_grid_FD.rds")


#### Ended of Script File ####



