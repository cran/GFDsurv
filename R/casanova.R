
#' CASANOVA: Cumulative Aalen survival analyis-of-variance
#'
#' The function \code{casanova} calculates the Wald-type statistic based on the
#' combination of differently weighted Nelson-Aalen-type integrals. Respective p-values
#' are obtained by a \eqn{\chi^2}-approximation and a permutation approach, respectively.
#' @param formula A model \code{formula} object. The left hand side contains the time variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param event The name of censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is "event"
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param cross logical. Should the crossing weight w(x) = 1 - 2x be included?
#'  The default is \code{TRUE}.
#' @param rg A list containing the exponents \code{c(r, g)} of the
#' weights
#'   \eqn{w(x) = x^r (1-x)^g}. Both exponents need to be natural numbers including 0.
#'  Default is \code{list( c(0, 0) )} corresponding to the log-rank weight.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#'
#' @details The \code{casanova} function calculates the Wald-type statistic of weighted
#' Nelson-Aalen type integrals
#' for general factorial survival designs. Crossed as well as hierachically nested designs are
#' implemented. Moreover, the approach allows the combination of
#' different weights into a
#' joint statistic. The user can choose between weights of the following form:
#' w(x) = 1 - 2x (\code{cross = TRUE}) and w(x) = x^r * (1-x)^g for natural numbers
#' r,g (including 0). The function automatically check whether the specified weights
#' fulfill the linear independence assumption and choose a subset of linearly independent
#' weights if the original weights violate the aforemention assumption.
#'
#'   The \code{casanova} function returns the test statistic as well as two
#'   corresponding p-values: the first is based on a \eqn{\chi^2} approximation and
#'   the second one is based on a permutation procedure.
#'
#' @return  A \code{casanova} object containing the following components:
#'  \item{pvalues_stat}{The p-values obtained by \eqn{\chi^2}-approximation}
#'  \item{pvalues_per}{The p-values of the permutation approach}
#'  \item{statistics}{The value of the casanova along with degrees of freedom of the
#'  central chi-square distribution and p-value, as well as the p-value of the
#'   permutation procedure.}
#'  \item{rg}{A list containing the exponents of the direction considered in the
#'  statistical analysis}
#'  \item{cross}{logical. Was the crossing direction considered in the statistical
#'  analysis}
#'  \item{indep}{logical. Were the directions specified by the user linearly
#'  independent?}
#'  \item{nperm}{The number of permutations used for calculating the permuted p-value.}
#'
#' @examples
#' \donttest{
#' library("survival")
#' data(veteran)
#' out <- casanova(formula ="time ~ trt*celltype",event = "status",
#'  data = veteran)
#'
#' ## Detailed informations:
#' summary(out)
#'}
#' @references Ditzhaus, M., Janssen, A. and Pauly, M. (2020). Permutation inference in factorial survival designs with the
#'           CASANOVA. ArXiv preprint (arXiv:2004.10818v2).
#'
#' @import stats
#' @importFrom magic adiag
#' @importFrom  MASS ginv
#' @importFrom survival survfit
#' @importFrom survminer ggsurvplot
#' @importFrom gridExtra grid.arrange


#' @export
#'
casanova <- function(formula, event ="event", data = NULL, nperm = 1999,
                      cross = TRUE, nested.levels.unique = FALSE, rg = list(c(0,0))){
  input_list <- list(formula = formula,event = event, data = data, nperm = nperm)
  #Zeit und in Formel einbinden
  formula2 <-  paste0(formula,"*",event)
  dat <- model.frame(formula2, data)


  #n
  subject <- 1:nrow(dat)
  n_all <- length(subject)

  formula <- as.formula(formula)
  nf <- ncol(dat) - 1 - 1
  nadat <- names(dat)

  if(anyNA(dat[,nadat])){
    stop("Data contains NAs!")
  }


  names(dat) <- c("time",nadat[2:(1+nf)],"event")

  dat2 <- data.frame(dat, subject = subject)

  ###Daten richtig ordnen

  nadat2 <- nadat[-c(1,nf+2)]
  dat2[,"time"] <- dat2[,"time"]  + runif(length(dat2[,"time"])) * 10^-7

  weight_names <- c("Combination")
  #Koeffizientencheck
            if (is.null(rg) == FALSE) {
              out <- coeff.check(cross = cross, rg = rg)
              cross <- out$cross
              indep <- out$indep
              rg <- out$rg
              weight_names <- c("Combination", rep(0,length(rg)))
              for(i in 1:length(rg)){
                weight_names[i+1] <- paste0("x^",rg[[i]][1],"(1-x)^",rg[[i]][2])
              }
            } else {
              indep <- TRUE
            }
            w <- list()

            w.funct <- function(x) {
              x <- unlist(x)
              w <- function(y) {
                y^x[1] * (1 - y)^x[2]
              }
            }
            match.fun(w.funct)
            w <- lapply(rg, w.funct)
            if (cross == TRUE) {
              w <- c(w, function(x) {
                1 - 2 * x
              })
            }

            if(cross == TRUE){
              weight_names <- c(weight_names,"1-2x")
            }
            ### Für prob
            if(sum(unlist(lapply(rg,function(x) sum(x == c(0,0))==2))==1)==1){
              weight_names[1+which(unlist(lapply(rg,function(x) sum(x == c(0,0))==2))==1)] <-"prop"
            }

  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)
  if (nf == 1) {
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    hypo_matrices <- list(diag(fl) - matrix(1/fl, ncol = fl, nrow = fl))
    group <- rep(1:length(n),n)
    dat2$group <- group

    ###############################
    dat2  <- dat2[order(dat2$time),]
    event <- dat2[,"event"]
    group <- dat2$group
    #print(hypo_matrices)
    results <- stat_factorial(hypo_matrices,group, event,n, n_all, w, nperm)


    #Ergebnis in Tabellenform
    m <- length(w)
    Stat_Erg <- results$Stat
    Stat_Erg  <- matrix(unlist(Stat_Erg),length(hypo_matrices),
                        m+1,byrow = TRUE)
    rank_C <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))
    # #Quantilmatrix
    # q_uncon <- sapply(1:length(hypo_matrices),function(x) qchisq(1-alpha, df = rank_C[x]))
    # q_uncon_c <- sapply(1:length(hypo_matrices),function(x) qchisq(1-alpha, df = rank_C[x]*m))
    # q_uncon <- matrix(c(q_uncon_c,rep(q_uncon,m)),length(hypo_matrices), m+1)
    #Tabelle der P-Werte
    pvalue_stat <-  round(t(sapply(1:length(hypo_matrices), function(x) c(1-pchisq(Stat_Erg[x,1],df=rank_C[x]*m),
                                                                          1-pchisq(Stat_Erg[x,2:(m+1)],df=rank_C[x])))),3)
    pvalue_stat <- matrix(unlist(pvalue_stat),length(hypo_matrices),m+1,
                          dimnames = list(fac_names, weight_names))


    #P-Werte für Perm
    per <- results$Perm
    per_unlist <-  matrix(unlist(lapply(per, t)),length(hypo_matrices)*(m+1),byrow = TRUE)
    stat_Erg_unlist <- unlist(Stat_Erg)
    pvalue_per <- sapply(1:(length(hypo_matrices)*(m+1)),function(x) mean(per_unlist[x,]>stat_Erg_unlist[x]))
    pvalue_per <- matrix(pvalue_per ,length(hypo_matrices),byrow = TRUE)
    pvalue_per <- matrix(unlist(pvalue_per),length(hypo_matrices),(m+1),dimnames = list(fac_names, weight_names))
    df <- sapply(1:length(hypo_matrices),function(x)rank_C[x]*m)

  }
  else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]),
                            ]
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    ###

    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure

    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is\n           not impemented!")
    }
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are\n           not implemented!")
    }
    if (length(fac_names) == nf) {
      TYPE <- "nested"
      if (nested.levels.unique) {
        n <- n[n != 0]
        blev <- list()
        lev_names <- list()
        for (ii in 1:length(levels[[1]])) {
          blev[[ii]] <- levels(as.factor(dat[, 3][dat[,
                                                      2] == levels[[1]][ii]]))
          lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
        }
        if (nf == 2) {
          lev_names <- as.factor(unlist(lev_names))
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev)
        }
        else {
          lev_names <- lapply(lev_names, rep, length(levels[[3]])/length(levels[[2]]))
          lev_names <- lapply(lev_names, sort)
          lev_names <- as.factor(unlist(lev_names))
          blev <- lapply(blev, rep, length(levels[[3]])/length(levels[[2]]))
          blev <- lapply(blev, sort)
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev,
                                        as.factor(levels[[3]]))
        }
        if (nf == 2) {
          fl[2] <- fl[2]/fl[1]
        }
        else if (nf == 3) {
          fl[3] <- fl[3]/fl[2]
          fl[2] <- fl[2]/fl[1]
        }
      }
      hypo_matrices <- HN(fl)
    }
    else {
      TYPE <- "crossed"
      hypo_matrices <- HC(fl, perm_names, fac_names)[[1]]
      fac_names <- HC(fl, perm_names, fac_names)[[2]]
    }
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    if (TYPE == "nested" & 0 %in% n & nested.levels.unique ==
        FALSE) {
      stop("The levels of the nested factor are probably labeled uniquely,\n           but nested.levels.unique is not set to TRUE.")
    }
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination\n           with less than 2 observations!")
    }
    group <- rep(1:length(n),n)
    dat2$group <- group
    #print(n)
    ###############################
    dat2  <- dat2[order(dat2$time),]
    event <- dat2[,"event"]
    group <- dat2$group

    #print(hypo_matrices)
    results <- stat_factorial(hypo_matrices,group, event,n, n_all, w, nperm)


    #Ergebnis in Tabellenform
    m <- length(w)
    Stat_Erg <- results$Stat
    Stat_Erg  <- matrix(unlist(Stat_Erg),length(hypo_matrices),
                        m+1,byrow = TRUE)
    rank_C <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))
    #Quantilmatrix
    # q_uncon <- sapply(1:length(hypo_matrices),function(x) qchisq(1-alpha, df = rank_C[x]))
    # q_uncon_c <- sapply(1:length(hypo_matrices),function(x) qchisq(1-alpha, df = rank_C[x]*m))
    # q_uncon <- matrix(c(q_uncon_c,rep(q_uncon,m)),length(hypo_matrices), m+1)
    #Tabelle der P-Werte
    pvalue_stat <-  round(t(sapply(1:length(hypo_matrices), function(x) c(1-pchisq(Stat_Erg[x,1],df=rank_C[x]*m),
                                                          1-pchisq(Stat_Erg[x,2:(m+1)],df=rank_C[x])))),4)
    pvalue_stat <- matrix(unlist(pvalue_stat),length(hypo_matrices),m+1,
                          dimnames = list(fac_names, weight_names))


    #P-Werte für Perm
    per <- results$Perm
    per_unlist <-  matrix(unlist(lapply(per, t)),length(hypo_matrices)*(m+1),byrow = TRUE)
    stat_Erg_unlist <- unlist(t(Stat_Erg))
    pvalue_per <- sapply(1:(length(hypo_matrices)*(m+1)),function(x) mean(per_unlist[x,]>stat_Erg_unlist[x]))
    pvalue_per <- matrix(pvalue_per ,length(hypo_matrices),byrow = TRUE)
    pvalue_per <- matrix(unlist(pvalue_per),length(hypo_matrices),(m+1),dimnames = list(fac_names, weight_names))
    df <- sapply(1:length(hypo_matrices),function(x)rank_C[x]*m)

  }
  output <- list()
  output$input <- input_list
  output$pvalues_stat  <- pvalue_stat
  output$pvalues_per <-  pvalue_per
  output$cross <- cross
  output$indep <- indep
  output$rg <- rg
  output$nperm <-nperm
  output$plotting <- list("dat" = dat,"nadat2" = nadat2)

  output$statistic <- cbind(Stat_Erg[,1],df,pvalue_stat[,1],round(pvalue_per[,1],4))
  rownames(output$statistic) <- fac_names
  colnames(output$statistic) <- c("Test statistic","df","p-value", "p-value perm")


   class(output) <- "casanova"
    return(output)


}

