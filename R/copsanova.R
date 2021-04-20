#' copSANOVA: concordance parameter survival analyis-of-variance
#'
#' The function \code{copanova} calculates the ANOVA-rank-type statistic for general
#' factorial
#' survival designs based on the (extended) concordance parameter. The respective
#'p-value is
#' obtained by a multiplier bootstrap approach.
#' @param formula A model \code{formula} object. The left hand side contains the
#' time variable
#' and the right hand side contains the factor variables of interest. An interaction
#' term must be
#'  specified.
#' @param event The name of censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is "event"
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param BSiter The number of bootstrap iterations; the default is 1999.
#' @param weights Character to specify the multiplier bootstrap approach. Either a
#' wild bootstrap
#' with centred Poisson ("pois", default) or standard normal ("norm") weights, or the
#' weird bootstrap ("weird") can be chosen. Moreover, both wild bootstrap strategies
#' can be selected
#' with a correcting factor for liberality by "corrLibPois" and "corrLibNorm".
#' @param tau The truncation time specifying the end of the relevant time window for
#' the analysis.
#' By default (\code{NULL}), the smallest 95\%-quantile of the times per group is
#' chosen.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#' @details
#' The \code{copsanova} function calculates the ANOVA-rank-type statistic for
#' general factorial
#' survival designs based on the (extended) concordance parameter. Crossed as well as
#' hierachically nested designs are implemented. The p-value is determined by a
#' multiplier bootstrap
#' approach. Here, a wild bootstrap with/without correcting factors for liberal
#' tests or the weird
#' bootstrap of Andersen et al. (1993) can be chosen. The concrete analysis is done
#' on the time window
#' [0,tau], where tau need to be chosen equal to (default) or smaller than the
#' smallest out of
#' the largest possible censoring times per group.
#'
#'   The \code{copsanova} function returns the test statistic as well as a
#'   corresponding p-value based on a the specified multiplier procedure.
#'
#' @return An  \code{copsanova} object containing the following components:
#'  \item{statistics}{The value of the copsanova along with the p-value
#'  of the specified multiplier bootstrap.}
#'  \item{Bsiter}{The number of bootstrap iterations.}
#'  \item{weights}{The chosen multiplier bootstrap method.}
#'  \item{tau}{The chosen truncation time specifying the end of the relevant time window for
#' the analysis.}
#'
#' @examples
#' \donttest{
#' library(condSURV)
#' data(colonCS)
#' out <- copsanova(formula ="Stime ~ rx*sex",event = "event",
#'                  data = colonCS, BSiter = 99)
#'
#' ##Detailed informations:
#'summary(out)}
#'
#' @references
#'  Dobler, D. and Pauly, M. (2020). Factorial analyses of treatment effects
#'  under independent right-censoring. Statistical Methods in Medical Research 29(2), 325-343. doi:10.1177/0962280219831316.
#'
#' @importFrom survival survfit
#' @importFrom survminer ggsurvplot
#' @importFrom gridExtra grid.arrange
#'
#' @export
#'

copsanova <- function(formula, event ="event", data = NULL, BSiter = 1999,
                      weights = "pois", tau = NULL,
                     nested.levels.unique = FALSE){
  input_list <- list(formula = formula,time = time, data = data, BSiter = BSiter,
                     weights = weights)
  #Zeit und in Formel einbinden
  formula2 <-  paste0(formula,"*",event)
  dat <- model.frame(formula2, data)
  #n
  subject <- 1:nrow(dat)
  n_all <- length(subject)

  formula <- as.formula(formula)
  nf <- ncol(dat) - 1 - 1
  nadat <- names(dat)

  if(anyNA(data[,nadat])){
    stop("Data contains NAs!")
  }

  names(dat) <- c("time",nadat[2:(1+nf)],"event")

  dat2 <- data.frame(dat, subject = subject)

  nadat2 <- nadat[-c(1,nf+2)]


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

    ##########

    event <- dat2[,"event"]
    group <- dat2$group
    diff_groups <- length(unique(group))

    dat2$exit <- dat2$time
    dat2$to <- ifelse(dat2$event == 1,"cens",1)
        data1 <- list()

    n <- numeric(0)
    tau_all <- numeric(0)
    for( k in 1:diff_groups){
      dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
      n <- c(n,length(dat_tmp$to))

      tau_all <- c(tau_all,quantile(dat_tmp$exit,0.95))
      # ind_tau <- dat_tmp$exit >= tau
      #
      # dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      # dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

    if(is.null(tau)){
      tau <- min(tau_all)
    }

    for( k in 1:diff_groups){
      dat_tmp <- data1[[k]]
      ind_tau <- dat_tmp$exit >= tau

      dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

    copsanova_erg <- test.data(data1, n, BSiter = BSiter,
                               Gewichte = weights, c.matrix = hypo_matrices)


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
    group <- rep(1:length(n),n)
    dat2$group <- group
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

    ###############################
    event <- dat2[,"event"]
    group <- dat2$group
    diff_groups <- length(unique(group))

    dat2$exit <- dat2$time
    dat2$to <- ifelse(dat2$event == 1,"cens",1)
    data1 <- list()


    n <- numeric(0)
    tau_all <- numeric(0)
      for( k in 1:diff_groups){
        dat_tmp <- dat2[dat2$group == k, c("exit","to","group")]
        n <- c(n,length(dat_tmp$to))

        tau_all <- c(tau_all,quantile(dat_tmp$exit,0.95))
        # ind_tau <- dat_tmp$exit >= tau
        #
        # dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
        # dat_tmp$exit[ind_tau] <- tau
         data1[[k]] <- dat_tmp
      }


    if(is.null(tau) || is.na(tau)){
      tau <- min(tau_all)

    }

        for( k in 1:diff_groups){
      dat_tmp <- data1[[k]]
      ind_tau <- dat_tmp$exit >= tau

      dat_tmp$to[ind_tau] <- "1"  # Alles was gr??er oder gleich tau ist, wird als unzensiert angesetzt, damit der Kaplan-Meier-Sch?tzer in tau auf Null f?llt.
      dat_tmp$exit[ind_tau] <- tau
      data1[[k]] <- dat_tmp
    }

      copsanova_erg <- test.data(data1, n, BSiter = BSiter,
                      Gewichte = weights, c.matrix = hypo_matrices)


  }
  output <- list()
  output$input <- input_list
  output$bsiter <- BSiter
  output$weights <- weights
  output$tau <- tau
  output$plotting <- list("dat" = dat,"nadat2" = nadat2)


  output$statistic <- cbind(copsanova_erg$test_statistics,round(copsanova_erg$value,4))
  rownames(output$statistic) <- fac_names
  colnames(output$statistic) <- c("Test statistic","p-value")

  class(output) <- "copsanova"
  return(output)


}
