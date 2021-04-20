#' medSANOVA: Median survival analyis-of-variance
#'
#' The function \code{medsanova} calculates the Wald-type test statistic for
#' inferring median survival differences in general factorial designs.
#' Respective p-values are obtain by a \eqn{\chi^2}-approximation and a permutation approach.
#' @param formula A model \code{formula} object. The left hand side contains the time variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param event The name of the censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is "event"
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#' @param var_method Method for the variance estimation of the sample medians. The default
#'  is the "one-sided" confidence interval approach. Additionally, the "two-sided" confidence
#'  interval approach can be used.
#' @param var_level A number between 0 and 1 specifying the confidence level for the
#'  variance estimation method; the default value is 0.9.
#'
#' @details
#' The \code{medsanova} function calculates the Wald-type statistic for median differences
#' in general factorial survival designs. Crossed as well as hierachically nested designs are
#' implemented. To estimate the sample medians' variances, a one-sided (resp. two-sided) confidence
#' interval approach is used and the level of this confidence interval can be specified by \code{var_level}.
#'
#'   The \code{medsanova} function returns the test statistic as well as two
#'   corresponding p-values: the first is based on a \eqn{\chi^2} approximation and
#'   the second one is based on a permutation procedure.
#' @return An  \code{medsanova} object containing the following components:
#'  \item{pvalues_stat}{The p-values obtained by \eqn{\chi^2}-approximation}
#'  \item{pvalues_per}{The p-values of the permutation approach}
#'  \item{statistics}{The value of the Wald-type test statistic along with the
#'  degrees of freedom of the \eqn{\chi^2}-distribution and the
#'  respective p-value, as well as the p-value of the
#'   permutation procedure.}
#'  \item{nperm}{The number of permutations used for calculating the permuted p-value.}
#' @examples
#' \donttest{
#' library("survival")
#' data(veteran)
#' out <- medsanova(formula ="time ~ trt*celltype",event = "status",
#'  data = veteran)
#'
#' ## Detailed informations:
#' summary(out)
#' }
#' @references Ditzhaus, M., Dobler, D. and Pauly, M.(2020). Inferring median survival
#'  differences in general factorial designs via permutation tests.
#'  Statistical Methods in Medical Research. doi:10.1177/0962280220980784.
#'
#' @importFrom  MASS ginv
#' @importFrom survival survfit
#' @importFrom survminer ggsurvplot
#' @importFrom gridExtra grid.arrange
#' @import plyr
#' @export
#'
#'

medsanova <-  function(formula, event ="event", data = NULL, nperm = 1999,
                       var_method = "twosided", var_level = 0.9,
                       nested.levels.unique = FALSE){
  input_list <- list(formula = formula, event ="event", data = data, nperm = nperm,
                     var_level = var_level)
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

  if(var_method == "twosided"){var_method = 2}
  if(var_method == "onesided"){var_method = 3}

  names(dat) <- c("Var",nadat[2:(1+nf)],"event")

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

    ###############################
    dat2  <- dat2[order(dat2$Var),]
    event <- dat2[,"event"]
    group <- dat2$group

    dat3 <- dat2[,c("Var","event","group")]


    erg_stat <-  wrap_sim2(dat3,group = dat3[,3],hypo_matrices,
                           var_method = var_method, var_level = var_level)
    out <- list()

    erg_perm <- perm_fun(dat3, nperm, hypo_matrices,
                         var_method = var_method, var_level = var_level)

    for(j in 1:length(hypo_matrices)){
      q_perm <- erg_perm$test_stat_erg
      t_int_perm <- mean(erg_stat[paste0("int_", j)] <= q_perm[paste0("int_", j), ], na.rm = TRUE)
      t_int_chi <- 1-pchisq(erg_stat[paste0("int_", j)], df = qr(hypo_matrices[[j]])$rank )

      t_int_perm <- ifelse(is.nan(t_int_perm), NA, t_int_perm)

      out1 <- c("perm" = t_int_perm, "chi" = t_int_chi)
      out[[j]] <- out1
    }

    out <- matrix(unlist(out),length(hypo_matrices),byrow=T)

    df <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))

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
    dat2  <- dat2[order(dat2$Var),]
    event <- dat2[,"event"]
    group <- dat2$group

    dat3 <- dat2[,c("Var","event","group")]


    erg_stat <-  wrap_sim2(dat3,group = dat3[,3],hypo_matrices,
                           var_method = var_method, var_level = var_level)
    out <- list()

    erg_perm <- perm_fun(dat3, nperm, hypo_matrices,
                         var_method = var_method, var_level = var_level)

    for(j in 1:length(hypo_matrices)){
      q_perm <- erg_perm$test_stat_erg
      t_int_perm <- mean(erg_stat[paste0("int_", j)] <= q_perm[paste0("int_", j), ], na.rm = TRUE)
      t_int_chi <- 1-pchisq(erg_stat[paste0("int_", j)], df = qr(hypo_matrices[[j]])$rank )

      t_int_perm <- ifelse(is.nan(t_int_perm), NA, t_int_perm)

      out1 <- c("perm" = t_int_perm, "chi" = t_int_chi)
      out[[j]] <- out1
    }

    out <- matrix(unlist(out),length(hypo_matrices),byrow=T)

    df <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))

  }


   output <- list()
   output$input <- input_list
   output$nperm <-nperm
   output$plotting <- list("dat" = dat,"nadat2" = nadat2)


   output$statistic <- cbind(erg_stat,df,round(out[,2],3),round(out[,1],3))
   rownames(output$statistic) <- fac_names
   colnames(output$statistic) <- c("Test statistic","df","p-value", "p-value perm")
   class(output) <- "medsanova"
   return(output)

}
