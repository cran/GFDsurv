#' medSANOVA: Median survival analysis-of-variance
#'
#' The function \code{medsanova} calculates the Wald-type test statistic for
#' inferring median survival differences in general factorial designs.
#' Respective p-values are obtain by a \eqn{\chi^2}-approximation and a permutation approach.
#' @param formula A model \code{formula} object. The left hand side contains the time variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param event The name of the censoring status indicator with values 0=censored and
#' 1=uncensored.
#' The default choice is \code{"event"}.
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param nonex_action One of the options \code{"redraw"} and \code{"setInf"}, 
#' which specifies how permutation samples with non-existing medians should be 
#' handled, see details. The default option is \code{"setInf"}.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is \code{FALSE}, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#' @param var_method Method for the variance estimation of the sample medians. The default
#'  is the "one-sided" confidence interval approach. Additionally, the "two-sided" confidence
#'  interval approach can be used.
#' @param var_level A number between 0 and 1 specifying the confidence level for the
#'  variance estimation method; the default value is \code{0.9}.
#' @param seed A single value, interpreted as an integer, or \code{NULL}; used for the 
#'  seed to guarantee reproducibility. The default value is \code{1}.
#'
#' @details
#' The \code{medsanova} function calculates the Wald-type statistic for median differences
#' in general factorial survival designs. Crossed as well as hierarchically nested designs are
#' implemented. To estimate the sample medians' variances, a one-sided (resp. two-sided) confidence
#' interval approach is used and the level of this confidence interval can be specified by \code{var_level}.
#'
#'   The \code{medsanova} function returns the test statistic as well as two
#'   corresponding p-values: the first is based on a \eqn{\chi^2} approximation and
#'   the second one is based on a permutation procedure.
#'   
#'    In the two-sample case, the \code{medsanova} function also provides the estimate for the
#'    difference between the medians, (asymptotic 
#'    and permutation based) 95\% confidence intervals for the difference as well as the standard error SE 
#'    of the median difference.
#'   
#'   For the argument \code{nonex_action}, \code{"redraw"} means that a new permutation 
#'   is drawn instead while \code{"setInf"} means that the permutation statistic is set 
#'   to \code{Inf}. The first option might lead to an incorrect type-I error 
#'   even under exchangeability while the latter yields conservative test decisions.
#'   
#' @return An  \code{medsanova} object contains the following components:
#'  \item{pvalues_stat}{The p-values obtained by \eqn{\chi^2}-approximation}
#'  \item{pvalues_per}{The p-values of the permutation approach}
#'  \item{statistics}{The value of the Wald-type test statistic along with the
#'  degrees of freedom of the \eqn{\chi^2}-distribution and the
#'  respective p-value, as well as the p-value of the
#'   permutation procedure.}
#'  \item{nperm}{The number of permutations used for calculating the permuted p-value.}
#'  \item{medians}{The calculated survival medians in all subgroups.}
#' @examples
#' \donttest{
#' library("survival")
#' data(veteran)
#' out <- medsanova(formula ="time ~ trt*celltype",event = "status",
#'  data = veteran)
#'
#' ## Detailed informations:
#' summary(out)
#' ## show the medians
#' out$medians
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

medsanova <-  function(formula, event ="event", data = NULL, nperm = 1999, nonex_action = "setInf",
                       var_method = "twosided", var_level = 0.9,
                       nested.levels.unique = FALSE, seed = 1){
  set.seed(seed)
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
  
  twosample <- FALSE
  
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
                         var_method = var_method, var_level = var_level, nonex_action = nonex_action)
    
    if(length(hypo_matrices) == 1 & identical(hypo_matrices[[1]],matrix(c(0.5,-0.5,-0.5,0.5),ncol=2))){
      twosample <- TRUE
    }
    
    for(j in 1:length(hypo_matrices)){
      q_perm <- erg_perm$test_stat_erg
      t_int_perm <- mean(erg_stat[paste0("int_", j)] <= q_perm[paste0("int_", j), ], na.rm = TRUE)
      t_int_chi <- 1-pchisq(erg_stat[paste0("int_", j)], df = qr(hypo_matrices[[j]])$rank )
      
      t_int_perm <- ifelse(is.nan(t_int_perm), NA, t_int_perm)
      
      out1 <- c("perm" = t_int_perm, "chi" = t_int_chi)
      
      if(twosample){
        t_mat_list <- list(hypo_matrices[[1]])
        
        values <- sort_data(dat3)
        values[, 3] <- group
        values_KME <- KME(values, group = values[, 3])
        
        var_int <- int_var_groups(values_KME,var_level = var_level, group = values[ ,3], var_method = var_method)
        
        est <- t(c(1,-1))%*%var_int$Median
        SE <- sqrt(sum(var_int$Variance))/sqrt(nrow(values))
        perm_quan <- sqrt(quantile(q_perm, prob = 0.95))
        chi_quan <- qnorm(0.975)
        out2 <- c("Estimate" = est, "Std. Error" = SE, "perm_lower" = est - SE*perm_quan, "perm_upper" = est + SE*perm_quan, 
                  "chi_lower" = est - SE*chi_quan, "chi_upper" = est + SE*chi_quan)
      }
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
                         var_method = var_method, var_level = var_level, nonex_action = nonex_action)
    
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
  
  if(twosample){
    output$difference <- round(out2[c(1,2,5,6,3,4)],3)
    names(output$difference) <- c("Diff. estimate","Std. Error","lower 95%","upper 95%",
                                    "perm lower 95%","perm upper 95%")
    
  }
  # output the medians
  df <- cbind(dat2, group_med = KME(dat3, group = dat3[, 3])[,5])
  df_u <- unique(df[, c(2:(nf + 1), which(names(df) == "group_med"))])
  nm <- apply(df_u[, 1:nf,drop=FALSE], 1, function(row) {
    paste0(
      names(df_u)[1:nf],
      ": ",
      row,
      collapse = "; "
    )
  })
  output$medians <- setNames(df_u$group_med, nm)[order(nm)]
  
  class(output) <- "medsanova"
  return(output)
  
}
