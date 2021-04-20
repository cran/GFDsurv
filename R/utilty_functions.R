#' @export
print.casanova <- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "CASANOVA: Cumulative Aalen survival analyis-of-variance:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.casanova <- function (object, ...) {
  x <- object
  if ( length(x$rg) == 0 ){
    cat("The chosen weights are linearly independent.", "\n",
        "The test is based on the crossing weight.", "\n","\n")
  }else{
    rg_rep <- paste0("c(", x$rg[[1]][1], ", ", x$rg[[1]][2], ")" )
    if ( length(x$rg) > 1 ){
      for (i in 2:length(x$rg)){
        rg_rep <- paste0( rg_rep, ", c(", x$rg[[i]][1], ", ", x$rg[[i]][2], ")" )
      }
    }
    cat("The chosen weights are", if ( x$indep == FALSE){ " not"},
        " linearly independent.", "\n", "The test is based on ",
        if ( x$cross == TRUE){ "the crossing weight and "},
        length(x$rg), " ","weight", if (length(x$rg) > 1){"s"}, " with exponents ", "\n",
        "     ", "c(r,g)=  ", rg_rep, ".", "\n","\n", sep="")
  }
  print(x)
}


#' @export
plot.casanova <- function (x, direction = "horizontal", by.group = FALSE, nr.group = FALSE,...) {
  plotting <- x$plotting
  requireNamespace("survival", quietly = TRUE)

  if (!("package:survival" %in% search())) {
    attachNamespace("survival")
  }
  requireNamespace("survminer", quietly = TRUE)

  if (!("package:survminer" %in% search())) {
    attachNamespace("survminer")
  }
  requireNamespace("gridExtra", quietly = TRUE)

  if (!("package:gridExtra" %in% search())) {
    attachNamespace("gridExtra")
  }



  if(length(plotting$nadat2)==1){

    fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1],", data = plotting$dat)")))

    plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
               censor = FALSE,
               #ggtheme = theme_bw(),
               pval = FALSE)
    return(plot_1)


  }
  if(length(plotting$nadat2)==2){

    if(by.group  ==  FALSE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                           censor = FALSE,
                           #ggtheme = theme_bw(),
                           pval = FALSE,
                          # surv.median.line = "hv",
                           facet.by = plotting$nadat2[1])
      plot_2 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                          censor = FALSE,
                          #ggtheme = theme_bw(),
                          pval = FALSE,
                         # surv.median.line = "hv",
                          facet.by = plotting$nadat2[2])
      if(direction == "horizontal"){
       gridExtra::grid.arrange(plot_1, plot_2, ncol=2)
      }
      if(direction == "vertical"){
       gridExtra::grid.arrange(plot_1, plot_2, ncol=1)
      }
    }
    if(by.group == TRUE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      # surv.median.line = "hv",
                                      facet.by = plotting$nadat2[nr.group])
      return(plot_1)
    }
  }

  if(length(plotting$nadat2)==3){
    stop("Plots for designs with MORE than two factors not implemented!")
  }
}






#' @export
print.medsanova<- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)

  cat("\n", "medSANOVA: Median survival analyis-of-variance:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.medsanova <- function (object, ...) {
  x <- object
  print(x)
}


#' @export
plot.medsanova  <- function (x, direction = "horizontal", by.group = FALSE, nr.group = FALSE,...) {
  plotting <- x$plotting
  requireNamespace("survival", quietly = TRUE)

  if (!("package:survival" %in% search())) {
    attachNamespace("survival")
  }
  requireNamespace("survminer", quietly = TRUE)

  if (!("package:survminer" %in% search())) {
    attachNamespace("survminer")
  }
  requireNamespace("gridExtra", quietly = TRUE)

  if (!("package:gridExtra" %in% search())) {
    attachNamespace("gridExtra")
  }




  if(length(plotting$nadat2)==1){

    fit <- eval(parse(text =paste0("survival::survfit(Surv(Var,event) ~ ",plotting$nadat2[1],", data = plotting$dat)")))

    plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                         censor = FALSE,
                         #ggtheme = theme_bw(),
                         pval = FALSE,
                         surv.median.line = "hv")
    return(plot_1)

  }

  if(length(plotting$nadat2)==2){
    if(by.group  ==  FALSE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(Var,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      surv.median.line = "hv",
                                      facet.by = plotting$nadat2[1])
      plot_2 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      surv.median.line = "hv",
                                      facet.by = plotting$nadat2[2])
      if(direction == "horizontal"){
        gridExtra::grid.arrange(plot_1, plot_2, ncol=2)
      }
      if(direction == "vertical"){
        gridExtra::grid.arrange(plot_1, plot_2, ncol=1)
      }
    }
    if(by.group == TRUE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(Var,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      surv.median.line = "hv",
                                      facet.by = plotting$nadat2[nr.group])
      return(plot_1)
    }


  }

  if(length(plotting$nadat2)==3){
    stop("Plots for designs with MORE than two factors not implemented!")
  }
}


#' @export
print.copsanova <- function(x, ...) {
  cat("Call:", "\n")
  print(x$input$formula)
  cat("\n", "multiplier bootstrap:",x$weights, "\n")

  cat("\n", "time window for analysis: [0,",paste0(x$tau,"]"), "\n")

  cat("\n", "copSANOVA: concordance probability SANOVA:","\n","\n", sep = "")
  print(x$statistic)
}

#' @export
summary.copsanova <- function (object, ...) {
  x <- object
  print(x)
}

#' @export
plot.copsanova  <- function (x, direction = "horizontal", by.group = FALSE, nr.group = FALSE,...) {
  plotting <- x$plotting
  requireNamespace("survival", quietly = TRUE)

  if (!("package:survival" %in% search())) {
    attachNamespace("survival")
  }
  requireNamespace("survminer", quietly = TRUE)

  if (!("package:survminer" %in% search())) {
    attachNamespace("survminer")
  }
  requireNamespace("gridExtra", quietly = TRUE)

  if (!("package:gridExtra" %in% search())) {
    attachNamespace("gridExtra")
  }




  if(length(plotting$nadat2)==1){

    fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1],", data = plotting$dat)")))

    plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                         censor = FALSE,
                         #ggtheme = theme_bw(),
                         pval = FALSE)
    return(plot_1)


  }
  if(length(plotting$nadat2)==2){

    if(by.group  ==  FALSE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      # surv.median.line = "hv",
                                      facet.by = plotting$nadat2[1])
      plot_2 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      # surv.median.line = "hv",
                                      facet.by = plotting$nadat2[2])
      if(direction == "horizontal"){
        gridExtra::grid.arrange(plot_1, plot_2, ncol=2)
      }
      if(direction == "vertical"){
        gridExtra::grid.arrange(plot_1, plot_2, ncol=1)
      }
    }
    if(by.group == TRUE){

      fit <- eval(parse(text =paste0("survival::survfit(Surv(time,event) ~ ",plotting$nadat2[1]," + ",plotting$nadat2[2],", data = plotting$dat)")))

      plot_1 <- survminer::ggsurvplot(fit, data = plotting$dat, fun = "pct",
                                      censor = FALSE,
                                      #ggtheme = theme_bw(),
                                      pval = FALSE,
                                      # surv.median.line = "hv",
                                      facet.by = plotting$nadat2[nr.group])
      return(plot_1)
    }
  }

  if(length(plotting$nadat2)==3){
    stop("Plots for designs with MORE than two factors not implemented!")
  }
}

