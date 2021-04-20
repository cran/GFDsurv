#Sort object of data_gen
#Input:
# values:   matrix - created by data_gen

sort_data <- function(values) values[order(values[, 1]), ]

#create Matrix C, containing the groups
#Input:
# values:   matrix - created by data_gen, possibly sorted by sort_data
#Output:
# Matrix in which a column represents a group and a row an individual
create_c_mat <- function(values) {
  group <- values[, 3]
  c_mat <- matrix(0, nrow = max(group), ncol = nrow(values))

  c_mat[(0:(nrow(values) - 1)) * 3 + group] <- 1
  return(t(c_mat))
}

#Example:
#create_c_mat(sort_data(test_dat))




#Kaplan-Meier estimator and median
#Input:
# values:   matrix created by data_gen and sorted by sort_data
# group:    integer vector containing the group of the observations. Default is
#           the third column of the values, the groups drawn by data_gen

#Output:
# matrix values with additional column containing the KME and the median
# for the corresponding group
KME <- function(values, group = values[, 3]) {
  #Add columns
  values <- cbind(values, 0, 0)
  colnames(values)[c(4, 5)] <- c("KME_VF", "group_med")

  for (i in 1:max(group)) {
    #calculate in group
    values2 <- values[group == i, ]
    n <- nrow(values2)
    vec_temp <- rep(1, n)

    #calculate the KME
    res <- 1 - cumprod(1 - values2[, 2] / (n + 1 - 1:n))

    #find median
    med_ind <- findInterval(0.5-1e-8, res) + 1
    if(med_ind > n) {
      #if median lies not in observations
      med <- Inf
    } else {
      med <- values2[med_ind, 1]
    }

    #add values to matri
    values[group == i, 4] <- res
    values[group == i, 5] <- rep(med, n)
  }
  return(values)
}

#Example
#test_KME <- KME(sort_data(test_dat))


#estimation of standard deviation by formula derived from confidence interval
#two-sided

#Help function to calculate the standard deviation out of the interval. Can
#only be used for one group
#Input:
# values:   - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Contains only one group
# n_all     - integer, total number of observations
# a         - numeric value between 0 and 1, alpha level, default 0.1
#Output:    Standard deviation estimated with the formula derived from interval

int_sd2 <- function(values, var_level  = var_level , n_all) {
  a <- 1 - var_level
  n <- nrow(values)
  indi <- values[, 1] <= values[, 5]

  Vi <- sum(indi * values[, 2]/(n:1)^2)
  z <- qnorm(1 - a/2)

  li <- max(0, 1/2*(1 - z * sqrt(Vi)))
  #lower quantile
  ui <- min(1, 1/2*(1 + z * sqrt(Vi)))
  #upper quantile

  if(ui > max(values[, 4])) {
    ui <- max(values[, 4])
    a <- 2 * dnorm((2*ui - 1)/sqrt(Vi/n))
    z <- qnorm(1-a/2)
    li <- max(0, 1/2*(1 - z * sqrt(Vi)))

    if(is.infinite(values[1, 5])){
      return(c("sigma" = NA))
    }
  }
  #upper quantile does not work

  fui_ind <- findInterval(ui, values[, 4])
  #index upper quantile, if the upper Interval is not in n, it is the last value
  fli_ind <- findInterval(li, values[, 4])
  #index lower quantile

  sig <- 1/(2*z) * sqrt(n_all) *
    (values[min(fui_ind+1, n), 1] - values[fli_ind+1, 1])

  return(c("sigma" = unname(sig)))
}





#estimation of standard deviation by formula derived from confidence interval
#one-sided

#Help function to calculate the standard deviation out of the interval. Can
#only be used for one group
#Input:
# values:   - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Contains only one group
# n_all     - integer, total number of observations
# a         - numeric value between 0 and 1, alpha level, default 0.1
#Output:    Standard deviation estimated with the formula derived from interval

int_sd3 <- function(values, var_level  = var_level, n_all) {
  a <- 1 - var_level
  n <- nrow(values)
  indi <- values[, 1] <= values[, 5]

  Vi <- sum(indi * values[, 2]/(n:1)^2)
  z <- qnorm(1 - a/2)

  li <- max(0, 1/2*(1 - z * sqrt(Vi)))
  #lower quantile

  if(is.infinite(values[1, 5])){
    return(c("sigma" = NA))
  }

  fli_ind <- findInterval(li, values[, 4])
  #index lower quantile

  sig <- 1/z * sqrt(n_all) *
    (values[1, 5] - values[fli_ind+1, 1])

  return(c("sigma" = unname(sig)))
}



#wrapper for int_sd several groups in data
#Input:
# values    - matrix created by dat_gen, sorted by sort_data and with added
#             columns by KME. Should contain several groups
# alpha     - numeric value between 0 and 1, alpha level
# group:    - integer vector containing the group of the observations. Default
#             is the third column of the values, the groups drawn by data_gen
#Output:    Standard deviation estimated with interval formula for each group

int_var_groups <- function(values, var_level = var_level, group = values[, 3],
                           var_method) {
  n_all <- nrow(values)
  n_group <- max(group)
  erg <- numeric(n_group)
  med <- numeric(n_group)

  for(i in 1:n_group) {
    values2 <- values[group == i,]

    #calculate for each group

    temp <- do.call(paste0("int_sd",var_method),
                    args = list(values = values2, var_level = var_level, n_all = n_all))
    erg[i] <- temp
    med[i] <- values2[1, 5]

  }

  names(erg) <- 1:n_group
  names(med) <- 1:n_group
  return(list("Median" = med, "Variance" = erg^2))
  #notice that ^2 changes from sd to var...
}

#int_var_groups(test_KME)



#Teststatistic
#Input:
# n         - integer, total number of observations
# t_mat     - matrix, null hypothesis matrix
# var_out   - Output from one of the variance functions boot_var_groups or
#             int_var_groups, containing median and variance
#
#Output: Teststatistic

test_stat <- function(values, t_mat, var_out) {
  n <- nrow(values)
  med_vec <- var_out$Median
  sig_vec <- var_out$Variance
  if(any(is.infinite(med_vec))) {
    stop("Median does not exist in all subgroups!")
  }
  if(any(is.na(sig_vec))) return(NA)

  return(n * t(t_mat %*% med_vec) %*%
           MASS::ginv(t_mat %*% diag(sig_vec) %*% t(t_mat)) %*%
           (t_mat %*% med_vec))
}


#wrapper to get the teststatistics directly from the sorted data
#Input
# values:   matrix. The data to start with
# group:    integer vector containing the group of the observations. Default is
#           the third column of the values, the groups drawn by data_gen
wrap_sim2 <- function(values, group = values[, 3], t_mat_list, var_method, var_level){

  C_mat <- function(x){
    t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
  }

  t_mat_list <- lapply(t_mat_list, C_mat)

  values <- sort_data(values)
  values[, 3] <- group
  values_KME <- KME(values, group = values[, 3])

  var_int <- int_var_groups(values_KME,var_level = var_level, group = values[ ,3], var_method = var_method)
  erg_int <- list()
  for( i in 1:length(t_mat_list)){
    erg_int[[i]] <- test_stat(values, t_mat_list[[i]], var_out = var_int)
  }


  out <- c( unlist(erg_int))
  names(out) <- c(paste0("int","_" , 1:length(t_mat_list)))

  return(out)
}


#wrap_sim(values = test_dat,
#         t_mat = null_mat_clas(3))


#Permutations
#Input:
# values    - Matrix, Data to be entered in Simulation
# nperm    - Integer, Number of permutations
# t_mat     - matrix, test matrix
#Output
# Vector with the quantiles of the test-statistics and number of permutations
#

perm_fun <- function(values, nperm, t_mat_list , var_method, var_level = var_level) {

  C_mat <- function(x){
    t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
  }

  t_mat_list <- lapply(t_mat_list, C_mat)

  values2 <- sort_data(values)
  group_org <- values2[, 3]
  group_new <- replicate(nperm, sample(group_org))
  test_stat_erg <- apply(group_new, 2,
                         function(x) wrap_sim2(values = values2, group = x,
                                              t_mat_list = t_mat_list, var_method = var_method, var_level = var_level) )
   if(length(t_mat_list)==1){
     test_stat_erg <- t(test_stat_erg)
     rownames(test_stat_erg) <- "int_1"
   }


  return(list(test_stat_erg = test_stat_erg ) )
}

