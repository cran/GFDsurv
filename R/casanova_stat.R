# function for calculation of the test statistic
y_fun <- function(k, group_mat,n){
  group_k <- group_mat[k,]
  n_k <- sum(group_k)
  #Y_j(X_i) :=
  c(n_k, n_k - cumsum(group_k[-n]))
}

int_fun <- function(k, event,  n.risk_vec, m, weight_Comb,len_sig, Comb_sig, n.risk_prod,weight_KME, group_mat){
  group_k <- group_mat[k,]
  n.risk <- n.risk_vec[,k]

  delta <- ifelse( n.risk == 0, 0, group_k*event/n.risk)

  delta_sig <- ifelse( n.risk == 0, 0, delta/n.risk)

  #Kombinieren der verschd. Gewichtsfunktionen W_tilde und W_n
  weight_combined <- apply(weight_Comb,1, function(x) weight_KME[[x[1]]]*n.risk_prod[[x[2]]])

  stat <- apply(weight_combined,2, function(x) sum(x*delta))


  ################################ Varianz und Kovarianz der Gewichtsfunktionen #####################################

  sig <- sapply(1:len_sig, function(x) sum(weight_combined[,Comb_sig[x,1]]*weight_combined[,Comb_sig[x,2]]*delta_sig))
  c(stat,sig)
}



#### extra stat funktion #####
stat <- function(C,group_mat, n, m, event, ngroup, Comb_sig,len_sig, weight_KME, weight_Comb){


  #n = n_all
  n.risk_vec <- sapply(1:ngroup, y_fun, group_mat = group_mat, n = n)

  #Gewichtsfkten W_n
  n.risk_prod <- list(apply( n.risk_vec, 1, function(x){prod(x)/max(1,sum(x)^{ngroup - 1})}))

  #Alternativen
  #    n.risk_prod2 <- apply( n.risk_vec, 1, function(x){prod(x)/max(sum(x),1)})
  #    n.risk_prod3 <- apply( n.risk_vec, 1, function(x){(max(x)*min(x))/max(1,sum(x))})

  # In Zeile 1-m: Z_wn ;  in Zeile (m+1)-2m:Sigma_Wn
  out <- sapply(1:ngroup, int_fun, event = event, n.risk_vec = n.risk_vec, m=m, weight_Comb=weight_Comb,len_sig = len_sig,
                Comb_sig = Comb_sig, n.risk_prod = n.risk_prod, weight_KME  = weight_KME , group_mat = group_mat)

  # print(out[1:m,])
  # print(C)
  ### C ist in Paper Matrix T
  # Erhalte T_z f?r jedes einzelne Z
  T_Z <- apply(out[1:m,],1,function(x) C%*% matrix(x , nrow = ngroup))

  # F?ge T_z zu T(m)_Z zusammen
  T_Z_m <- matrix(c(T_Z),m*ngroup)

  # Erhalte T_Sigma_t
  out_sig <- split(t(out[(m+1):(m+len_sig),]),rep(1:ncol(t(out[(m+1):(m+len_sig),])), each = nrow(t(out[(m+1):(m+len_sig),]))))
  T_Sigma_T <- lapply(out_sig, function(x) C%*% diag(x)%*% t(C))

  ### Diagonalmatrix f?r Standardvarianzen
  diag_T_Sigma_T <- Reduce(magic::adiag,T_Sigma_T[1:m])
  ### Kovarianzen in Matrix einf?gen
  for(i in (m+1):len_sig){
    diag_T_Sigma_T[(ngroup*(Comb_sig[i,1]-1)+1):(ngroup*Comb_sig[i,1]),(ngroup*(Comb_sig[i,2]-1)+1):(ngroup*Comb_sig[i,2])] <- T_Sigma_T[[i]]
    diag_T_Sigma_T[(ngroup*(Comb_sig[i,2]-1)+1):(ngroup*Comb_sig[i,2]),(ngroup*(Comb_sig[i,1]-1)+1):(ngroup*Comb_sig[i,1])] <- T_Sigma_T[[i]]
  }

 out_c <-  t(T_Z_m)%*%MASS::ginv(diag_T_Sigma_T)%*% T_Z_m

 out <-  sapply(1:m,function(x) T_Z[,x]%*%MASS::ginv(T_Sigma_T[[x]])%*% T_Z[,x])


 out <- c(out_c,out)
}


stat_1 <- function(C,group_mat, n, m, event, ngroup, Comb_sig,len_sig, weight_KME, weight_Comb){

  n.risk_vec <- sapply(1:ngroup, y_fun, group_mat = group_mat, n = n)

  #Gewichtsfkten W_n
  n.risk_prod <- list(apply( n.risk_vec, 1, function(x){prod(x)/max(1,sum(x)^{ngroup - 1})}))

  # In Zeile 1-m: Z_wn ;  in Zeile (m+1)-2m:Sigma_Wn
  out <- sapply(1:ngroup, int_fun, event = event, n.risk_vec = n.risk_vec, m=m, weight_Comb=weight_Comb,len_sig = len_sig,
                Comb_sig = Comb_sig, n.risk_prod = n.risk_prod, weight_KME  = weight_KME , group_mat = group_mat)




  ### C ist in Paper Matrix T
  # Erhalte T_z f?r jedes einzelne Z
  T_Z <- apply(t(out[1:m,]),1,function(x) C%*% matrix(x , nrow = ngroup))

  # F?ge T_z zu T(m)_Z zusammen
  T_Z_m <- matrix(c(T_Z),m*ngroup)



  # Erhalte T_Sigma_t
  out_sig <- split(t(out[(m+1):(m+len_sig),]),rep(1:len_sig, each = ngroup))
  T_Sigma_T <- lapply(out_sig, function(x) C%*% diag(x)%*% t(C))

  ### Diagonalmatrix f?r Standardvarianzen
  diag_T_Sigma_T <- Reduce(magic::adiag,T_Sigma_T[1:m])
  ### Kovarianzen in Matrix einf?gen

  out_c <-  t(T_Z_m)%*%MASS::ginv(diag_T_Sigma_T)%*% T_Z_m

  out <- out_c


  out <- c(out_c,out)
}


stat_factorial <- function(hypo_matrices,group, event,n,  n_all, w, nperm){
KME <- (1 - c(1, cumprod(1- event/((n_all+1)-1:n_all))))[1:n_all]

wFn <- list()

for(i in 1:length(w)){
  wFn[[i]] <- mapply(w[[i]],KME)
}

# Kaplan-Meier + Gewichtungen, die unabh?ngig von Gruppe
weight_KME <- wFn

#Anzahl an Gewichten
m <- length(wFn)
#Anzahl an Sigma Kobinationen
Comb_sig <- t(matrix(rep(1:m,each = 2),2))
if(m > 1){
  Comb_sig <- rbind(Comb_sig ,t(combn(1:m,2)))
}
#L?nge der Sigmakombinationen
len_sig <- length(Comb_sig[,1])

#Kombinationen aller Gewichtsfkten
weight_Comb <- expand.grid(1:m,1)
ngroup <- length(n)


matrix1 <- matrix(rep(1:ngroup, each = n_all), nrow = ngroup, byrow = TRUE)
group_mat <- matrix1 == matrix(rep(group, ngroup), nrow = ngroup, byrow = TRUE)

C_mat <- function(x){
  t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
}

hypo_matrices <- lapply(hypo_matrices, C_mat)

      ##Statistik berechnen. F?r ein Gewicht andere Vorgehensweise
      if(m == 1){
        stat_Erg <-lapply(hypo_matrices, function(x) stat_1(x,group_mat, n= n_all, m, event, ngroup, Comb_sig, len_sig,
                                                     weight_KME, weight_Comb))
      } else {
        stat_Erg <-lapply(hypo_matrices, function(x) stat(x,group_mat, n=n_all, m, event, ngroup, Comb_sig, len_sig,
                                                   weight_KME, weight_Comb))
      }

      group_perm <- function(group){
        perm <- sample(1:n_all)
        group <- group[perm]
        group_mat <- matrix1 == matrix(rep(group, ngroup), nrow = ngroup, byrow = TRUE)
      }

      #Permutations
      if(m == 1){
        per <- lapply(hypo_matrices, function(x) replicate(nperm, stat_1(x,group_mat = group_perm(group), n=n_all, m, event, ngroup, Comb_sig, len_sig,
                                                                  weight_KME, weight_Comb)))
      } else {
        per <- lapply(hypo_matrices, function(x) replicate(nperm, stat(x,group_mat = group_perm(group), n=n_all, m, event, ngroup, Comb_sig, len_sig,
                                                                weight_KME, weight_Comb)))
      }


results <- list("Stat" = stat_Erg,"Perm"=per)

}
