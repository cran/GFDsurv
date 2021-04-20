tr <- function(A){
  sum(diag(A))
}



p.func <- function(times1, times2, KME1, KME2){
  # Berechnet die p's fuer die einzelnen Gruppenpaare
  times <- sort( unique( c(times1, times2)))  #Vektor, der die entscheidenden Zeitpunkte der 1 UND 2 SP enth?lt
  # Weise den Zeitpunkten aus time das Intervall aus timej zu, in dem sie stecken:
  # (-infty,timej_1) -> 0, [timej_1,timej_2) -> 1 usw.
  find1 <- findInterval( times, times1)
  find2 <- findInterval( times, times2)
  zeros1 <- sum(find1 == 0)
  zeros2 <- sum(find2 == 0)

  # Vektor der Ueberlebensw. (KME) zu ALLEN aufeinander folgenden Zeitpunkten
  ntot <- length(times)

  if(zeros1==ntot){
    surv.vec1 <- rep(1,zeros1)
  }else{
    surv.vec1 <- c( rep(1,zeros1) , KME1[find1[(zeros1+1):ntot]] )
  }
  if(zeros2==ntot){
    surv.vec2 <- rep(1,zeros2)
  }else{
    surv.vec2 <- c( rep(1,zeros2) , KME2[find2[(zeros2+1):ntot]] )
  }

  # Bilde nun die normalisierte Version des Integranden
  surv.vec1.norm <- 0.5 * ( surv.vec1 + c(1, surv.vec1[-ntot]) )

  # Zuwaechse der KMEs
  surv.diff2 <- diff(c(1,surv.vec2))
  - surv.vec1.norm %*% surv.diff2
}


cov.matrix <- function(n, times, KME, at.risk, oneMinusDeltaA, G.list=NULL, finds=NULL, zeros=NULL){
  ## times: falls G.list == NULL, dann ist das eine Liste mit den jeweiligen Ereigniszeiten
  ##        sonst ist es ein Vektor der sortierten, einmaligen Ereigniszeiten von allen Gruppen.
  ## finds und zeros sind fuer die Teststat NULL, sonst (beim Wild BS) verschieden von NULL
  no.groups <- length(n)

  if(is.null(G.list)){
    all.times <- times[[1]]
    for(j in 2:no.groups){
      all.times <- c(all.times, times[[j]])
    }
    all.times <- sort(unique(all.times))
    n.times <- length(all.times)

    V <- matrix(0, no.groups, no.groups)
    Gamma.norm <- list()
    length(Gamma.norm) <- no.groups
    # Gesamtzahl Individuen in allen Gruppen
    n.tot <- sum(n)

    # Folgende arithmetischen Mittel werden spaeter benoetigt
    weighted.Gamma.mean <- matrix(0, n.times, n.times)
    Surv.mean <- rep(0, n.times)

    # Speichere "find1" usw in einer Liste ab.
    # Gamma ist eine Liste von gruppenspezifischen Cov-Matrizen
    finds <- zeros <- Gamma <- diffs.Surv <- list()
    length(finds) <- length(zeros) <- length(Gamma) <- length(diffs.Surv) <- no.groups

    for(j in 1:no.groups){
      finds[[j]] <- findInterval(all.times, times[[j]])
      zeros[[j]] <- sum(finds[[j]] == 0)
      if(zeros[[j]] == n.times){
        KME.temp <- rep(1, zeros[[j]])
        int <- rep(0, n.times)
      }else{
        KME.temp <- c(rep(1, zeros[[j]]), KME[[j]][finds[[j]]])
        diffs.Surv[[j]] <- diff(c(1, KME.temp))
        int <- n[[j]] * cumsum( (1-oneMinusDeltaA[[j]]) / (at.risk[[j]] * oneMinusDeltaA[[j]]))

        if(length(oneMinusDeltaA[[j]]) == 1){
          # In dem Fall wuerde die Max-Findung im anderen Fall leer ausgehen
          # und generell wuerde das Integral direkt auf unendlich springen
          int <- rep(0, n.times)
        }else{
          int[which(oneMinusDeltaA[[j]] == 0)] <- max(int[which(oneMinusDeltaA[[j]] > 0)])
          int <- c(rep(0, zeros[[j]]), int[finds[[j]]])
        }
      }

      int.mat <- pmin( replicate(n.times, int), t(replicate(n.times, int)))
      Gamma[[j]] <- (KME.temp %*% t(rep(1, n.times))) * int.mat * t(KME.temp %*% t(rep(1, n.times)))
      Gamma.norm[[j]] <- 0.25 * ( Gamma[[j]] + cbind(rep(0, n.times), Gamma[[j]][,-n.times])
                                  + rbind(rep(0, n.times), Gamma[[j]][-n.times,])
                                  + rbind(rep(0, n.times), cbind(rep(0, n.times-1), Gamma[[j]][1:(n.times-1),1:(n.times-1)])) )
      weighted.Gamma.mean <- weighted.Gamma.mean + n.tot / n[[j]] * Gamma.norm[[j]]
      Surv.mean <- Surv.mean + KME.temp
    }
    Surv.mean <- Surv.mean / no.groups
    diffs.Surv.mean <- diff(c(1,Surv.mean))
  }else{
    n.times <- length(times)
    all.times <- times

    V <- matrix(0, no.groups, no.groups)
    Gamma.norm <- list()
    length(Gamma.norm) <- no.groups
    # Gesamtzahl Individuen in allen Gruppen
    n.tot <- sum(n)

    # Folgende arithmetischen Mittel werden spaeter benoetigt
    weighted.Gamma.mean <- matrix(0, n.times, n.times)
    Surv.mean <- rep(0, n.times)

    # Speichere "find1" usw in einer Liste ab.
    # Gamma ist eine Liste von gruppenspezifischen Cov-Matrizen
    Gamma <- diffs.Surv <- list()
    length(Gamma) <- length(diffs.Surv) <- no.groups
    for(j in 1:no.groups){
      if(zeros[[j]] == n.times){
        KME.temp <- rep(1, zeros[[j]])
        int <- rep(0, n.times)
      }else{
        KME.temp <- KME[[j]]
        ## IN WBS.TESTSTAT WURDE FINDS SCHON ANGEWANDT!
        ## ABER NUR AUF DEN KME!

        diffs.Surv[[j]] <- diff(c(1, KME.temp))
        at.risk.temp <- at.risk[[j]][finds[[j]]]
        oMDA.temp <- oneMinusDeltaA[[j]][finds[[j]]]

        int <- n[[j]] * cumsum( G.list[[2*j]] / (at.risk[[j]]^2 * oneMinusDeltaA[[j]]) )

        if(length(oneMinusDeltaA[[j]]) == 1){
          # In dem Fall wuerde die Max-Findung im anderen Fall leer ausgehen
          # und generell wuerde das Integral direkt auf unendlich springen
          int <- rep(0, n.times)
        }else{
          int[which(oneMinusDeltaA[[j]] == 0)] <- max(int[which(oneMinusDeltaA[[j]] > 0)])
          int <- c(rep(0, zeros[[j]]), int[finds[[j]]])
        }
      }

      int.mat <- pmin( replicate(n.times, int), t(replicate(n.times, int)))
      Gamma[[j]] <- (KME.temp %*% t(rep(1, n.times))) * int.mat * t(KME.temp %*% t(rep(1, n.times)))
      Gamma.norm[[j]] <- 0.25 * ( Gamma[[j]] + cbind(rep(0, n.times), Gamma[[j]][,-n.times])
                                  + rbind(rep(0, n.times), Gamma[[j]][-n.times,])
                                  + rbind(rep(0, n.times), cbind(rep(0, n.times-1), Gamma[[j]][1:(n.times-1),1:(n.times-1)])) )
      weighted.Gamma.mean <- weighted.Gamma.mean + n.tot / n[[j]] * Gamma.norm[[j]]
      Surv.mean <- Surv.mean + KME.temp
    }
    Surv.mean <- Surv.mean / no.groups
    diffs.Surv.mean <- diff(c(1,Surv.mean))
  }

  # Berechne nun die einzelnen Eintrage von V
  for(j in 1:no.groups){
    for(k in j:no.groups){
      if(is.null(diffs.Surv[[k]])){
        # wenn der KME NULL ist... brechen wir ab und geben die Identitaet zurueck
        V <- diag(rep(1,no.groups))
        break
      }
      if(j == k){
        int1 <- n.tot / n[[j]] * ( t(diffs.Surv.mean - 2 / no.groups * diffs.Surv[[j]])
                                   %*% Gamma.norm[[j]] %*% diffs.Surv.mean )
        int2 <- 1 / no.groups^2 * ( t(diffs.Surv[[j]]) %*% weighted.Gamma.mean %*% diffs.Surv[[j]] )
        # Negative Varianzen sind nicht erlaubt...
        V[j,j] <- max(int1 + int2, 0)
      }else{
        int1 <- - n.tot / n[[k]] / no.groups * ( t(diffs.Surv[[j]]) %*% Gamma.norm[[k]] %*% diffs.Surv.mean )
        int2 <- - n.tot / n[[j]] / no.groups * ( t(diffs.Surv[[k]]) %*% Gamma.norm[[j]] %*% diffs.Surv.mean )
        int3 <- 1 / no.groups^2 * ( t(diffs.Surv[[j]]) %*% weighted.Gamma.mean %*% diffs.Surv[[k]] )
        V[j,k] <- V[k,j] <- int1 + int2 + int3
      }
    }
  }
  list(V=V, all.times=all.times, finds=finds, zeros=zeros)
}


####3.25.Schritt: Erzeuge DDMB-Gewichte
DDMB_Weights <- function(n, ntime, Gewichte, at.risk, n.ev){
  # Erzeugt die gewuenschten DDMB-Gewichte.
  # Und zwar als Summe: Die Gewichte derselben Zeitpunkte werden summiert.
  # Dieselbe Summe der quadrierten Gewichte wird ebenso berechnet.
  out <- list()
  length(out) <- 2 * length(n)

  for(j in 1:length(n)){
    if( Gewichte == "pois" ){
      G <- sapply(n.ev[[j]], function(k) c(sum(z <- rpois(k,1) - 1), sum(z^2)))
    }else{ if( Gewichte == "norm" ){
      G <- sapply(n.ev[[j]], function(k) c(sum(z <- rnorm(k)), sum(z^2)))
    }else{ if( Gewichte == "rade" ){
      G <- sapply(n.ev[[j]], function(k) c(sum(z <- 2*rbinom(k, 1, prob=0.5) - 1), sum(z^2)))
    }else{
      at.risk.tmp <- sapply(at.risk[[j]], max, 1)
      G <- sapply(n.ev[[j]], function(k) c(sum(z <- rbinom(n = k, size = at.risk[[j]], prob = 1/at.risk[[j]]) - 1), sum(z^2)))
    }}}

    # Hier hatte ich zunaechst einen Indizierungsfehler begangen:
    # Durch das sapply entsteht eigentlich eine Matrix...
    # Man kann es aber auch als Vektor indizieren,
    # dabei wechseln sich allerdings die quadrierten und nicht-quadrierten Werte ab!
    out[[2*j-1]] <- G[2*(1:ntime[[j]])-1]
    out[[2*j]] <- G[2*(1:ntime[[j]])]
  }
  out
}


wild.bs.teststat <- function(at.risk, n, ntime, times, Gewichte, deltaN, oneMinusDeltaA, KME, finds, zeros, T.matrix){
  # Erzeugt die WBS-Version der Teststatistik.
  # Zuerst werden pro Gruppe n zufaellige G erzeugt nach der Verteilung, die durch Gewichte angegeben wird.
  # times enthaelt alle verschiedenen Ereigniszeiten.
  # finds enthaelt alle Ergebnnisse von findInterval() fuer alle Gruppen
  # analog fuer die zeros
  no.groups <- length(n)
  n.times <- length(times)
  wbs.kme <- list()
  length(wbs.kme) <- no.groups
  dphi <- matrix(0, no.groups, no.groups)
  G.list <- DDMB_Weights(n, ntime, Gewichte, at.risk, deltaN)

  for(j in 1:no.groups){
    wbs.kme[[j]] <- KME[[j]] * sqrt(sum(n)) * cumsum(G.list[[2*j-1]] / sqrt(oneMinusDeltaA[[j]] * at.risk[[j]]^2))
    wbs.kme[[j]][which(oneMinusDeltaA[[j]] == 0)] <- 0
    times.temp <- length(wbs.kme[[j]])

    if(zeros[[j]]==length(times)){
      wbs.kme[[j]] <- rep(0, zeros[[j]])
      KME[[j]] <- rep(1, zeros[[j]])
    }else{
      wbs.kme[[j]] <- c(rep(0, zeros[[j]]), wbs.kme[[j]][finds[[j]]])
      KME[[j]] <- c(rep(1, zeros[[j]]), KME[[j]][finds[[j]]])
    }

    # Bilde nun die Normalisierung.
    wbs.kme[[j]] <- (0.5 * (c(0,wbs.kme[[j]]) + c(wbs.kme[[j]],1))[-n.times])
  }

  for(j in 1:(no.groups-1)){
    for(k in (j+1):no.groups){
      dphi[j,k] <- wbs.kme[[j]] %*% diff(c(1,KME[[k]])) - wbs.kme[[k]] %*% diff(c(1,KME[[j]]))
      dphi[k,j] <- - dphi[j,k]
    }
  }

  # Vektorisiere nun dphi... laufe die Zeilen ab. Daher: transponieren!

  temp <- cov.matrix(n, times, KME, at.risk, oneMinusDeltaA, G.list, finds, zeros)
  V <- temp$V

  Edphi <- 1 / no.groups * (diag(rep(1, no.groups)) %x% t(rep(1, no.groups))) %*% c(t(dphi))
  dphiETEdphi <- t(Edphi) %*% T.matrix %*% Edphi
  dphiETEdphi / tr(T.matrix %*% V)
}


DDMB.quant.discrete <- function(at.risk, n, ntime, times, BSiter, Gewichte, deltaN, oneMinusDeltaA, KME, finds, zeros, T.matrix){
  #Input:
  # 1. Eine Liste U erzeugt durch UMatrix(siehe 2.Schritt) bestehend aus der Matrix U und dem Vektor time
  # 2. n Stichprobenumfang
  # 3. Gewichte, um anzugegeben, wie die Gewichte G erzeugt werden. Zurzeit gibt es diese Moeglichkeiten
  # (i) 	pois :  Die G's sind i.i.d. Poissonverteilt zum Parameter 1 verschoben um -1
  #                 (kurz: G ~ Poisson(1) - 1)
  # (ii) 	norm: Die G's sind i.i.d. standardnormalverteilt
  #	(iii)	weird: Die G's sind binomialverteilt mit zufaelligen Param.; siehe ABGK "weird bootstrap"
  #	(iv)	corrLibPois: wie (i) mit korrigierendem Faktor fuer zu liberale Tests
  #	(v)	corrLibNorm: wie (ii) mit korrigierendem Faktor fuer zu liberale Tests
  # (vi)	corrLibWeird: wie (iii) mit korrigierendem Faktor fuer zu liberale Tests
  replicate( BSiter, wild.bs.teststat(at.risk, n, ntime, times, Gewichte, deltaN, oneMinusDeltaA, KME=KME, finds, zeros, T.matrix))
}

#copSANOVA: concordance probability SANOVA
####5.Schritt: Konfidenzintervalle
test.data <- function(data, n, BSiter, Gewichte, c.matrix){
  #Input:
  # 0. data ist eine Liste von data.frames mit ueblichem Format:
  #     $exit fuer die Ereignis-/Zensierzeit
  #     $to fuer den Zensierindikator ("cens") oder den Ereignisindikator ("1")
  # 1. n steht fuer den Stichprobenumfang: davon sind links-trunkierte ggf. abzuziehen!
  # 2. BSiter steht fuer die Anzahl der Monte-Carlo-Schritt um den krit. Wert zu bestimmen
  # 3. alpha entspricht dem Signifikanzniveau
  # 4. Gewichte wird beim BS benoetigt, um die Verteilung der Gewichte G anzugeben (m?glich ist Gewichte = "pois", "rade", "norm", siehe oben es geht auch weird ...
  # 5. c.matrix ist die Kontrastmatrix zum Testen von Cp=0

  C_mat <- function(x){
    t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
  }

  T.matrix <- lapply(c.matrix, C_mat)


  no.groups <- length(n)
  p.matrix <- matrix(0, no.groups, no.groups)
  diag(p.matrix) <- 0.5

  # Speichere alles in Listen ab.
  dummy <- list()
  length(dummy) <- no.groups
  times <- at.risk <- n.event <- oneMinusDeltaA <- KME <- KMEleft <- list()
  # Gesamtvektor aller Ereigniszeitpunkte
  times.tot <- (data[[1]]$exit)[which((data[[1]]$to) != "cens")]

    ntime <- numeric(no.groups)

  anyGroupEmpty <- FALSE

  for(j in 1:no.groups){
    # von den unzensierten Ereigniszeiten die einzigartigen hernehmen
    times[[j]] <- sort(unique((data[[j]]$exit)[which((data[[j]]$to) != "cens")]))
    times.tot <- c(times.tot,(data[[1]]$exit)[which((data[[1]]$to) != "cens")])
    ntime[j] <- length(times[[j]])
    at.risk[[j]] <- sapply(1:(ntime[j]), function(i){sum((data[[j]])$exit >= (times[[j]])[i])})
    n.event[[j]] <- sapply(1:(ntime[j]), function(i){sum(((data[[j]])$exit == (times[[j]])[i]) & ((data[[j]])$to != "cens"))})
    oneMinusDeltaA[[j]] <- 1 - n.event[[j]] / at.risk[[j]]
    KME[[j]] <- cumprod(oneMinusDeltaA[[j]])
    if(length(times[[j]])==0){
      anyGroupEmpty <- TRUE
    }
  }

  times.tot <- sort(unique(times.tot))


  for(j in 1:(no.groups-1)){
    for(k in (j+1):no.groups){
      p.matrix[j,k] <- p.func(times1 = times[[j]], times2 = times[[k]], KME1 = KME[[j]], KME2 = KME[[k]])
      p.matrix[k,j] <- 1 - p.matrix[j,k]
    }
  }


  # Versuche mit der naechsten Abfrage zu verhindern, dass es irgendwo sonst zu Fehlern kommt.
  if(anyGroupEmpty){
    0 #list(test = 0, FNT = 0, DDMB.teststats = 0)
  }else{

    # Berechne die Vergleiche "Einzelne vs. Gruppen-Mittel"
    # p = E w   im Papier
    p.vec <- 1 / no.groups * (diag(rep(1, no.groups)) %x% t(rep(1, no.groups))) %*%  c(t(p.matrix))


    temp <- cov.matrix(n, times, KME, at.risk, oneMinusDeltaA, NULL, NULL, NULL)
    V <- temp$V

    FNT_mat <- function(x){
      c(sum(n) * t(p.vec) %*% x %*% p.vec / tr(x %*% V))
    }

    FNT <-  lapply(T.matrix, FNT_mat)



    ## Bestimme das Quantil mittels DDMB


    DDMB.quant.discrete_mat <- function(x){
      DDMB.quant.discrete(at.risk, n, ntime, times = temp$all.times, BSiter, Gewichte, deltaN = n.event, oneMinusDeltaA = oneMinusDeltaA, KME=KME, finds = temp$finds, zeros = temp$zeros, x)
      }
    DDMB.teststats <- lapply(T.matrix, DDMB.quant.discrete_mat)
    #q.alpha.discrete <- quantile(DDMB.teststats, probs = 1-alpha)

    erg <- c()

    for(i in 1:length(DDMB.teststats)){
      erg[i] <- sum(c(DDMB.teststats[[i]], FNT[[i]]) >= FNT[[i]]) / (BSiter+1)
    }

    return(list("value" = erg,"test_statistics" = unlist(FNT)))

  }
}


