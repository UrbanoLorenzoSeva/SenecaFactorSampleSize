##############################################################
#
#          S E N E C A    E S T I M A T E     R   C O D E
#
#                 AUTHORS:   URBANO LORENZO-SEVA
#                            PERE J. FERRANDO
#
#                     URV, TARRAGONA (SPAIN)
#                       DATE: 14/04/2023
#
# DATA EXAMPLE: 6 ITEMS AND N = 459
#
# DATA TEXT FILE USED:    example.dat
#
# EXECUTE THE CODE AS:    source("SenecaEstimate.r")
#
##############################################################

##############################################################
# UPDATE HERE THE NAME OF YOUR INPUT DATA FILE
# A TEXT FILE WITH NO LABELS IS EXPECTED
#
filein <- "example.dat"
# 
# UPDATE HERE THE NUMBER OF RANDOM SAMPLES (AT LEAST 100)
# 
K <- 1000
# 
#
# UPDATE HERE THE NUMBER OF FACTOR EXPECTED IN THE FACTOR MODEL
# IF THE VALUE IS ZERO, THEN A FULLY EXPLORATORY APPROACH IS
# ASSUMED
# 
p <- 1
#
# UPDATE HERE THE THRESHOLD VALUE FOR RMSR 
# IT MUST BE A VALUE BETWEEN 0.001 and 0.100
# 
rmsr_t <- 0.03
#
##############################################################

############## SENECA ESTIMATE CODE ##########################
# FUNCTIONS
PrincipalAxes<-function(R,p) {
m <- ncol(R)
IR <-solve(R)
RM <- matrix(0,m,m)
for (i in 1:m){
     RM[i,i] <- 1-1/IR[i,i]
}
RC <- R - diag(diag(R)) + RM
s <- svd(RC)
Kq <- as.matrix(s$v[1:m,1:p])
if (p == 1){
   Dq <- as.matrix(s$d[1:p])
} else{
   Dq <- as.matrix(diag(s$d[1:p]))
}
B <- Kq %*% solve(Dq)^(0.5)
A <- t(solve(t(B)%*%R%*%B) %*% t(B) %*% R )
A <- A * as.matrix(numeric(m)+1)%*%sign(colSums(A))
return(A)
}
KaiserRule<-function(R) {
m <- ncol(R)
eigen <- eigen(R)
ld <-0
  for (i in 1:m){
   if ( eigen$values[i] > 1) {
     ld <- ld +1
   }
  }
return(ld)  
}
rmsr<-function(A,B) {
m <- ncol(A)
h<-0
r<-0
 for (i in 1:(m-1)){
   for (j in (i+1):m){
       r <- r + (A[i,j] - B[i,j])^2
       h <- h+1
  }
 }
 r <- sqrt(r/h)
return(r)
}
izfisher<-function(z) {
r <- (exp(2*z)-1)/(exp(2*z)+1)
return(r)
}
zfisher<-function(r) {
z <- (log(1+r) - log(1-r))/2
return(z)
}
npd<-function(R) {
h<-0
m <- ncol(R)
eigen <- eigen(R)
if ( eigen$values[m] < 0.001) h <- 1 
return(h)
}
######################################
# IN PARAMETERS CHECK 
if (K<100) k <- 100
if (p<0) p <- 0
if (rmsr_t < 0.001) rmsr_t <- 0.001
if (rmsr_t > 0.100) rmsr_t <- 0.100

max_iter <- 20

X <- read.table(filein)
N <- nrow(X)
m <- ncol(X)
if (p>m) p <- m

# SIZE OF THE SMALLEST SAMPLE CONSIDERED
Nmin <- 100 + m*2

# OBSERVED PEARSON CORRELATION MATRIX 
R <- cor(X)
######################################

again <-1
while (again == 1) {
 zse <- 1.96 * 1/sqrt(N-3)
 RP <- R
 for (i in 1:(m-1)){
   for (j in (i+1):m){
     zi <- zfisher(RP[i,j])
     RP[i,j] <- izfisher(min(zi - zse,zi+ zse))
     RP[j,i] <-RP[i,j]
   } 
 }
 if (npd(RP) == 1){
   N <- round(N+N*0.5)
 } else{
   again <- 0
 }
}

if (p==0) {
 KRule <- KaiserRule(RP) 
} else {
 KRule <- p 
}

AP <-PrincipalAxes(RP,KRule)
RRP <- AP %*% t(AP)
L <- chol(RP)
NP <- 100000
Na <- 0
Nf <- matrix(0,K,1)
repi <- 1
NCONV <- 0   
alive=1
alive2=1
while (repi <= K) {

  if (alive2==4) {
     cat("\014")
     print(sprintf(" The size of sample estimated so far is %i (%i%%)",Na,round(repi/K*100)))
     alivestr <- switch(alive," Computing... \U002B                                  "," Computing... \U00D7                                  ")
     print(alivestr)
     alive <- alive+1
     alive2 <- 1
     if (alive>2) alive <- 1
  } else alive2 <- alive2+1

  again <- 1
  iter <- 0
  while (again == 1) {
    Z = matrix(0,NP,m)
    for (i in 1:m){
      Z[,i]=rnorm(NP,0,1)
    }
    XP <- Z %*% L
    #FIT IN THE LARGEST SAMPLE
    RPi <- cor(XP)
    APi <-PrincipalAxes(RPi,KRule)      
    RMSR_max<- rmsr(RRP, APi %*% t(APi))
    if (RMSR_max > rmsr_t) {
      NP <- round(NP+NP*0.5)
    } else {
      again <- 0
    }
    iter <- iter+1
    if (iter > 10) {
      #THE ESTIMATE SIZE OF SAMPLE IS TOO LARGE, NO CONVERGENCE MUST BE REPORTED
       Na <- NP
       NCONV <- 1
       repi <- K   
    }
  }
  if (NCONV == 0) {
    #FIT IN THE SMALLEST SAMPLE
    Xi <- as.matrix(XP[1:Nmin,1:m])
    Ri <- cor(Xi)
    Ai <-PrincipalAxes(Ri,KRule)   
    RMSR_min<- rmsr(RRP, Ai %*% t(Ai))    
    if (RMSR_min < rmsr_t) {
        #THE SMALLEST SAMPLE IS ALREADY GOOD ENOUGHT, THERE IS NOT NEED TO LOOK FOR ANYMORE
        Nm <- Nmin
    } else {
        #WE ARE LOOKING FOR THE FIRST SAMPLE SIZE THAT ACHIEVES THAT RMSR  < RMSR_T
        Nm <- Nmin
     
        inc <- round((NP-Nm)/4)
        again0 <- 1
        iter <- 0
        while (again0 == 1) {
        
            #LET'S FIND IN FIVE POINTS THE FIRST TIME THAT RMSR < RMSR_T
            again <- 1
            while (again == 1) {           
                cat("\014")
                print(sprintf(" The size of sample estimated so far is %i (%i%%)",Na,round(repi/K*100)))
                alivestr <- switch(alive," Computing... \U002B                                  "," Computing... \U00D7                                  ")

                print(alivestr)
                alive <- alive+1
                if (alive>2) alive <- 1
            
                Xi <- as.matrix(XP[1:Nm,1:m])
                if ((Nm == Nmin) | (Nm == NP)) {
                    if (Nm == Nmin) RMSR <- RMSR_min
                    else RMSR <- 0
                } else {
                    Ri <- cor(Xi)
                    Ai <-PrincipalAxes(Ri,KRule)                       
                    RMSR <- rmsr(RRP, Ai %*% t(Ai))
                }
                if (RMSR < rmsr_t) {
                    again <- 0
                } else {
                  Nm <- Nm + inc
                  if (Nm > NP) {
                      Nm <- NP
                  }
                }
            }
            
            iter <- iter + 1
            if (iter > max_iter) {
              again0 <- 0
            } else {
              #LET'S MAKE THE INTERVAL SMALLER AND LET'S LOOK FOR IN THE FIVE PREVIOUS POINTS
              inc <- round(inc/2)
              if (inc < 5) {
                  again0 <- 0
              } else {
                if (round(Nm - inc*4) < Nmin) {
                    inc <- round((Nm-Nmin)/4)
                    if (inc < 2) {
                        inc <- 2
                    }
                    Nm <- Nmin    
                } else {
                    Nm = round(Nm - inc*4)
                }
              }
                   
            } 
        
        }
    }
    #LET'S KEEP DATA SORTENED
    if (repi == 1) {
        Nf[repi,1] <- Nm
        repi <- repi +1
    } else {
      j <- 0
      for (i in 1:repi){
          if (j == 0) {
              if (Nf[i,1] > Nm) {
                  temp <- Nf[i,1]
                  Nf[i,1] <- Nm
                  j <- 1
              } 
          } else {
              temp2 <- Nf[i,1]
              Nf[i,1] <- temp
              temp <- temp2
          }
      }

      if (j==0) {
         Nf[repi,1] <- Nm
      } 
      repi <- repi + 1
      #ESTIMATE SIZE OF SAMPLE AS THE MEDIAN ROUNDED UP TO THE TEN
      Na <- Nf[round(repi/2),1]
      Na <- ceiling(Na/10)*10
    }
         
  }
}

cat("\014")
print("Job done!                                                                ")
print("                                                                         ")


   
####################   END OF LOSEFER CODE  ####################################

print("#########################################################################")
print("#                                                                        ")
print("#        S E N E C A    E S T I M A T E    R   C O D E                   ")
print("#                                                                        ")
print("#               AUTHORS:   URBANO LORENZO-SEVA                           ")
print("#                          PERE J. FERRANDO                              ")
print("#                                                                        ")
print("#                    URV, TARRAGONA (SPAIN)                              ")
print("#                                                                        ")
print("#                       DATE: 14/04/2023                                 ")
print("#                                                                        ")
print("#########################################################################")
print("#                                                                        ")
print(paste("# Filein                : ",filein," "))
print(paste("# Cases                 : ",N," "))
print(paste("# Variables             : ",m," "))
print(paste("# Random samples        : ",K," "))
if (p>0) {
  print("# Aim of the analysis   :  Confirmatory")
  print(paste("# Number of factors     : ",p," "))
} else {
  print("# Aim of the analysis   :  Fully exploratory")
}
print(sprintf("# RMSR threshold        :  %5.3f",rmsr_t))
print("#                                                                        ")
print("#                                                                        ")
print(sprintf("#Seneca estimate sample size =   %5i                                   ",Na))
print("#                                                                        ")
print("#                                                                        ")
print("#########################################################################")
