############################################################ MACROECONOMETRIA - TP 2 ############################################################ 

rm(list=ls())
gc()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("LP.R")
source("TP2_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap2.R")
source("PS2_SVAR_Plots.R")

############################################################################################################################################
########################################################## PUNTO 1 #########################################################################
############################################################################################################################################

# 1 A----
# Data 

Yl.f <- cbind(gc, yk) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 04))
Yd <- window(Yd.f, start = c(2004, 02), end = c(2019, 04))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
VAR <- VAR(Ydt, p = p, type = "const")
summary(VAR)

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

# SVAR t0 impact matrix (Cholesky decomposition)
S <- t(resid(VAR)) %*% resid(VAR) / N
P.chol <- t(chol(S)) # Cholesky decomposition
S

# SVAR t0 impact matrix (implied by AB model)
P <- solve(SVAR$A, SVAR$B) # inv(A) %% B
S.SVAR <- P %*% t(P)
S.SVAR

# Other SVAR parameters
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) # SVAR
pars.R
pars.S

# SVAR analysis

H <- 12 # Horizon
H_ERPT <- 50 # Horizon for Multiplicador fiscal... OJO QUE TODOS LOS NOMBRES APARECEN COMO ERPT, PERO ES PARA NO MODIFICAR TODO EL CODIGO...SOLO MODIFICAMOS LA ETIQUETA DEL PLOT!!

# Bootstrap inference

a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications

# Bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
  } else {
    p <- popt$selection[3] # SC
  }
  Y0 <- Y[1:pmax, ]
  Yt <- Y[(pmax - p + 1):nrow(Y), ]
  VAR <- VAR(Yt, p = p, type = "const")
  #VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
  print(roots(VAR, modulus = TRUE))
  print(serial.test(VAR, lags.bg = r, type = "ES"))
  SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
  Yb <- boot.rb.replicate(VAR, Y0, pmax, R)
  if (cumulative == FALSE) {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = FALSE)
  } else {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = TRUE)
  }
  return(SVAR.ERPT.boot)
}

# ERPT in log-levels
SVAR.ERPT.boot.lvl <- assess.erpt(Amat, Bmat, Yl, pmax, 12, H_ERPT, R, a)
plot.erpt.boot(SVAR.ERPT.boot.lvl, H_ERPT)




# 1 B----
# Data 

Yl.f <- cbind(gs, yk) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 04))
Yd <- window(Yd.f, start = c(2004, 02), end = c(2019, 04))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
VAR <- VAR(Ydt, p = p, type = "const")
summary(VAR)

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

# SVAR t0 impact matrix (Cholesky decomposition)
S <- t(resid(VAR)) %*% resid(VAR) / N
P.chol <- t(chol(S)) # Cholesky decomposition
S

# SVAR t0 impact matrix (implied by AB model)
P <- solve(SVAR$A, SVAR$B) # inv(A) %% B
S.SVAR <- P %*% t(P)
S.SVAR

# Other SVAR parameters
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) # SVAR
pars.R
pars.S

# SVAR analysis

H <- 12 # Horizon
H_ERPT <- 50 # Horizon for Multiplicador fiscal... OJO QUE TODOS LOS NOMBRES APARECEN COMO ERPT, PERO ES PARA NO MODIFICAR TODO EL CODIGO...SOLO MODIFICAMOS LA ETIQUETA DEL PLOT!!

# Bootstrap inference

a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications

# Bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
  } else {
    p <- popt$selection[3] # SC
  }
  Y0 <- Y[1:pmax, ]
  Yt <- Y[(pmax - p + 1):nrow(Y), ]
  VAR <- VAR(Yt, p = p, type = "const")
  #VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
  print(roots(VAR, modulus = TRUE))
  print(serial.test(VAR, lags.bg = r, type = "ES"))
  SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
  Yb <- boot.rb.replicate(VAR, Y0, pmax, R)
  if (cumulative == FALSE) {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = FALSE)
  } else {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = TRUE)
  }
  return(SVAR.ERPT.boot)
}

# ERPT in log-levels
SVAR.ERPT.boot.lvl <- assess.erpt(Amat, Bmat, Yl, pmax, 12, H_ERPT, R, a)
plot.erpt.boot(SVAR.ERPT.boot.lvl, H_ERPT)




# 1 C----
# Data 

Yl.f <- cbind(gk, yk) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2007, 01), end = c(2019, 04))
Yd <- window(Yd.f, start = c(2007, 02), end = c(2019, 04))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
VAR <- VAR(Ydt, p = p, type = "const")
summary(VAR)

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

# SVAR t0 impact matrix (Cholesky decomposition)
S <- t(resid(VAR)) %*% resid(VAR) / N
P.chol <- t(chol(S)) # Cholesky decomposition
S

# SVAR t0 impact matrix (implied by AB model)
P <- solve(SVAR$A, SVAR$B) # inv(A) %% B
S.SVAR <- P %*% t(P)
S.SVAR

# Other SVAR parameters
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) # SVAR
pars.R
pars.S

# SVAR analysis

H <- 12 # Horizon
H_ERPT <- 50 # Horizon for Multiplicador fiscal... OJO QUE TODOS LOS NOMBRES APARECEN COMO ERPT, PERO ES PARA NO MODIFICAR TODO EL CODIGO...SOLO MODIFICAMOS LA ETIQUETA DEL PLOT!!

# Bootstrap inference

a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications

# Bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
  } else {
    p <- popt$selection[3] # SC
  }
  Y0 <- Y[1:pmax, ]
  Yt <- Y[(pmax - p + 1):nrow(Y), ]
  VAR <- VAR(Yt, p = p, type = "const")
  #VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
  print(roots(VAR, modulus = TRUE))
  print(serial.test(VAR, lags.bg = r, type = "ES"))
  SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
  Yb <- boot.rb.replicate(VAR, Y0, pmax, R)
  if (cumulative == FALSE) {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = FALSE)
  } else {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = TRUE)
  }
  return(SVAR.ERPT.boot)
}

# ERPT in log-levels
SVAR.ERPT.boot.lvl <- assess.erpt(Amat, Bmat, Yl, pmax, 12, H_ERPT, R, a)
plot.erpt.boot(SVAR.ERPT.boot.lvl, H_ERPT)



