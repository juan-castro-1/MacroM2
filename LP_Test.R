remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)

Data_FM <- read_csv("Data_FM.csv")

gc <- ts(Data_FM$CP_K / 4, start = c(2004, 01), frequency = 4)
yk <- ts(Data_FM$PIB_K / 4, start = c(2004, 01), frequency = 4)

gc <- window(gc, end = c(2019, 04))
yk <- window(yk, end = c(2019, 04))

source("LP.R")

rto <- gc / yk
rto <- lag1(rto, 1, 1, length(rto))

gc <- diff(log(gc))
yk <- diff(log(yk))

# Example
LP.irf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95)
LP.cirf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, cumulative = TRUE)
LP.mult.adj <- irf.lp.iv.multiplier(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, adj = rto)

# IRFs
LP.irf$pe
LP.cirf$pe
LP.mult.adj$pe

# Plot
plot(0:4, LP.mult.adj$pe[, 1], type = "l", col = "red", xlab = "", ylab = "", ylim = c(min(LP.mult.adj$lb[, 1]), max(LP.mult.adj$ub[, 1])))
lines(0:4, LP.mult.adj$ub[, 1])
lines(0:4, LP.mult.adj$lb[, 1])


################
### PARA GS ####
################

gc <- ts(Data_FM$GS_K / 4, start = c(2004, 01), frequency = 4)
yk <- ts(Data_FM$PIB_K / 4, start = c(2004, 01), frequency = 4)

gc <- window(gc, end = c(2019, 04))
yk <- window(yk, end = c(2019, 04))

source("LP.R")

rto <- gc / yk
rto <- lag1(rto, 1, 1, length(rto))

gc <- diff(log(gc))
yk <- diff(log(yk))

# Example
LP.irf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95)
LP.cirf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, cumulative = TRUE)
LP.mult.adj <- irf.lp.iv.multiplier(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, adj = rto)

# IRFs
LP.irf$pe
LP.cirf$pe
LP.mult.adj$pe

# Plot
plot(0:4, LP.mult.adj$pe[, 1], type = "l", col = "red", xlab = "", ylab = "", ylim = c(min(LP.mult.adj$lb[, 1]), max(LP.mult.adj$ub[, 1])))
lines(0:4, LP.mult.adj$ub[, 1])
lines(0:4, LP.mult.adj$lb[, 1])

################
### PARA GK ####
################

gc <- ts(Data_FM$GK_K / 4, start = c(2004, 01), frequency = 4)
yk <- ts(Data_FM$PIB_K / 4, start = c(2004, 01), frequency = 4)

gc <- window(gc, start = c(2007, 01), end = c(2019, 04))
yk <- window(yk, start = c(2007, 01), end = c(2019, 04))

source("LP.R")

rto <- gc / yk
rto <- lag1(rto, 1, 1, length(rto))

gc <- diff(log(gc))
yk <- diff(log(yk))

# Example
LP.irf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95)
LP.cirf <- irf.lp.iv(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, cumulative = TRUE)
LP.mult.adj <- irf.lp.iv.multiplier(y = yk, s = gc, Z = gc, py = 2, ps = 2, H = 4, a = 0.95, adj = rto)

# IRFs
LP.irf$pe
LP.cirf$pe
LP.mult.adj$pe

# Plot
plot(0:4, LP.mult.adj$pe[, 1], type = "l", col = "red", xlab = "", ylab = "", ylim = c(min(LP.mult.adj$lb[, 1]), max(LP.mult.adj$ub[, 1])))
lines(0:4, LP.mult.adj$ub[, 1])
lines(0:4, LP.mult.adj$lb[, 1])

