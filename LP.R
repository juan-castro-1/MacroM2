###################
# FUNCTION INPUTS #
###################

#   y:           endogenous variable.
#   s:           shock variable(s).
#   Z:           instrumental variables.
#   W:           exogenous variables. Default is NULL.
#  py:           number of lags of endogenous variable. Default is zero.
#  ps:           number of lags of shock variable(s). Default is zero.
#  pw:           number of lags of exogenous variables. Default is zero.
#   H:           horizon.
#   a:           confidence level.
# adj:           shock adjustment series. Default is NULL.
# cumulative:    whether or not cumulative IRFs. Default is FALSE.


# Auxiliary Functions

lag1 <- function(X, m, p, N) {
  L1 <- matrix(NA, N - p, m * p)
  for (l in 1:p) {
    L1[, (m * (l - 1) + 1):(m * l)] <- as.matrix(X)[(p + 1 - l):(N - l), ]
  }
  return(L1)
}

rhs.y <- function(y, py, pmax, N, h) {
  if (py == 0) {
    Y1 <- NULL
  } else {
    Y1 <- lag1(y, 1, h + pmax, N)[, (h + 1):(h + py)]
  }
  return(Y1)
}

rhs.s <- function(s, ms, ps, pmax, N, h) {
  if (ps == 0) {
    S1 <- NULL
  } else {
    S1 <- lag1(s, ms, h + pmax, N)[, (h + 1):(ms * (h + ps))]
  }
  return(S1)
}

rhs.W <- function(W, mw, pw, pmax, N, h) {
  if (is.null(W) == TRUE) {
    Wt <- NULL
  } else {
    if (pw == 0) {
      Wt <- W[(pmax + 1):(N - h), ]
    } else {
      W0 <- W[(pmax + 1):(N - h), ]
      W1 <- lag1(W, mw, h + pmax, N)[, (h + 1):(mw * (h + pw))]
      Wt <- cbind(W0, W1)
    }
  }
  return(Wt)
}

rhs <- function(y, s, W, ms, mw, py, ps, pw, pmax, N, h) {
  Y1 <- rhs.y(y, py, pmax, N, h)
  S1 <- rhs.s(s, ms, ps, pmax, N, h)
  Wt <- rhs.W(W, mw, pw, pmax, N, h)
  return(cbind(Y1, S1, Wt))
}

OLS <- function(y, X) {
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  return(as.vector(solve(XtX, Xty)))
}

OLS.EHW <- function(y, X, j) {
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  b <- as.vector(solve(XtX, Xty))
  u <- as.vector(y - X %*% b)
  iXtX <- chol2inv(chol(XtX))
  XtOX <- t(X) %*% diag(u ^ 2) %*% X
  VEHW <- iXtX %*% XtOX %*% iXtX
  return(list(pe = b[1:j], se = sqrt(diag(VEHW)[1:j])))
}

# IRF by Local Projections, cumulative or not

irf.lp.iv <- function(y, s, Z, W = NULL, py = 0, ps = 0, pw = 0, H, a, cumulative = FALSE) {
  #
  y <- as.vector(y)
  s <- as.matrix(s)
  Z <- as.matrix(Z)
  if (is.null(W) == FALSE) { W <- as.matrix(W) }
  #
  pmax <- max(py, ps, pw)
  #
  N <- length(y)
  #
  ms <- ncol(s)
  if (is.null(W) == FALSE) { mw <- ncol(W) }
  #
  irf.pe <- array(NA, c(H + 1, ms))
  irf.ub <- array(NA, c(H + 1, ms))
  irf.lb <- array(NA, c(H + 1, ms))
  #
  z <- qnorm(1 - (1 - a) / 2)
  for (h in 0:H) {
    st <- s[(pmax + 1):(N - h), ]
    Xt <- rhs(y, s, W, ms, mw, py, ps, pw, pmax, N, h)
    Zt <- Z[(pmax + 1):(N - h), ]
    if ((cumulative == TRUE) & (h > 0)) {
      st <- tail(st, -h)
      Xt <- tail(Xt, -h)
      Zt <- tail(Zt, -h)
      yh <- diff(cumsum(c(0, y[(pmax + h + 1):N])), h + 1)
    } else {
      yh <- y[(pmax + h + 1):N]
    }
    # IV Auxiliary Regression
    Pt <- cbind(Zt, Xt)
    alpha <- OLS(st, Pt)
    st <- as.matrix(Pt %*% alpha)
    # IRF Estimation
    Xt <- cbind(st, Xt)
    beta <- OLS.EHW(yh, Xt, ms)
    irf.pe[h + 1, ] <- beta$pe
    irf.ub[h + 1, ] <- irf.pe[h + 1, ] + z * beta$se
    irf.lb[h + 1, ] <- irf.pe[h + 1, ] - z * beta$se
  }
  rownames(irf.pe) <- 0:H
  rownames(irf.ub) <- 0:H
  rownames(irf.lb) <- 0:H
  colnames(irf.pe) <- paste("Rsp.", "S", 1:ms, sep = "")
  colnames(irf.ub) <- paste("Rsp.", "S", 1:ms, sep = "")
  colnames(irf.lb) <- paste("Rsp.", "S", 1:ms, sep = "")
  return(list(ub = irf.ub, pe = irf.pe, lb = irf.lb))
}

# Ratio of cumulative IRF by Local Projections, with the possibility of adjustment

irf.lp.iv.multiplier <- function(y, s, Z, W = NULL, py = 0, ps = 0, pw = 0, H, a, adj = NULL) {
  #
  y <- as.vector(y)
  s <- as.matrix(s)
  Z <- as.matrix(Z)
  if (is.null(W) == FALSE) { W <- as.matrix(W) }
  if (is.null(adj) == FALSE) { adj <- as.vector(adj) }
  #
  pmax <- max(py, ps, pw)
  #
  N <- length(y)
  #
  ms <- ncol(s)
  if (is.null(W) == FALSE) { mw <- ncol(W) }
  #
  irf.pe <- array(NA, c(H + 1, ms))
  irf.ub <- array(NA, c(H + 1, ms))
  irf.lb <- array(NA, c(H + 1, ms))
  #
  z <- qnorm(1 - (1 - a) / 2)
  for (h in 0:H) {
    Xt <- rhs(y, s, W, ms, mw, py, ps, pw, pmax, N, h)
    Zt <- Z[(pmax + 1):(N - h), ]
    if (is.null(adj) == FALSE) { adjt <- adj[(pmax + 1):(N - h)] }
    if (h > 0) {
      Xt <- tail(Xt, -h)
      Zt <- tail(Zt, -h)
      if (is.null(adj) == FALSE) { adjt <- tail(adjt, -h) }
      if (ms > 1) {
        sh <- diff(rbind(0, apply(s[(pmax + h + 1):N, ], 2, cumsum)), h + 1)
      } else {
        sh <- diff(c(0, cumsum(s[(pmax + h + 1):N, ])), h + 1)
      }
      yh <- diff(cumsum(c(0, y[(pmax + h + 1):N])), h + 1)
    } else {
      sh <- s[(pmax + 1):N, ]
      yh <- y[(pmax + 1):N]
    }
    if (is.null(adj) == FALSE) { sh <- sh * adjt }
    # IV Auxiliary Regression
    Pt <- cbind(Zt, Xt)
    alpha <- OLS(sh, Pt)
    sh <- as.matrix(Pt %*% alpha)
    # Multiplier Estimation
    Xt <- cbind(sh, Xt)
    beta <- OLS.EHW(yh, Xt, ms)
    irf.pe[h + 1, ] <- beta$pe
    irf.ub[h + 1, ] <- irf.pe[h + 1, ] + z * beta$se
    irf.lb[h + 1, ] <- irf.pe[h + 1, ] - z * beta$se
  }
  rownames(irf.pe) <- 0:H
  rownames(irf.ub) <- 0:H
  rownames(irf.lb) <- 0:H
  colnames(irf.pe) <- paste("Rsp.", "S", 1:ms, sep = "")
  colnames(irf.ub) <- paste("Rsp.", "S", 1:ms, sep = "")
  colnames(irf.lb) <- paste("Rsp.", "S", 1:ms, sep = "")
  return(list(ub = irf.ub, pe = irf.pe, lb = irf.lb))
}