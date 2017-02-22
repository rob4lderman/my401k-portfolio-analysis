#
# Functions for Mean-Variance Portfolio Analysis 
# HI!
#
# TODO
#
# gmv.weights.vec <- MeanVar$getMinVarPortfolio( mu.vec, Sigma.mat )
# tan.weights.vec <- MeanVar$getTangencyPortfolio( mu.vec, Sigma.mat, risk.free.rate )
# eff.weights.vec <- MeanVar$getEfficientPortfolio( mu.vec, Sigma.mat, target.return )
# eff.weights.mat <- MeanVar$buildEfficientFrontier(mu.vec, Sigma.mat, alpha.vec=seq(from=-1,to=2,by=0.1), gmv.portfolio.weights=NULL, eff.portfolio.weights=NULL) 
#
# gmv.weights.vec <- MeanVar$getMinVarPortfolio.noShort( mu.vec, Sigma.mat )
# tan.weights.vec <- MeanVar$getTangencyPortfolio.noShort( mu.vec, Sigma.mat, risk.free.rate )
# eff.weights.vec <- MeanVar$getEfficientPortfolio.noShort( mu.vec, Sigma.mat, target.return )
# eff.weights.mat <- MeanVar$buildEfficientFrontier.noShort(mu.vec, Sigma.mat, target.return.vec )
#
# portfolioWeightsBarplot(x.mat, ...)
#

library(quadprog)   # solve.QP
library(Rsolnp)     # solnp


#
# List object holds all the functions
#
MeanVar <- list()


#
# CTOR for portfolio objects.
# 
# Portfolio objects are just a list of data elements:
#       port$asset.vec
#       port$mu.vec
#       port$stdev.vec
#       port$Sigma.mat
#       port$tan.rf.rate
#
#       port$w.vec
#       port$er
#       port$stdev
#
#       port$gmv.w.vec
#       port$gmv.er
#       port$gmv.stdev
#
#       port$tan.w.vec
#       port$tan.er
#       port$tan.stdev
#
#       port$eff.w.mat
#       port$eff.er.vec
#       port$eff.stdev.vec
#
# 
MeanVar$initPortfolio <- function( asset.vec,
                                   mu.vec,
                                   stdev.vec,
                                   Sigma.mat,
                                   tan.rf.rate ) {

    list(asset.vec=asset.vec,
         mu.vec=mu.vec,
         stdev.vec=stdev.vec,
         Sigma.mat=Sigma.mat,
         tan.rf.rate=tan.rf.rate)
}


#
# @param mu.vec: expected returns for all assets
# @param Sigma.mat: covariance matrix
#
# @return x.vec: a portfolio (vector of asset weights) that represents the
#                global minimum variance portfolio
#
#
# Constrained optimization: Lagrangian Multipliers
#
# A.mat %*% z.vec = b.vec
# z.vec = A.mat^-1 %*% b.vec
#
# A.mat = [ 2*Sigma  1 ]
#         [ 1'       0 ]
#
# b.vec = [ 0 ]
#         [ 0 ]
#         [ 1 ]
#
# z.vec = [ x.vec  ]
#         [ lambda ]
#
MeanVar$getMinVarPortfolio <- function( mu.vec, 
                                        Sigma.mat) {

    N <- length(mu.vec)

    A.mat.top = cbind( 2*Sigma.mat, rep(1,N) )
    A.mat = rbind( A.mat.top, c(rep(1,N),0) )

    b.vec = c(rep(0,N),1)

    z.vec = solve(A.mat) %*% b.vec
    x.vec = z.vec[1:N,1]
    x.vec
}



#
# Parameters:
#   mu.vec:     expected returns for all assets
#   Sigma.mat:  covariance matrix
#
# Return Value:
#   x.vec: vector of weights for each asset in the global-min-var portfolio
#
# Quadratic Programming problem with inequality constraints.
# Numerical solutions only! 
# Makes use of solve.QP
#
#       [ A_eq'  ]
# A' =  [ A_neq' ]
#
#       [ b_eq  ]
# b  =  [ b_neq ]
#
#                         [ 1  1  0  0 ]
# A' = [ 1'  ]        A = [ 1  0  1  0 ]
#      [ I_n ]            [ 1  0  0  1 ]
#
# b  = [   1    ]
#      [   0    ]
#      [   0    ]
#      [   0    ]
#
# D = 2*Sigma
#
# d = (0,...,0)'
#
MeanVar$getMinVarPortfolio.noShort <- function(mu.vec, Sigma.mat) {

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec = rep(0,N)

    A.mat.t <- rbind( rep(1,N),
                      diag(N) )
    A.mat <- t(A.mat.t)

    b.vec <- c(1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=1)       # 1 equality constraints (x'*1=1)
    x.vec <- qp.out$solution
    names(x.vec) <- names(mu.vec)
    x.vec
}


#
# @return portfolio object with gmv.w.vec, gmv.er, gmv.stdev
#
MeanVar$getMinVarPortfolioP.noShort <- function( portfolio ) {

    portfolio$gmv.w.vec <- MeanVar$getMinVarPortfolio.noShort(mu.vec = portfolio$mu.vec,
                                                              Sigma.mat = portfolio$Sigma.mat )

    portfolio$gmv.er <- ( t(portfolio$gmv.w.vec) %*% portfolio$mu.vec )[1]
    portfolio$gmv.stdev <- ( sqrt( t(portfolio$gmv.w.vec) %*% portfolio$Sigma.mat %*% portfolio$gmv.w.vec) )[1]

    portfolio
}


#
# @return portfolio with er and stdev
#
MeanVar$calcPortfolioMeanVar <- function( portfolio ) {

    portfolio$er <- ( t(portfolio$w.vec) %*% portfolio$mu.vec )[1]
    portfolio$stdev <- ( sqrt( t(portfolio$w.vec) %*% portfolio$Sigma.mat %*% portfolio$w.vec) )[1]

    portfolio
}


#
# Compute the tangency portfolio (aka the Sharpe Optimal Portfolio).
# The tangency portfolio is the point along the efficient frontier that
# is tangent to a line drawn from the point [0,risk.free.rate].
#
# The tangency portfolio is computed using constrained optimization
# to maximize sharpe ratio:
#
#               Sigma.mat^-1 %*% (mu.vec - risk.free.rate * 1)
#     t.vec = --------------------------------------------
#               1' %*% Sigma.mat^-1 %*% (mu.vec - risk.free.rate * 1)
#
#
# @param mu.vec: asset expected returns
# @param Sigma.mat: asset covariance matrix
# @param risk.free.rate: risk-free rate
#
# @return t.vec: tangency portfolio weights
#
#
MeanVar$getTangencyPortfolio <- function(mu.vec, 
                                         Sigma.mat, 
                                         risk.free.rate) {

    rf <- risk.free.rate
    one.vec <- rep(1,length(mu.vec))
    Sigma.mat.inv <- solve(Sigma.mat)

    tan.top <- Sigma.mat.inv %*% (mu.vec - rf * one.vec)
    tan.bot <- t(one.vec) %*% Sigma.mat.inv %*% (mu.vec - rf * one.vec) 

    t.vec =  tan.top[,1] / tan.bot
    t.vec
}


#
# Parameters:
#   mu.vec: asset expected returns
#   Sigma.mat: asset covariance matrix
#   risk.free.rate: risk-free rate
#
# Return Value:
#   t.vec: tangency portfolio weights
#
# Constrained optimization, maximize sharpe ratio:
#
#                 mu_p - rf
# max(x) SR_p = -------------
#                  sigma_p
# 
# mu_p = x' * mu
# 
# sigma_p^2 = x' * Sigma * x
# 
# s.t:
#       x' * 1 = 1
#
#       x_i >= 0
#
#         
#       [ A_eq'  ]
# A' =  [ A_neq' ]
#
#       [ b_eq  ]
# b  =  [ b_neq ]
#
#                                  [ mu_1-rf  1  0  0 ]
# A' = [ (mu_p - rf)' ]        A = [ mu_2-rf  0  1  0 ]
#      [ I_n          ]            [ mu_n-rf  0  0  1 ]
#
#      
# b  = [   1    ]       # equality constraint       (x'*1  = 1)
#      [   0    ]       # inequality constraints    (x_i >= 0 )
#
# I don't understand this solution...
#
MeanVar$getTangencyPortfolio.noShort <- function(mu.vec, Sigma.mat, risk.free.rate) {

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec <- rep(0, N)

    er.excess <- mu.vec - risk.free.rate
    A.mat <- cbind(er.excess, diag(1,N))
    b.vec <- c(1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=1)
    t.vec <- qp.out$solution / sum(qp.out$solution)
    names(t.vec) <- names(mu.vec)
    t.vec
}


#
#
# @param portfolio - must have tan.rf.rate defined
#
# @return portfolio with 
#         tan.w.vec containing tangency portfolio weights
#         tan.er, tan.stdev
#
MeanVar$getTangencyPortfolioP.noShort <- function(portfolio) {

    t.vec <- MeanVar$getTangencyPortfolio.noShort( portfolio$mu.vec,
                                                   portfolio$Sigma.mat,
                                                   portfolio$tan.rf.rate)

    portfolio$tan.w.vec <- t.vec
    portfolio$tan.er <- t(portfolio$tan.w.vec) %*% portfolio$mu.vec
    portfolio$tan.stdev <- sqrt( t(portfolio$tan.w.vec) %*% portfolio$Sigma.mat %*% portfolio$tan.w.vec )

    portfolio
}



#
# Parameters:
#   mu.vec: asset expected returns
#   Sigma.mat: asset covariance matrix
#   target.return: target expected return of the efficient portfolio
#
# Return Value:
#   x.vec: asset weights for the efficient portfolio
#
# Constrained optimization: minimize portfolio variance subject to constraints:
#   1) asset weights must sum to 1 
#   2) portfolio expected return must equal target return.
#
# min(xA,xB,xC) sigma2_p,x = x' * Sigma * x
#
#           s.t:
# 
#       mu_p,x = x' * mu = mu_p,t = target return
#
#       x'*1 = 1
#
# L(x,lambda1, lambda2) = x'*Sigma*x + 
#                         lambda1 * [ x'*mu - mu_p,t ] +
#                         lambda2 * [ x'*1 - 1 ]
#
# Take partial derivatives wrt x, lambda1, lambda2. 
# End up with 5 equations, 5 unknowns.
# 
# [ 2*Sigma  mu  1 ]   [    x    ]    [   0    ]
# [ mu'      0   0 ] * [ lambda1 ] =  [ mu_p,t ]
# [ 1'       0   0 ]   [ lambda2 ]    [   1    ]
# 
#         Ax         *     zx      =      b0
# 
# zx = Ax^1 * b0
#
# 
# TODO: PERF: separately solve for A.mat.inv and allow user to pass in
#
MeanVar$getEfficientPortfolio <- function(mu.vec, 
                                         Sigma.mat, 
                                         target.return) {

    N <- length(mu.vec)

    A.mat.top <- cbind( 2*Sigma.mat, mu.vec, rep(1,N) )
    A.mat.mid <- c( mu.vec, 0, 0 )
    A.mat.bot <- c( rep(1,N), 0, 0 )

    A.mat <- rbind(A.mat.top, 
                   A.mat.mid,
                   A.mat.bot)
    A.mat.inv <- solve(A.mat)

    b0.vec <- c(rep(0,N), target.return, 1)

    z.vec <- A.mat.inv %*% b0.vec
    x.vec <- z.vec[1:N,1]
    x.vec
}


#
# Parameters:
#   mu.vec: asset expected returns
#   Sigma.mat: asset covariance matrix
#   target.return: target expected return of the efficient portfolio
#
# Return Value:
#   x.vec: asset weights for the efficient portfolio, no shorting allowed
#
# Quadratic Programming problem with inequality constraints.
# Numerical solutions only! 
# Makes use of solve.QP
#
#        [ A_eq^T  ]
# A^T =  [ A_neq^T ]
#
#       [ b_eq  ]
# b  =  [ b_neq ]
#
#       [ mu' ]            [ mu_1  1  1  0  0 ]
# A^T = [ 1'  ]        A = [ mu_2  1  0  1  0 ]
#       [ I_n ]            [ mu_n  1  0  0  1 ]
#
#      [ mu_targ ]       # first equality constraint (x'*mu = target.return)
# b  = [   1     ]       # 2nd equality constraint   (x'*1  = 1)
#      [   0     ]       # inequality constraints    (x_i >= 0 )
#
# Constraint: A^T * x >= b
#
# D = 2*Sigma
#
# d = (0,...,0)'
#
MeanVar$getEfficientPortfolio.noShort <- function(mu.vec, Sigma.mat, target.return) {

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec = rep(0,N)

    A.mat.t <- rbind( t(mu.vec),
                      rep(1,N),
                      diag(N) )
    A.mat <- t(A.mat.t)

    b.vec <- c(target.return, 1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=2)       # 2 equality constraints
    x.vec <- qp.out$solution
    names(x.vec) <- names(mu.vec)
    x.vec
}


#
#
# All efficient frontier portfolios can be represented as a convex ("sum-to-one") 
# combination of two other eff frontier portfolios.
#
# @param mu.vec: asset expected returns
# @param Sigma.mat: asset covariance matrix
# @param alpha.vec: sequence of alpha values (eff frontier = convex combo of two frontier portfolios)
# @param gmv.portfolio.weights: weights for the global min var portfolio. can be null
# @param eff.portfolio.weights: weights for another eff portfolio. can be null
#
# @return x.mat: matrix of asset weights for the efficient portfolio
#                each column vector is a portfolio             
#                each row is an asset in the portfolio
#
# Computing means and variances for the portfolios in x.mat:
#
# eff.frontier.means <- apply(x.mat, 2, function(x.vec) { t(x.vec) %*% mu.vec })
# eff.frontier.means <- t(x.mat) %*% mu.vec 
#
# eff.frontier.sigmas <- apply(x.mat, 2, function(x.vec) { sqrt(t(x.vec) %*% Sigma.mat %*% x.vec) })
#
# eff.frontier.means <- apply(x.mat,
#                             1,
#                             function(x.vec) { t(x.vec) %*% mu.vec })
#
MeanVar$buildEfficientFrontier <- function(mu.vec, 
                                          Sigma.mat, 
                                          alpha.vec=seq(from=-1,to=2,by=0.1),
                                          gmv.portfolio.weights=NULL,
                                          eff.portfolio.weights=NULL) {

    if (is.null(gmv.portfolio.weights)) {
        gmv.portfolio.weights = MeanVar$getMinVarPortfolio(mu.vec, Sigma.mat)
    }

    if (is.null(eff.portfolio.weights)) {
        eff.portfolio.weights = MeanVar$getEfficientPortfolio(mu.vec, Sigma.mat, target.return=max(mu.vec))
    }

    x.mat = sapply( as.array(alpha.vec),
                    function(alpha) { alpha * gmv.portfolio.weights + (1 - alpha) * eff.portfolio.weights }
                  )
    x.mat
}


#
# @param mu.vec: asset expected returns
# @param Sigma.mat: asset covariance matrix
# @param target.return.vec: sequence of target returns. Typically between gmv portfolio return and max individual asset return
#
# @return x.mat: matrix of asset weights for the efficient portfolio
#                each column vector is a portfolio             
#                each row is an asset in the portfolio
#
#
# Quadratic Programming problem with inequality constraints.
# Numerical solutions only! 
# Makes use of solve.QP
#
MeanVar$buildEfficientFrontier.noShort <- function(mu.vec, Sigma.mat, target.return.vec) {

    x.mat = sapply( as.array(target.return.vec),
                    function(target.return) { MeanVar$getEfficientPortfolio.noShort(mu.vec, Sigma.mat, target.return) }
                  )
    x.mat
}


#
#
# Generate EFF portfolios from the GMV mean/var to max(asset-mean)
#
# @param gmv.portfolio - a portfolio with gmv.er filled in (global min var portfolio)
#
# @return portfolio with eff.w.mat added
#                   each column vector in eff.w.mat is a set of portfolio weights
#                   each row is an asset in the portfolio
#
#
MeanVar$buildEfficientFrontierP.noShort <- function(gmv.portfolio) {

    target.return.vec <- seq(from=gmv.portfolio$gmv.er,
                             to=max(gmv.portfolio$mu.vec),
                             length.out=16)

    #
    # Note: having issues with the final target return.
    #       Removing it for now
    #
    w.mat <- MeanVar$buildEfficientFrontier.noShort(gmv.portfolio$mu.vec,
                                                    gmv.portfolio$Sigma.mat,
                                                    target.return.vec[-16]) 

    gmv.portfolio$eff.w.mat <- w.mat
    rownames(gmv.portfolio$eff.w.mat) <- gmv.portfolio$asset.vec

    #
    # Compute portfolio exp. return and stdev
    #
    gmv.portfolio$eff.er.vec <- apply(gmv.portfolio$eff.w.mat,
                                      2,   # cols
                                      function(w.vec) { t(w.vec) %*% gmv.portfolio$mu.vec } )

    gmv.portfolio$eff.stdev.vec <- apply(gmv.portfolio$eff.w.mat,
                                         2, 
                                         function(w.vec) { sqrt( t(w.vec) %*% gmv.portfolio$Sigma.mat %*% w.vec ) } )
    gmv.portfolio
}


#
# Plot the efficient frontier
#
# @param eff.portfolio - must contain eff.er.vec, eff.stdev.vec
#
MeanVar$plotEfficientFrontierP <- function(eff.portfolio, ...) {

    plot(x=eff.portfolio$eff.stdev.vec,
         y=eff.portfolio$eff.er.vec,
         pch=16, 
         col=1:length(eff.portfolio$eff.er.vec),
         xlim=c(0, max(eff.portfolio$eff.stdev.vec) *1.10), 
         ylim=c(0, max(eff.portfolio$eff.er.vec) * 1.10), 
         xlab=expression(sigma), 
         ylab=expression(mu),
         ...)

}

#
# Plot the efficient frontier
#
# @param eff.portfolio - must contain eff.er.vec, eff.stdev.vec
#
MeanVar$pointsEfficientFrontierP <- function(eff.portfolio, ...) {

    points(x=eff.portfolio$eff.stdev.vec,
           y=eff.portfolio$eff.er.vec,
           pch=16, 
           type="b",
           col="darkgreen",
           ...)
}


#
# Flesh out portfolio stats:
#   a) GMV portfolio
#   b) Efficient Frontier
#   c) Tangency portfolio
#
# @return portfolio object TODO
# 
MeanVar$buildPortfolio.noShort <- function(portfolio) {

    gmv.portfolio <- MeanVar$getMinVarPortfolioP.noShort( portfolio )

    eff.portfolio <- MeanVar$buildEfficientFrontierP.noShort( gmv.portfolio )

    portfolio <- MeanVar$getTangencyPortfolioP.noShort( eff.portfolio )

    portfolio
}


#
#
# Plot:
#   1. the efficient frontier
#   2. tangency portfolio
#   3. ?
#
MeanVar$plotPortfolio <- function(portfolio, ...) {

    #
    # Mean-Var chart 
    #
    MeanVar$meanVarPlot( portfolio$stdev.vec,
                         portfolio$mu.vec,
                         portfolio$asset.vec,
                         main="Mean-Var Plot")

    MeanVar$pointsEfficientFrontierP(portfolio, ...)


    #
    # Add tangency portfolio
    #
    points(x=c(0, portfolio$tan.stdev),
           y=c(portfolio$tan.rf.rate, portfolio$tan.er),
           type="b",
           pch=16,
           col="blue")

    text(x=c(0, portfolio$tan.stdev),
         y=c(portfolio$tan.rf.rate, portfolio$tan.er),
         labels=c(paste("rf=",round(portfolio$tan.rf.rate,6),sep=""),"TAN"),
         pos=4)

}


#
#
# Plot eff portfolio asset weights
#
#
MeanVar$plotEffAssetWeights <- function(portfolio, ...) {

    matplot( t(portfolio$eff.w.mat), 
             type="p", 
             pch=16, 
             col=1:length(portfolio$asset.vec),
             xlab="Portfolios",
             ylab="Asset weights",
             ...)
    legend("top", legend=portfolio$asset.vec, col=1:length(portfolio$asset.vec), pch=16)
}


#
# List out all efficient portfolio weights, er, and stdev
#
MeanVar$printEffPortfolios <- function( portfolio ) {
    rbind( round(portfolio$eff.w.mat, 2), 
           er=round(portfolio$eff.er.vec,4), 
           stdev=round(portfolio$eff.stdev.vec, 4) )
}



#
# Generate a mean-var plot for the given sigma and mu data.
#
# 
MeanVar$meanVarPlot <- function( sigma.vec,
                                 mu.vec,
                                 labels.vec, ... ) {

    plot(x=sigma.vec,
         y=mu.vec,
         pch=16, 
         col=1:length(sigma.vec),
         xlim=c(0, max(sigma.vec) *1.10), 
         ylim=c(0, max(mu.vec) * 1.10), 
         xlab=expression(sigma), 
         ylab=expression(mu),
         ...)

    text(x=sigma.vec,
         y=mu.vec,
         labels=labels.vec,
         pos=4,
         col=1:length(sigma.vec),
         cex=0.75)
}



portfolioWeightsBarplot <- function(x.mat, ...) {
    # Parameters:
    #   x.mat: N x m matrix of portfolio weights
    #          N = number of portfolios
    #          m = number of weights (assets) per portfolio
    # 
    # Generates barplot of portfolio weights
    #

    test1 <- test2 <- t(x.mat)
    test1[test1<0] <- 0
    test2[test2>0] <- 0
    myrange <- c(min(colSums(test2)),max(colSums(test1)))
    barplot(test1,ylim=myrange,legend.text=T,...)
    barplot(test2,add=TRUE,ylim=rev(myrange),...)
}



#
#
# 1. Downloads SBUX,MSFT,IBM price data from yahoo
# 2. Computes monthly returns
# 3. Plots assets on a mean-var plot
# 4. Plots global min-var portfolio
# 5. Computes and plots efficient portfolios that match asset returns 
# 6. Plots the full efficient frontier
# 7. Plots the tangency portfolio
# 8. Plots RF + tangency line (what's it called again??)
# 9. Plots min-var, tangency, efficient frontier with NO SHORTING
# 10. Stacked bar chart of eff portfolios (asset weights)
#
#
my.portfolio.test <- function() {


    library(tseries)    # get.hist.quote
    library(zoo)        # coredata
    library(quadprog)   # solve.QP

    #
    # Load price data from yahoo
    #
    SBUX_prices <- get.hist.quote(instrument="sbux", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)
    MSFT_prices <- get.hist.quote(instrument="msft", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)
    IBM_prices <-  get.hist.quote(instrument="ibm", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)

    #
    # Compute simple returns, means, sd, cov
    # Portfolio theory assumes simple returns (as opposed to cc returns)
    # 
    all_prices <- merge(IBM_prices, MSFT_prices, SBUX_prices)

    # diff: computes pt1 - pt0
    # lag: shifts all_prices by k=-1 (so that pt-1 -> pt)
    simple_returns <- diff(all_prices) / lag(all_prices,k=-1)
    simple_returns <- coredata(simple_returns)

    asset_names <- c("IBM", "MSFT", "SBUX")
    colnames(simple_returns) <- asset_names

    mu.vec <- apply(simple_returns, 2, mean)
    sigma.vec <- apply(simple_returns, 2, sd)
    sigma.mat <- cov(simple_returns)


    par(mfrow=c(3,1))

    #
    # Mean-Var Plot: Plot and label the assets 
    #
    plot(x=sigma.vec,
         y=mu.vec,
         pch=16, 
         ylim=c(0, max(mu.vec) * 1.5), 
         xlim=c(0, max(sigma.vec) *1.5), 
         xlab=expression(sigma[p]), 
         ylab=expression(mu[p]))


    # ..OR..:
    text(x=sigma.vec,
         y=mu.vec,
         labels=names(sigma.vec),
         pos=4)

    #
    # Compute global min var portfolio and plot it.
    #
    x_p.gmv <- MeanVar$getMinVarPortfolio(mu.vec=mu.vec,
                                         Sigma.mat=sigma.mat)
    mu_p.gmv <- t(x_p.gmv) %*% mu.vec
    sigma_p.gmv <- sqrt( t(x_p.gmv) %*% sigma.mat %*% x_p.gmv )

    points(x=sigma_p.gmv,
           y=mu_p.gmv,
           pch=16,
           col="orange")
    text(x=sigma_p.gmv,
         y=mu_p.gmv,
         labels="GMV",
         pos=2)


    #
    # Compute efficient portfolio w/ same return as individual asset returns, 
    # and plot them
    #
    for (i in seq_along(asset_names)) {

        x_p.eff <- MeanVar$getEfficientPortfolio(mu.vec=mu.vec, 
                                                Sigma.mat=sigma.mat, 
                                                target.return=mu.vec[i])
        mu_p.eff <- t(x_p.eff) %*% mu.vec
        sigma_p.eff <- sqrt( t(x_p.eff) %*% sigma.mat %*% x_p.eff)

        # x_p.eff
        # mu_p.eff
        # sigma_p.eff

        points(x=sigma_p.eff,
               y=mu_p.eff,
               pch=16,
               col="green")
        text(x=sigma_p.eff,
             y=mu_p.eff,
             labels=paste("EFF-",asset_names[i]),
             pos=2)
    }


    #
    # Now we want to plot the efficient frontier
    #
    eff.frontier.weights <- MeanVar$buildEfficientFrontier(mu.vec=mu.vec,
                                                          Sigma.mat=sigma.mat,
                                                          gmv.portfolio.weights=x_p.gmv)
    eff.frontier.means <- apply(eff.frontier.weights,
                                2,
                                function(x.vec) { t(x.vec) %*% mu.vec })

    eff.frontier.sigmas <- apply(eff.frontier.weights,
                                 2,
                                 function(x.vec) { sqrt(t(x.vec) %*% sigma.mat %*% x.vec) })
    
    points(x=eff.frontier.sigmas,
           y=eff.frontier.means,
           type="b",
           pch=16,
           col="green")


    #
    # Highlight the tangency portfolio
    # 
    rf <- 0.03/12   # monthly
    x_p.tan <- MeanVar$getTangencyPortfolio(mu.vec=mu.vec,
                                           Sigma.mat=sigma.mat,
                                           risk.free.rate=rf)
    mu_p.tan <- t(x_p.tan) %*% mu.vec
    sigma_p.tan <- sqrt( t(x_p.tan) %*% sigma.mat %*% x_p.tan )

    points(x=c(0, sigma_p.tan),
           y=c(rf, mu_p.tan),
           pch=16,
           cex=2,
           col="orange")
    text(x=c(0, sigma_p.tan),
         y=c(rf, mu_p.tan),
         labels=c(paste("rf=",rf),"TAN"),
         pos=4)


    #
    # Compute and plot tangency + risk-free portfolio combinations
    #
    tan.portfolio.alpha <- seq(from=0,to=2,by=0.1)
    tan.frontier.means <- rf + tan.portfolio.alpha * ( mu_p.tan - rf )
    tan.frontier.sigmas <- tan.portfolio.alpha * sigma_p.tan

    points(x=tan.frontier.sigmas,
           y=tan.frontier.means,
           pch=16,
           col="orange")


    #
    # NO SHORTING! Global Min Var Portfolio
    #
    x_p.gmv.ns <- MeanVar$getMinVarPortfolio.noShort(mu.vec, sigma.mat)
    mu_p.gmv.ns <- t(x_p.gmv.ns) %*% mu.vec
    sigma_p.gmv.ns <- sqrt(  t(x_p.gmv.ns) %*% sigma.mat %*% x_p.gmv.ns )

    points(x=sigma_p.gmv.ns,
           y=mu_p.gmv.ns,
           pch=16,
           col="red")
    text(x=sigma_p.gmv.ns,
         y=mu_p.gmv.ns,
         labels="GMV-NS",
         pos=4)

    #
    # NO SHORTING! Highlight the tangency portfolio
    # 
    rf <- 0.03/12   # monthly
    x_p.tan.ns <- MeanVar$getTangencyPortfolio.noShort(mu.vec=mu.vec,
                                                      Sigma.mat=sigma.mat,
                                                      risk.free.rate=rf)
    mu_p.tan.ns <- t(x_p.tan.ns) %*% mu.vec
    sigma_p.tan.ns <- sqrt( t(x_p.tan.ns) %*% sigma.mat %*% x_p.tan.ns )

    points(x=c(0, sigma_p.tan.ns),
           y=c(rf, mu_p.tan.ns),
           pch=16,
           cex=2,
           col="red")
    text(x=c(0, sigma_p.tan.ns),
         y=c(rf, mu_p.tan.ns),
         labels=c(paste("rf=",rf),"TAN-NS"),
         pos=4)


    #
    # NO SHORTING! Eff frontier from GMV portfolio to MAX(mu.vec)
    # 
    target.return.vec <- seq(from=mu_p.gmv.ns, to=max(mu.vec), length.out=10)

    eff.frontier.ns.weights <- MeanVar$buildEfficientFrontier.noShort( mu.vec=mu.vec,
                                                                      Sigma.mat=sigma.mat,
                                                                      target.return.vec=target.return.vec)

    eff.frontier.ns.means = apply(eff.frontier.ns.weights,
                                  2,
                                  function(x.vec) { t(x.vec) %*% mu.vec } )

    eff.frontier.ns.sigmas = apply(eff.frontier.ns.weights,
                                   2,
                                   function(x.vec) { sqrt( t(x.vec) %*% sigma.mat %*% x.vec ) } )
    points(x=eff.frontier.ns.sigmas,
           y=eff.frontier.ns.means,
           type="b",
           pch=16,
           col="red")


    #
    # Plot efficient frontier portfolio weights in a stacked bar plot.
    #
    portfolioWeightsBarplot( eff.frontier.weights, beside=T)


    #
    # Plot efficient tangent + risk-free portfolio weights in a stacked bar plot.
    #
    tan.frontier.weights <- t(sapply(tan.portfolio.alpha,function(alpha) { alpha * x_p.tan} ))

    # add risk-free weight
    tan.frontier.weights = cbind(1-tan.portfolio.alpha, tan.frontier.weights  )
    colnames(tan.frontier.weights)[1] = "RF"

    portfolioWeightsBarplot( tan.frontier.weights, beside=T)

    par(mfrow=c(1,1))

}


#
# TODO
#
# @return list(x=x.eff.vec, mu=mu.eff, sigma=sigma.eff, sol=sol)
#
MeanVar$getEfficientPortfolio.targetRisk <- function(mu.vec,
                                                    Sigma.mat,
                                                    target.risk) {

    #
    # ----- solnp: Parameters
    #
    #     pars: The starting parameter vector.
    # 
    #      fun: The main function which takes as first argument the parameter
    #           vector and returns a single value.
    # 
    #    eqfun: (Optional) The equality constraint function returning the
    #           vector of evaluated equality constraints.
    # 
    #      eqB: (Optional) The equality constraints.
    # 
    #  ineqfun: (Optional) The inequality constraint function returning the
    #           vector of evaluated inequality constraints.
    # 
    #   ineqLB: (Optional) The lower bound of the inequality constraints.
    # 
    #   ineqUB: (Optional) The upper bound of the inequality constraints.
    # 
    #       LB: (Optional) The lower bound on the parameters.
    # 
    #       UB: (Optional) The upper bound on the parameters.
    # 
    #
    # -----  Return Value: List containing the following values:
    #
    #     pars: Optimal Parameters.
    # 
    # convergence : Indicates whether the solver has converged (0) or not (1
    #           or 2).
    # 
    #   values: Vector of function values during optimization with last one
    #           the value at the optimal.
    # 
    # lagrange: The vector of Lagrange multipliers.
    # 
    #  hessian: The Hessian of the augmented problem at the optimal solution.
    # 
    #   ineqx0: The estimated optimal inequality vector of slack variables
    #           used for transforming the inequality into an equality
    #           constraint.
    # 
    # nfuneval: The number of function evaluations.
    # 
    #  elapsed: Time taken to compute solution.
    #

    #
    # Use solnp to compute an efficient portfolio for the given
    # level of risk.
    #
    #

    d <- length(mu.vec)
    x.vec <- rep(1/d,d)
    sigma.p <- target.risk

    sol <- solnp( pars=x.vec,
                  fun=function(x.vec) { -1 * t(x.vec) %*% mu.vec },     # (-1) because solnp minimizes this function
                  eqfun=function(x.vec) { sum(x.vec) },
                  eqB=c(1),
                  ineqfun=function(x.vec) { sqrt( t(x.vec) %*% Sigma.mat %*% x.vec ) },
                  ineqUB=sigma.p,
                  ineqLB=sigma.p-0.001,
                  control=list(trace=0))
    # warnings()

    x.eff.vec <- sol$pars
    mu.eff <- t(x.eff.vec) %*% mu.vec
    sigma.eff <- sqrt( t(x.eff.vec) %*% Sigma.mat %*% x.eff.vec )

    list(x=x.eff.vec, mu=mu.eff, sigma=sigma.eff, sol=sol)
}

