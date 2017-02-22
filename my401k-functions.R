#
# @rob4lderman
# aug2016
# 
#


library(tseries)    # get.hist.quote
library(zoo)        # coredata
library(quadprog)   # solve.QP
library(dplyr)      # mutate, arrange, filter

source("mean-var.R")

My401k <- list()


#
# TODO: why simple returns and not CC returns?  check comp finance notes
#
# compute simple returns 
# simple.return: (p_t - p_t-1) / p_t-1
# diff: computes p_t - p_t-1
# lag: shifts by k=-1 (so that p_t = p_t-1)
#
# @param quotes.zoo: zoo object containing all merged ticker quotes
#
# @return zoo object of simple returns
#
My401k$computeSimpleReturns <- function(quotes.zoo) {

    simple.returns.zoo <- diff(quotes.zoo) / stats::lag(quotes.zoo,k=-1)

    simple.returns.zoo <- na.omit(simple.returns.zoo)       

    simple.returns.zoo
}


#
# Get return data with fees subtracted
#
# @param simple.returns.zoo - as returned by computeSimpleReturns
# @param funds.df - funds list, contains $fees column
#
# @return zoo object of simple returns w/ fees subtracted
#
My401k$computeSimpleReturnsWithFees <- function(simple.returns.zoo,
                                                funds.df) {

    fees <- funds.df[ colnames(simple.returns.zoo),]$fees

    num.funds <- ncol(simple.returns.zoo)
    num.months <- nrow(simple.returns.zoo)

    # 
    # Expand fees vector into a matrix (every column vector is
    # just the associated fund's fee repeated in every row)
    #
    fees.mat <- matrix( rep(fees, each=num.months), 
                        nrow=num.months,
                        ncol=num.funds )

    simple.returns.zoo - fees.mat / 100 / 12
}


#
# Compute mean, weighted mean, stdev, and cumulative return for all funds
# in given simple.returns.zoo.
#
# @param simple.returns.zoo - zoo object with simple monthly returns
#
# @return data.frame(mean=, 
#                    weighted.mean=,
#                    cum.return=,
#                    sd=)
# 
My401k$computeSimpleReturnsStats <- function(simple.returns.zoo) {

    returns.mean <- apply(simple.returns.zoo, 2, mean)

    num.months <- nrow( simple.returns.zoo )
    weights <- 1:num.months / num.months        # higher weights for later (more recent) months
    returns.weighted.mean <- apply(simple.returns.zoo, 2, weighted.mean, w=weights)

    returns.sd <- apply(simple.returns.zoo, 2, sd)

    returns.cum <- apply( 1 + simple.returns.zoo, 2, cumprod )[num.months,]   # select only the final entry

    data.frame(ticker=colnames(simple.returns.zoo),
               mean=returns.mean,
               sd=returns.sd,
               weighted.mean=returns.weighted.mean,
               cum.return=returns.cum)
}


#
# @param cluster.tickers.vec
#
# @return a zoo object containing timeseries data for the given cluster.tickers.vec only.
#
My401k$filterReturnsForCluster <- function( simple.returns.zoo,
                                            cluster.tickers.vec ) {

    cluster.returns.zoo <- simple.returns.zoo[ , cluster.tickers.vec ]

    #
    # convert to data frame in case cluster.tickers has only
    # 1 entry (which causes zoo to return a vector instead of
    # a data.frame/matrix).
    #
    df <- data.frame(cluster.returns.zoo)
    colnames(df) <- cluster.tickers.vec
    zoo(df, order.by=index(cluster.returns.zoo))
}


#
# Generate cumulative returns charts for all clusters.
# 
My401k$chartCumulativeReturns <- function( simple.returns.zoo,
                                           simple.returns.stats.df,
                                           cluster.tickers.list,
                                           cluster.funds.list ) {

    num.clusters <- length(cluster.tickers.list)

    for (cluster.index in 1:num.clusters) {

        cluster.tickers.vec <- cluster.tickers.list[[cluster.index]]

        #
        # Filter for the returns in this cluster and chart them
        #
        cluster.returns.zoo <- My401k$filterReturnsForCluster( simple.returns.zoo, cluster.tickers.vec )

        #
        # Mean-Var chart of the funds in the cluster
        #
        cluster.returns.stats.df <- subset(simple.returns.stats.df, ticker %in% cluster.tickers.vec)

        #
        # List the funds
        #
        cluster.funds.df <- cluster.funds.list[[cluster.index]]

        My401k$generateReturnsCharts( cluster.returns.zoo,
                                      cluster.returns.stats.df,
                                      cluster.funds.df,
                                      cumreturns.main = paste("Cumulative Returns for cluster: ", cluster.index, "/", num.clusters, sep="") )
    }
}


#
# Generate:
#   a) cumulative return chart
#   b) mean-var chart
#   c) list of funds
#
My401k$generateReturnsCharts <- function(simple.returns.zoo,
                                         simple.returns.stats.df,
                                         funds.df,
                                         cumreturns.main = "Cumulative Returns",
                                         ...) {

    My401k$plotCumReturnsAndMeanVar(simple.returns.zoo,
                                    simple.returns.stats.df,
                                    cumreturns.main,
                                    ...);
    #
    # List the funds
    #
    funds.df <- merge(funds.df, simple.returns.stats.df, by="ticker", all.x=T)
    print( arrange( funds.df, desc(cum.return) ) )

}


#
#
# @param simple.returns.zoo - zoo object containing monthly returns
# @param simple.returns.stats.df - data frame containing $sd, $mean, and $ticker 
#
# Generate:
#   a) cumulative return chart
#   b) mean-var chart
#
My401k$plotCumReturnsAndMeanVar<- function(simple.returns.zoo,
                                           simple.returns.stats.df,
                                           cumreturns.main = "Cumulative Returns",
                                           ...) {
    #
    # Cumulative returns chart
    #
    chart.CumReturns( simple.returns.zoo,
                      wealth.index=T, 
                      colorset = 1:ncol(simple.returns.zoo),   # show all
                      legend.loc="topleft",
                      main=cumreturns.main,
                      cex.legend=0.65,
                      ...
                    )

    #
    # Mean-Var chart 
    #
    MeanVar$meanVarPlot( simple.returns.stats.df$sd,
                         simple.returns.stats.df$mean,
                         simple.returns.stats.df$ticker,
                         main="Mean-Variance Plot")

}



