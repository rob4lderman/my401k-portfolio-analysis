#
# @rob4lderman
# feb 2017
#
#

library(tseries)                # get.hist.quote
library(zoo)                    # coredata
library(dplyr)                  # mutate, arrange, filter


My401kEtl <- list()


#
# Load and clean fund list
#
# @return data frame containing fund list
#
My401kEtl$loadFundList <- function() {

    funds <- read.csv("funds2.csv")

    #
    # add fod field (foreign vs domestic):
    # note: looks like I already did this and wrote it back to funds2.csv
    #
    # fod <- ifelse( grepl("foreign|emer|world|int'l|international|global|asia|europe|latin", funds2$category, ignore.case=T), "Foreign", "Domestic" )
    # funds2 <- mutate(funds2, fod = as.factor(fod))
    #

    funds <- funds %>% 
                filter(ticker != "TODO") %>%
                mutate(fees = as.numeric( gsub( "%", "", as.character(gross.expense.ratio) ) ) ) %>%
                select(name,ticker,asset.class,category,fod,fees)

    rownames(funds) <- funds$ticker

    funds
}

#
# Download ticker quotes from yahoo.
#
# @return zoo object containing monthly quotes for the given ticker symbol
#
My401kEtl$downloadTickerQuotes <- function(symbol, 
                                        start="2001-01-01",
                                        end="2017-12-31", 
                                        quote="AdjClose",
                                        provider="yahoo", 
                                        origin="1970-01-01",
                                        compression="m", 
                                        retclass="zoo", 
                                        quiet = TRUE) {

    quote <- tryCatch( { 
            get.hist.quote(instrument=symbol,
                           start=start,
                           end=end,
                           quote=quote,
                           provider=provider,
                           origin=origin,
                           compression=compression,
                           retclass=retclass,
                           quiet = quiet)
        },
        warning = function(war) {
            print(war)
            return( zoo() )
        },
        error = function(err) {
            print(err)
            return( zoo() )
        })

    quote
}


#
# usage: mapply(My401kEtl$writeTickerQuotes, quotes, names(quotes) )
# 
# write zoo data to file
# note: write.zoo and read.zoo are NOT symmetric (c'mon guys..)
#
My401kEtl$writeTickerQuotes <- function(tickerQuotes,ticker) {
    write.zoo(tickerQuotes, file=paste("data/",trimws(ticker),".zoo",sep=""))
}


#
# usage: quotes2 <- lapply(funds$ticker, My401kEtl$readTickerQuotes)
#
# read zoo data from file
#
# @return zoo object containing ticker quotes
#
My401kEtl$readTickerQuotes <- function(ticker) {

    tickerQuotes <- tryCatch( {
            read.zoo(paste("data/",trimws(ticker),".zoo",sep=""),
                     header = TRUE, 
                     format = "%Y-%m-%d", 
                     drop=F)          
        },
        warning = function(war) {
            print(war)
            return( zoo() )
        },
        error = function(err) {
            print(err)
            return( zoo() )
        })

    tickerQuotes
}


#
# @param symbols - vector of symbol names
# 
# @return historical quotes for all symbols merged into a single zoo data frame
#
My401kEtl$downloadUniverse <- function(symbols) {
    retMe <- lapply(symbols, My401kEtl$downloadTickerQuotes)
    retMe.merged <- do.call(merge, retMe)
    names(retMe.merged) <- symbols
    retMe.merged
}


#
# Downloads quotes for all tickers in the given fund list.
#
# Each ticker's quote data is written to the file data/{ticker}.zoo.
#
# @return list of zoo objects containing ticker quotes
# 
My401kEtl$downloadAllTickerQuotes <- function(funds.df) {

    #
    # Note: this download takes a few minutes...
    #
    # quotes.list <- lapply(funds.df$ticker, My401kEtl$downloadTickerQuotes, start="1970-01-01")
    quotes.list <- lapply(funds.df$ticker, My401kEtl$downloadTickerQuotes)
    names(quotes.list) <- funds.df$ticker
    
    #
    # write zoo data to file
    # note: write.zoo and read.zoo are NOT symmetric (c'mon guys..)
    # note: assigning to dummy "x" just to avoid printing output to console
    #
    x <- mapply(My401kEtl$writeTickerQuotes, quotes.list, names(quotes.list) )

    quotes.list
}


#
# Read ticker data from cache.
#
# @return a list of zoo objects containing ticker quotes
#
My401kEtl$readAllTickerQuotes <- function(funds.df) {

    #
    # read zoo data from file
    # note: write.zoo and read.zoo are NOT symmetric (c'mon guys..)
    #
    quotes.list <- lapply(funds.df$ticker, My401kEtl$readTickerQuotes)
    names(quotes.list) <- funds.df$ticker

    quotes.list
}


#
# Loads all ticker quotes, by either:
#   a) downloading from yahoo
#   b) loading from the local cache
#   
# @return 
#
My401kEtl$loadAllTickerQuotes <- function(funds.df) {

    # Uncomment to refresh the cache.
    # My401kEtl$downloadAllTickerQuotes(funds.df)

    quotes <- My401kEtl$readAllTickerQuotes(funds.df)

    quotes <- My401kEtl$filterAndMerge(quotes)

    quotes
}


#
# Filter out funds with insufficient data and merge all tickers into a
# single zoo object.
#
# @return single zoo object containing all ticker quotes 
#
My401kEtl$filterAndMerge <- function(quotes) {

    #
    # skip the first entry in each list since it tends to be on a non-month
    # boundary and therefore screws up the merge for that month
    #
    quotes <- lapply(quotes, function(quote) { quote[-1] })
    quotes <- do.call( merge, quotes )

    #
    # correct the names
    #
    names(quotes) <- sub("AdjClose.","",names(quotes))

    #
    # remove funds with less than 5 years of data
    # 
    quote.counts <- apply( !is.na(quotes),2,sum)
    keep.symbols <- (quote.counts >= 60)
    quotes <- quotes[, names(keep.symbols)[ keep.symbols ] ]

    quotes
}




