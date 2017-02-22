
# Portfolio Analysis of My 401K


## Financial Cluster Analysis in academia / prior art

* Categorizing Mutual Funds using Clusters
    * marathe and shawky
    * [http://www.academia.edu/23287051/CATEGORIZING\_MUTUAL\_FUNDS\_USING\_CLUSTERS\_by](http://www.academia.edu/23287051/CATEGORIZING_MUTUAL_FUNDS_USING_CLUSTERS_by)
    * 1995 data
    * uses 28 variables
        * 1 year return
        * 3 year return
        * 5 year return
        * turnover
        * % bonds, % stocks, % cash
        * etc
    * uses hartigan rule-of-thumb to determine number of clusters
    * checks robustness of clustering using PCA
* Clustering Financial Time Series
    * Pattarin, Minerva, Patterlini
    * [https://www.researchgate.net/publication/222011409\_Clustering\_financial\_time\_series\_An\_application\_to\_mutual\_funds\_style\_analysis](https://www.researchgate.net/publication/222011409_Clustering_financial_time_series_An_application_to_mutual_funds_style_analysis)
    * uses monthly return data
    * uses PCA to reduce data
    * uses evolutionary clustering algorithm for clustering
    * shows results for 2 thru 10 clusters
* Patent: [https://www.google.com/patents/US20030065602](https://www.google.com/patents/US20030065602)
    * number of clusters is determined by 
        * the number of largest principal components
        * that explain most of the variance
        * in the sample covariance matrix of returns
    * uses covariance matrix to determine clusters
        * uh.. no
        * why? cuz funds with low correlations with other funds
            * can be un-correlated in different ways
            * and un-correlated with each other
            * but still clustered together
* [http://xiforofinanzas.ua.es/trabajos/1077.pdf](http://xiforofinanzas.ua.es/trabajos/1077.pdf)
    * classifies based on fund characterstics, not time-series returns
    * interesting discussion of self-organizing maps


## What does "done" look like?

* all funds clustered into a dozen or so groups, pick a “representative” fund from all groups
* compute rolling “growth-of-$1” returns
    * don’t forget to factor in fees!
* compute correlations between “representative” funds
    * compute rolling correlations
    * find funds that are particularly UN-correlated with representative funds
* compute optimum mean-variance portfolios based on rolling correlations
    * change over time!
* write program that does this every once in a while
    * re-balance portfolio


## TODO

* Is there a machine-learning way to determine the optimal portfolio allocations?
    * Need a cost function
    * Need a response variable?
        * Monthly return: target `>> 0`
        * which can be simulated with logistic regression
        * essentially all target y's are 1
        * x is the monthly return vector across selected funds (or all funds?)
        * BUT!! additional constraint: weights (theta) must sum to 1
        * additional constraint: weights (theta) must be positive!
            * could use regularization to try to enforce 
            * heavy penalty for negative theta
            * heavy penalty for sum(theta) larger than 1


TODO: Research all funds
      Full return history
      Recent return history
      Recession return history


1. get all fund ticker symbols from netbenefits.com
2. download historical fund price data from yahoo
3. Preliminary analysis
    * determine correlations between funds
    * group funds in broad categories based on correlations
        * large cap
        * small cap
        * commodities
            * metals
        * bonds
           * govt
           * investment grade corp
           * junk corp
        * real estate
    * analyze "normality" of monthly and daily returns

4. do basic portfolio analysis using CER model and sample stats
    * subtract fees from fund performance
    * note: some funds overlap.
        * this breaks IID assumption
    * adjust for divs? inflation?
5. determine sharpe ratios for all funds
6. determine MCR for all funds



## My Investments

* LARGE COMPANY IDX
    * Tracks the S&P500 (^GSPC)
    * 28.70%
* LARGE-CAP GROWTH IDX
    * Tracks Russell 1000 Growth Index (symbol: ^RLG, VRGWX)
    * 28.55%
* VANG ST TREASURY ADM (VFIRX)
    * Short-term Govt Bonds
    * 23.97%
* MODERATE LS
    * Weight: 18.78%
    * 12% Inflation Protected Bond Fund
        * Tracks Barclays U.S. Treasury Inflation Protected Securities - Series L Index.
    * 20% Total Bond Market & Interest Income Funds
        * Total Bond Market: Tracks Barclays U.S. Aggregate Bond Index
        * Interest Income Fund: ?
    * 3% High Yield & Emerging Markets Bond Fund
        * ?
    * 6% Real Estate Investment Trust Index Fund
    * 2% International Real Estate Index Fund
    * 10% Balanced Exposure Fund
    * 3% Commodities Fund
    * 31% Total Stock Market Index Fund 
    * 13% Total International Stock Market Index Fund. 
    * State Street Global Advisors periodically rebalances the fund to the target allocations.


## Russell Indices:


* Russell 3000 Index 
    * 1-3000
    * Top 3000 U.S. companies (98% of investable market)
    * ^RUA, ^RAG, ^RAV  (all, growth, value)
* Russell 2500 Index
    * 501-3000
    * bottom 2,500 stocks in the Russell 3000 Index.
* Russell 2000 Index 
    * 1001-2000
    * The small-cap benchmark index of the bottom 2,000 stocks in the Russell 3000 Index.
    * ^RUT, ^RUO, ^RUJ 
* Russell 1000 Index 
    * 1-1000
    * The large-cap index of the top 1,000 stocks in the Russell 3000 Index.
    * ^RUI, ^RLG, ^RLV 
* Russell Top 200 Index 
    * 1-200
    * The mega-cap index of the very largest 200 stocks in the Russell 3000 Index.
    * ^RT200, ^RT200G, ^RT200V: 
    * IWL, IWY, IWX
* Russell Top 50 Index 
    * 1-50
    * Measures the performance of the 50 largest companies in the Russell 3000 Index.
    * ^RU50 
* Russell Midcap Index 
    * 201-1000
    * The bottom 800 stocks in the Russell 1000 Index. The Russell Top 200 Index plus the Russell Midcap Index yields the Russell 1000 Index.
    * ^RMCC, ^RMCCG, ^RMCCV
    * IWR, IWP, IWS
* Russell Microcap Index
    * 2001-4000
    * A micro-cap index of the stocks ranked from 2,001-4,000 in the Russell indexing universe, consisting of capitalizations ranging from about $50 million to $2.5 billion. Hence, this is an index of the 1,000 smallest Russell 3000 stocks, plus the 1,000 smaller stocks.
    * ^RUMIC, ^RUMICG, ^RUMICV
    * IWC
* Russell Small Cap Completeness Index 
    * ~501-3000
    * Top 3000 minus SP500


## Other Funds

* IPE: Barclays U.S. Treasury Inflation Protected Securities - Series L Index.
    * SPDR Barclays TIPS ETF
* AGG: Barclays U.S. Aggregate Bond Index, a widely recognized measure of the investment-grade taxable U.S. bond market. 
    * Shares Core US Aggregate Bond (AGG)
    * The index consists of more than 5,000 U.S. Treasury, federal agency, mortgage-backed, and corporate securities.
* Barclays Corp High-Yield
* Dow Jones U.S. Total Stock Market Index
    * which covers all regularly traded U.S. stocks.
    * ^DJI
* ACWX: MSCI® All Country World ex USA Investable Market Index, 
    * iShares MSCI ACWI ex US (ACWX)
    * which includes equities in Europe, Asia/Pacific, Canada, Latin America, Africa, and the Middle East.
* Global Real Estate Stock Index Fund
    * MSCIUSREIT INDEX (^RMZ) 
    * ^REIT: Dow Jones Equity All REIT Tota
    * ENDEM.L
    * composite benchmark of 70% MSCI US REIT Index and 30% FTSE EPRA/NAREIT Developed ex-US Rental Index.
* Barclays U.S. Long Credit Index
* IEUR: MSCI® Europe Index
    * iShares Core MSCI Europe (IEUR) 
    * which is made up of stocks from 16 European countries,
* IPAC: MSCI® Pacific Index, 
    * iShares Core MSCI Pacific (IPAC)
    * which is made up of stocks from companies in Australia, Hong Kong, Japan, New Zealand, and Singapore. 
    * Japanese stocks represent a significant majority of the index and, thus, of the fund's assets. 
* FTSE Emerging Markets All Cap China A Transition Index, 
    * which is made up of stocks from emerging market countries in Europe, Asia, Africa, and Latin America.
* VWO: Vanguard FTSE Emerging Markets ETF (VWO)
* ^RMZ: MSCI US REIT Index 
    * MSCIUSREIT INDEX (^RMZ) 
    * a portfolio of equities whose total rate of return will approximate the capitalization weighted total rate of return of about 90% of the market.
    * the investment universe consists of the United States market for publicly traded real estate equity properties.
* FTSE EPRA/NAREIT Developed ex-U.S. Rental Index 
    * ENDEM.L
    * a portfolio of equities whose total rate of return will approximate the capitalization weighted total rate of return of about 90% of the market.
    * The investment universe consists of the international market for securities of companies principally engaged in the real estate industry and other real estate related investments. 
* FNMIX: Fidelity New Markets Income Fund
    * Normally investing at least 80% of assets in securities of issuers in emerging markets and other investments that are tied economically to emerging markets. 
    * Normally investing primarily in debt securities of issues in emerging markets
* PELBX: PIMCO Emerging Local Bond Fund Institutional Class
    * the fund invests at least 80% of its assets in Fixed Income Instruments denominated in currencies of countries with emerging securities markets, 
    * which may be represented by forwards or derivatives such as options, futures contracts or swap agreements. 
    * It may invest without limitation in Fixed Income Instruments that are economically tied to emerging market countries. 
    * It is non-diversified.
* PEBIX: PIMCO Emerging Markets Bond Fund Institutional Class
    * The fund normally invests at least 80% of its assets in Fixed Income Instruments that are economically tied to emerging market countries, 
    * which may be represented by forwards or derivatives such as options, futures contracts or swap agreements. 
    * It may invest in both investment-grade securities 
    * and junk bonds subject to a maximum of 15% of its total assets in securities rated below B by Moody's, 
        * or equivalently rated by S&P or Fitch, or, if unrated, determined by PIMCO to be of comparable quality.
* FAGIX: Fidelity Capital & Income Fund
    * Investing in equity and debt securities, including defaulted securities, with an emphasis on lower-quality debt securities. 
    * Investing in companies in troubled or uncertain financial condition.
* FFRHX: Fidelity Floating Rate High Income Fund
    * Normally investing at least 80% of assets in floating rate loans, which are often lower-quality debt securities, and other floating rate debt securities. 
    * Investing in companies in troubled or uncertain financial condition. 
    * Investing in money market and investment grade debt securities, and repurchase agreements.
* SPHIX: Fidelity® High Income Fund
    * Normally investing primarily in income-producing debt securities, preferred stocks, and convertible securities, with an emphasis on lower-quality debt securities. 
    * Investing in companies in troubled or uncertain financial condition. 
    * Potentially investing in non-income producing securities, including defaulted securities and common stocks. 
"PHIYX","PIMCO High Yield Fund Institutional Class"

* FID GNMA (FGMNX)
* FID INFLAT PROT BOND (FINPX)
* FID INTM GOVT INCOME (FSTGX)
* PIM GNMA INST (PDMIX)
* PIM MORTGAGE BCKD IS (PTRIX)
* VANG GNMA ADM (VFIJX)
* VANG INFL PROT INST (VIPIX)
* VANG INTM TREAS ADM (VFIUX)
* DODGE & COX INCOME (DODIX)
* FID INTERMED BOND (FTHRX)
* FID TOTAL BOND (FTBFX)
* FIDELITY GOVT INCOME (FGOVX)
* PIM INVT GRD BD INST (PIGIX)
* PIM MOD DURAT INST (PMDRX)
* PIM TOT RT III INST (PTSAX)
* PIM TOTAL RT INST (PTTRX)
* VANG INTM BOND INST (VBIMX)
* VANG INTM INV GR ADM (VFIDX)
* PIM LT US GOVT INST (PGOVX)
* PIM REAL RETURN INST (PRRIX)
* PIM RL RT ASSET INST (PRAIX)
* VANG LT TREASURY ADM (VUSUX)
* VANG LT BOND IDX IS (VBLLX)
* VANG LT INV GR ADM (VWETX)
* FID STRATEGIC INCOME (FSICX)
* PIM DIVERS INC INST (PDIIX)
* PIM UNCONSTRNED BD I (PFIUX)
* FID LTD TERM GOVT (FFXSX)
* VANG ST FEDERAL ADM (VSGDX)
* VANG ST TREASURY ADM (VFIRX)
* FID SHORT TERM BOND (FSHBX)
* PIM LOW DUR III INST (PLDIX)
* PIM LOW DUR INST (PTLDX)
* VANG ST BOND IDX IS (VBITX)
* VANG ST INVT GR INST (VFSIX)
* AF CAP WORLD BOND R6 (RCWGX)
* DODGE & COX GLBL BND (DODLX)
* PIM FOR BD UNHG INST (PFUIX)
* PIM FOR BD US$HG I (PFORX)
* PIM GLOB BD UNHG I (PIGLX)
* PIM GLOB BD US$HG I (PGBIX)
* VANG WELLESLEY ADM (VWIAX)
* FID CONVERTIBLE SEC (FCVSX)
* VANG CONVERTIBLE SEC (VCVSX)
* AF BALANCED R6 (RLBGX)
* DODGE & COX BALANCED (DODBX)
* FID BALANCED K (FBAKX)
* PIM ALL A ALL AUTH I (PAUIX)
    * blend
* PIM ALL ASSET INST (PAAIX)
* VANG STAR (VGSTX)
* VANG WELLINGTON ADM (VWENX)
* FID GLOBAL BALANCED (FGBLX)
* FID FREEDOM K 2005 (FFKVX)
* FID FREEDOM K 2010 (FFKCX)
* FID FREEDOM K 2015 (FKVFX)
* FID FREEDOM K 2020 (FFKDX)
* FID FREEDOM K 2025 (FKTWX)
* FID FREEDOM K 2030 (FFKEX)
* FID FREEDOM K 2035 (FKTHX)
* FID FREEDOM K 2040 (FFKFX)
* FID FREEDOM K 2045 (FFKGX)
* FID FREEDOM K 2050 (FFKHX)
* FID FREEDOM K 2055 (FDENX)
* FID FREEDOM K INCOME (FFKAX)
* VANG INST TR 2010 (VIRTX)
* VANG INST TR 2015 (VITVX)
* VANG INST TR 2020 (VITWX)
* VANG INST TR 2025 (VRIVX)
* VANG INST TR 2030 (VTTWX)
* VANG INST TR 2035 (VITFX)
* VANG INST TR 2040 (VIRSX)
* VANG INST TR 2045 (VITLX)
* VANG INST TR 2050 (VTRLX)
* VANG INST TR 2055 (VIVLX)
* VANG INST TR INCOME (VITRX)
* AF FUNDMNTL INV R6 (RFNGX)
* AF INV CO OF AMER R6 (RICGX)
* DFA US CORE EQ 1 I (DFEOX)
* FID DIVIDEND GR K (FDGKX)
* FID FUND K (FFDKX)
* PIMCO STKPLUS INST (PSTKX)
* VANG DIV GROWTH INV (VDIGX)
* VANG FTSE SOC IDX IS (VFTNX)
* VANG GRTH & INC ADM (VGIAX)
* AF AMCAP R6 (RAFGX)
* AF GRTH FUND AMER R6 (RGAGX)
* AF NEW ECONOMY R6 (RNGGX)
* FID BLUE CHIP GR K (FBGKX)
* FID CAP APPREC K (FCAKX)
* FID CONTRAFUND K (FCNKX)
* FID EXPORT & MULTI K (FEXKX)
* FID FOCUSED STOCK (FTQGX)
* FID GROWTH CO K (FGCKX)
* FID GROWTH DISC K (FGDKX)
* FID INDEPENDENCE K (FDFKX)
* FID LARGE CAP STOCK (FLCSX)
* FID OTC K (FOCKX)
    * Fidelity OTC Portfolio is an open-end fund incorporated in the USA. The
      Fund's objective is capital appreciation. The Fund normally invests at
      least 80% of assets in securities principally traded on NASDAQ or another
      over-the-counter market, which has more small and medium-sized companies
      than other markets. The Fund invests more than 25% of total assets in the
      technology sector.
* FID TREND (FTRNX)
* VANG MORGAN GRTH ADM (VMRAX)
* VANG PRIMECAP CORE (VPCCX)
* AF AMER MUTUAL R6 (RMFGX)
* AF WASH MUTL INV R6 (RWMGX)
* DODGE & COX STOCK (DODGX)
* FID EQUITY INCOME K (FEIKX)
* VANG EQUITY INC ADM (VEIRX)
* VANG WINDSOR II ADM (VWNAX)
* VANGUARD WINDSOR ADM (VWNEX)
* FID LEVERGD CO STK K (FLCKX)
* FID LOW PRICED STK K (FLPKX)
* FID VALUE STRAT K (FVSKX)
* VANG MIDCAP IDX INST (VMCIX)
* VANG STRATEGIC EQ (VSEQX)
* FID MID CAP STOCK K (FKMCX)
* VANG MIDCAP GRTH INV (VMGRX)
* DFA US VECTOR EQ I (DFVEX)
* FID VALUE K (FVLKX)
* VANG SELECTED VALUE (VASVX)
* DFA US SMALL CAP I (DFSTX)
* FID SMALL CAP STOCK (FSLCX)
* VANG SM CAP IDX INST (VSCIX)
* AF SMALLCAP WORLD R6 (RLLGX)
* FID SM CAP DISCOVERY (FSCRX)
* FID STK SEL SM CAP (FDSCX)
* VANG EXPLORER ADM (VEXRX)
* DFA US TARGET VALUE (DFFVX)
* PIM COM REAL RET I (PCRIX)
* FID REAL ESTATE INC (FRIFX)
* FID REAL ESTATE INVS (FRESX)
* PIM RE REAL RET INST (PRRSX)
* AF NEW WORLD R6 (RNWGX)
* DFA EMERG MKTS VALUE (DFEVX)
* DFA EMERGING MARKETS (DFEMX)
* FID PACIFIC BASIN (FPBFX)
* FID EUROPE (FIEUX)
* AF EUROPAC GROWTH R6 (RERGX)
* DFA LARGE CAP INTL I (DFALX)
* FID INTL DISCOVERY K (FIDKX)
* FID OVERSEAS K (FOSKX)
* VANG INTL GROWTH ADM (VWILX)
* VANGUARD INTL VALUE (VTRIX)
* FID DIVERSIFD INTL K (FDIKX)
* DFA INTL VALUE I (DFIVX)
* DODGE & COX INTL STK (DODFX)
* FID INTL SMALL CAP (FISMX)
* VANG INTL EXPLOR INV (VINEX)
* DFA INTL SMALL CO I (DFISX)
* DFA INTL VECTOR EQ I (DFVQX)
* DFA GLOB REAL ESTATE (DFGEX)
* FID INTL REAL ESTATE (FIREX)
* FID JAPAN (FJPNX)
* FID LATIN AMERICA (FLATX)
* FID CANADA (FICDX)
* FID CHINA REGION (FHKCX)
* FID EMERGING ASIA (FSEAX)
* AF CAP WORLD G&I R6 (RWIGX)
* AF INTL GTH & INC R6 (RIGGX)
* AF NEW PERSPECT R6 (RNPGX)
* DODGE & COX GLOBAL (DODWX)
* VANG GLB MIN VOL ADM (VMNVX)
* VANG GLOBAL EQ INV (VHGEX)
"https://en.wikipedia.org/wiki/IShares","TODO: iShares ETFs: [https://en.wikipedia.org/wiki/IShares]"
    * analyze all


## Simple vs. CC (log) Returns

* Distributions:
    * monthly cc returns tend to approximate normal distribution
    * better than daily cc returns
    * normal distribution may not be suitable for simple returns
        * since normal distribution allows for impossible simple returns
            * e.g. `R < -1`
* cc return ~= simple return, for small values of R
* multi-period cc returns are ADDITIVE
    * linear combinations of normally distributed random variables
        * are ALSO normally distributed random variables
    * aggregation of multi-period cc returns is the same as...
        * linear combination of normally-distributed random variables
* multi-period simple returns are NOT additive
    * however they are approximately additive for small values of R
* CER model assumes IID normally distributed CC returns
* Portfolio theory typically uses simple returns
    * since portfolio return is linear combo of simple returns


