# E-backtesting
R code for the numerical study in "E-backtesting" by Q. Wang, R. Wang, and J.F. Ziegel (2022)

* Compare_GREEGREL.R: compare GREE, GREL, and GREM methods in Example 7
* E-backtest.R: e-backtest via stationary time series in Section 7.1
* Compare_HD22: compare e-backtesting method with the monitoring method in Hoga and Demetrescu (2022) on detecting structral change in Section 7.2
* E-backtest_Real.R: financial data analysis via NASDAQ index in Section 8.1
* E-backtest_Portfolio.R: financial data analysis via optimized portfolios in Section 8.2
* Rfns.R: all relevant functions for e-backtest methods
* Rfns_HD22.R: functions used to replicate the monitoring method in Hoga and Demetrescu (2022) for comparison
* nasdaq2021.csv: data file of NASDAQ Composite index (Jan 16,1996 - Dec 31,2021)
* portfolio.csv: data file of portfolio negated percentage log returns (Jan 5, 2001 - Dec 31, 2021)

R code for additional numerical result in "Simulation and data analysis for e-backtesting" by Q. Wang, R. Wang, and J.F. Ziegel (2022)

* IID.R: e-test with iid observations in Section 2.1
* IID_traditional.R: traditional testing method in McNeil and Frey (2000) with iid observation in Section 2.1 for comparison
