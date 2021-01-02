* Sample is only the first 50 observations of chapter 10 data
clear all
set seed 10101

cd /Users/lininghui/Desktop/Microeconometrics/Code_Data
use mus10data.dta
quietly keep if year02 == 1
quietly drop if _n > 50
quietly keep docvis chronic age
quietly save bootdata.dta, replace

* Option vce (bootstrap) to compute bootstrap standard errors
poisson docvis chronic
poisson docvis chronic , vce(boot, reps(400) seed(10101) nodots)
/*
nodots and dots(#) specify whether to display replication dots. By default, one dot character is displayed for each successful replication. A red ‘x’ is displayed if command returns an error or if
any value in exp list is missing. Y
*/

* How many bootstraps?
* Bootstrap standard errors for different reps and seeds
quietly poisson docvis chronic, vce(boot, reps(50) seed(10101))
estimates store boot50
quietly poisson docvis chronic, vce(boot, reps(50) seed(20202))
estimates store boot50diff
quietly poisson docvis chronic, vce(boot, reps(2000) seed(10101))
estimates store boot2000
quietly poisson docvis chronic, vce(robust)
estimates store robust
estimates table boot50 boot50diff boot2000 robust, b(%8.5f) se(%8.5f)

* clustered bootstraps
* Option vce(boot, cluster) to compute cluster-bootstrap standard errors
poisson docvis chronic, vce(boot, cluster(age) reps(400) seed(10101) nodots)

poisson docvis chronic, vce(cluster age)

* Bootstrap confidence intervals example
* Bootstrap confidence intervals : normal-based, percentile, BC , and BCa
quietly poisson docvis chronic, vce(boot, reps(999) seed(10101) bca)
estat bootstrap, all

matrix list e(b_bs)

* Bootstrap command applied to Stata estimation command
bootstrap, reps(400) seed(10101) nodots noheader: poisson docvis chronic

*Bootstrap standard error estimate of the standard error of a coeff estimate
bootstrap _b _se, reps(400) seed(10101) nodots: poisson docvis chronic
bootstrap _b _se, reps(400) seed(10101) nodots: poisson docvis chronic, vce(robust)

* Program to return b and robust estimate V of the VCE
program poissrobust, eclass
version 10.1
tempname b V
poisson docvis chronic , vce(robust)
matrix `b' = e(b)
matrix `V' = e(V)
ereturn post `b' `V'
end

* Check preceding program by running once
poissrobust
ereturn display

* Bootstrap standard-error estimate of robust standard errors
bootstrap _b _se , reps(400) seed(10101) nodots nowarn: poissrobust

* Set up the selection model two-step estimator data of chapter 16
use mus16data.dta, clear
generate y = ambexp
generate dy = y > 0
generate lny = ln(y)
global xlist age female educ blhisp totchr ins

* Program to return b for Heckman 2-step estimator of selection model
program hecktwostep, eclass
version 10.1
tempname b V
tempvar xb
capture drop invmills
probit dy $xlist
predict `xb', xb
generate invmills = normalden(`xb')/normal(`xb')
regress lny $xlist invmills
matrix `b' = e(b)
ereturn post `b'
end

hecktwostep
ereturn display

*Bootstrap for Heckman two-step estimator using chapter 16 example
bootstrap _b, reps(400) seed(10101) nodots nowarn: hecktwostep

* Program to return (b1-b2) for Hausman test of endogeneity
program hausmantest, eclass
version 10.1
tempname b bols biv
regress ldrugexp hi_empunion totchr age female blhisp linc, vce(robust)
matrix `bols' = e(b)
ivregress 2sls ldrugexp (hi_empunion = ssiratio) totchr age female blhisp linc , vce(robust)
matrix `biv' = e(b)
matrix `b' = `bols' - `biv'
ereturn post `b'
end

* Bootstrap estimates for Hausman test using chapter 6 example
use mus06data.dta , clear
bootstrap _b, reps(400) seed(10101) nodots nowarn: hausmantest

* Perform Hausman test on the potentially endogenous regressor
test hi_empunion

use bootdata.dta , clear
summarize docvis

bootstrap coeffvar=(r(sd)/r(mean)), reps(400) seed(10101) nodots nowarn: summarize docvis

* Percentile-t for a single coefficient: Bootstrap the t statistic
use bootdata.dta, clear
quietly poisson docvis chronic, vce(robust)
local theta = _b[chronic]
local setheta = _se[chronic]
bootstrap tstar = ((_b[chronic]-`theta')/_se[chronic]), seed(10101) ///
rep(999) nodots saving(percentilet, replace): poisson docvis chronic, ///
vce(robust)

* Percentile-t p-valuc for symmetric two-sided Wald test of H O : theta = 0
use percentilet, clear
quietly count if abs(`theta'/`setheta') < abs(tstar)
display "p-value = " r(N)/_N

* Percentile-t critical values and confidence interval
_pctile tstar, p(2.5, 97.5)
scalar lb = `theta' + r(r1)*`setheta'
scalar ub = `theta' + r(r2)*`setheta'
display "2.5 and 97.5 percentiles of t* distn: " r(r1) " , " r(r2) _n "95 percent percentile-t Confidence interval is (" lb "," ub ")"




