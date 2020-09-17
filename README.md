# keplerspline
Wrapper to fit a basis spline to 1d data with outlier rejection. The spline code itself is from the SDSS IDL utils library (http://www.sdss3.org/dr9/software/idlutils_doc.php). 

Usage: 

t = dindgen(500)/50         

f = jrandomn(seed, 500) + sin(t)

plot, t,f,/yno, ps = 3

s = keplerspline(t, f, ndays = 1.0); Fit a spline with a specified break-point spacing of 1.0 day. 

oplot, t, s

include = where(t lt 10 or t gt 12);

s2 = keplersplinetest(t, f, ndays = 1.0, include = include); Fit a spline with a break-point spacing of 1.0 day. Exclude data taken between t = 10 and t = 12. 

oplot, t, s2

s3 = choosekeplerspline(t, f); Fit a spline, while choosing the best break-point spacing to use based on the Bayesian Information Criterion. 

oplot, t, s3

s4 = choosekeplerspline(t, f); Fit a spline, while choosing the best break-point spacing to use based on the Bayesian Information Criterion and 
;excluding data taken between t = 10 and t = 12.

oplot, t, s4


