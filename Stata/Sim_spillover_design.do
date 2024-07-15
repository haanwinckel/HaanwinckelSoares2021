
set seed 44
set sortseed 4444

*Gen cities
clear
set obs 4000
gen city = _n
gen J_city = floor(runiform()*(50-10+1) +10)
gen city_effect = 0.25*rnormal()
gen city_trend_ineq = 0.1*rnormal()
expand J_city
bysort city: gen firm = _n
gen firm_effect = 0.5*rnormal()
gen firm_mean_worker_fe = 0.5*firm_effect + 0.25*rnormal()
gen firm_N = floor(exp(firm_effect)*runiform()*(50-10+1) +10)
expand firm_N
bysort city firm: gen worker = _n
gen worker_effect = firm_mean_worker_fe + rnormal()
gen lw = city_effect + firm_effect + worker_effect + 0.1*rnormal()

//Identify low firm effect
summ firm_effect, det
gen rel_firm_effect_real = firm_effect - r(p5)
gen rel_worker_effect_real = worker_effect + r(p5)

//Measures with error
gen hat_firm_effect = rel_firm_effect_real + rnormal()/sqrt(firm_N)
gen hat_worker_effect = rel_worker_effect_real + rnormal()/sqrt(firm_N)

//Min wage
summ lw, det
gen lmw = r(p25)

//Affected workers
gen direct_effect = lw < lmw
gen affected_workers_real = rel_worker_effect_real < lmw
gen affected_workers_hat = hat_worker_effect < lmw

collapse (mean) city_* firm_N firm_effect ///
    rel_* hat_* direct_effect affected_workers_*, by(city firm)
	
bysort city: egen city_affected = mean(direct_effect)
	
//Real change
gen D_ineq = -0.5*direct_effect + 0.3 * affected_workers_real + ///
	0.1 * affected_workers_real * city_affected + city_trend_ineq + 0.1*rnormal()

*Estimate
areg D_ineq direct_effect affected_workers_hat ///
	c.affected_workers_hat#c.city_affected, abs(city)