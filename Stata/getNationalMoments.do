args bootstrapId

use AlmeidaCarneiroData, clear

if `bootstrapId' == 0 {
	local outputFolder "../CleanData"
}
else {
	local outputFolder "../CleanData/Bootstrap"
	local outputSuffix "_`bootstrapId'"
	local seed = 123 + `bootstrapId'
	set seed `seed'
	bsample, strata(uf_code)
}

***NOTE:
* Many lines of code below were copied from the original code from 
* Rita Almeida and Pedro Carneiro. We thank those authors for making their
* code and data available.

*Drop useless dummies:
drop duf3 duf5 duf27

*Choose base category: sp
drop duf20

keep if keep == 1


*1) Table 1 - there is a typo for total_target_workers_pc

gen inf_rate = sh_pop_status4/(sh_pop_status3 + sh_pop_status4)

ge ltotal_target_workers_pc = ln((total_target_workers+1)/emp)

ge piball_uf = lpiballpc
ge linsp_uf = linspectors_pf_ufb

local variables " linsp piball cg acesso HH "

foreach var in `variables' {
	egen m`var'_uf = mean(`var'_uf)
	local m`var'_uf = m`var'_uf
}

egen mminhor1 = mean(minhor1)
local minhor1 = mminhor1

global controls minhor1 minhor1_2 minhor1_cg minhor1_acesso minhor1_piball minhor1_HH caphor1 caphor1_2 caphor1_linsp caphor1_cg caphor1_acesso caphor1_piball caphor1_HH  ltransp_cost95_linsp ltransp_cost95 ltransp_cost95_2 ltransp_cost95_cg ltransp_cost95_acesso ltransp_cost95_piball ltransp_cost95_HH city_latitude city_longitude city_altitude larea2000 lpop91_amc lrenda_pc_91 lpop80_amc lpib80pc_amc lpop70_amc lpib70pc_amc
qui: ivregress 2sls inf_rate (ltotal_target_workers_pc = minhor1_linsp) $controls duf*
local dInf_dlogEnf = _b[ltotal_target_workers_pc]
di "#`bootstrapId': beta = `dInf_dlogEnf'"


clear
set obs 1
local momNames laborShare share_500_100
local momValues 0.528 0.72235
local momSE 0.0028 0.02
forvalues i = 1/2 {
	local varname : word `i' of `momNames'
	local val : word `i' of `momValues'
	local se : word `i' of `momSE'
	if `bootstrapId' == 0 {
		gen z_`varname' = 0
	}
	else {
		gen z_`varname' = rnormal()
	}
	gen `varname' = `val' + z_`varname' * `se'
}
drop z_*
gen infElast = `dInf_dlogEnf'
outsheet using "`outputFolder'/national_2003`outputSuffix'.csv", ///
	replace comma nolabel
