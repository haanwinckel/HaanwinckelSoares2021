args year bootstrapId

local mdFolder "./Panel"
if `year' == 2003 {
	use "`mdFolder'/pmenova_painel_C_rs.dta", clear
	foreach panel in D E F G H {
		append using "`mdFolder'/pmenova_painel_`panel'_rs.dta"
	}
}
else if `year' == 2012 {
	use "`mdFolder'/pmenova_painel_Q_rs.dta", clear
	foreach panel in R S T {
		append using "`mdFolder'/pmenova_painel_`panel'_rs.dta"
	}
}

rename v070 month
rename v075 year
rename v035 region

gen t = (year-2003)*12+month

*Keeping data from April through December
keep if (year == `year' & month >= 4)

rename vD3 worforce
rename vD23 usual_income_main_job
rename vD27 hours
rename vD17 sector
rename vD1 empl_status
rename vD19 firm_size_cat
rename v215 samp_weight
rename v407A occupation
rename v408A industry
rename v406 n_jobs

if `bootstrapId' > 0 { //Bootstrap reps
	local seed = 123 + `bootstrapId'
	set seed `seed'
	bsample, strata(region) cluster(idind) idcluster(newidind)
	drop idind
	rename newidind idind
}

*Keep people in the labor force
keep if worforce==1


*Generate minimum wage value, as hourly wage
gen mw = 240/(44*(31/7))
//Inflation adjustnment
if `year' == 2012 {
	replace mw = mw * 1.6073
}

//Sample selection
replace hours = . if hours < 5 | hours > 70
_pctile vD25, percentiles(0.5 99.5)
replace vD25 = . if vD25 < r(r1)
replace vD25 = . if vD25 > r(r2) 
gen wh=vD25/(hours*(30.4/7))
replace wh = wh/mw //normalize so 1 == mw
gen lwh = log(wh)

*Employment status
gen wage_worker = sector<=2
gen informal = sector==2 if wage_worker==1
gen unemp=empl_status==2 if empl_status==1 | empl_status==2

gen binding=wh<=1.2 if informal == 0 & unemp == 0

*Firm size
gen fs_6_10 = v412 == 2 if v412 != .
gen fs_11p = v412 == 3 if v412 != .
*Some interactions
gen fs_1_5_X_for = (1-informal)*(1-fs_6_10-fs_11p)
gen fs_6_10_X_inf = fs_6_10*informal
gen fs_6_10_X_for = fs_6_10*(1-informal)
gen fs_11p_X_inf = fs_11p*informal
gen fs_11p_X_for = fs_11p*(1-informal)

*Educational groups
gen college = v303==9 | v307==9 | ///
	(v307==6 & v311==1)
gen college_drop = (v307==6 | v303==5) ///
	& college==0
rename vDAE2 ed_group
replace ed_group=3 if college_drop==1
replace ed_group=3 if college==1

*Age groups
gen birth_year = v224 if v224 != 9999
gen birth_month = v214 if v214 != 99
replace birth_month = 6.5 if v214 == 99
gen age = ((year-birth_year)*12 + (month - birth_month))/12
drop if age < 16 | age > 60
gen age_group = 1
replace age_group = 2 if age >= 20
replace age_group = 3 if age >= 25
replace age_group = 4 if age >= 30


***INTERMISSION***
*Getting summary data for Table 1*
if `bootstrapId' == 0 {
preserve
replace ed_group = 4 if college == 1
gen lwh_formal = lwh if informal == 0
gen lwh_informal = lwh if informal==1
gen fs_1_5_inf = 1-fs_6_10-fs_11p if informal==1
gen fs_11p_inf = fs_11p if informal==1
gen fs_1_5_for = 1-fs_6_10-fs_11p if informal==0
gen fs_11p_for = fs_11p if informal==0
gen ones = 1
collapse (mean) informal lwh_* unemp ///
	fs_1_5_inf fs_11p_inf fs_1_5_for fs_11p_for (sum) N = ones ///
	[aw = samp_weight], by(ed_group)
egen sumN = total(N)
gen share = N/sumN
drop N sumN
gen year = `year'
tempfile tmp
save `tmp'
collapse (mean) informal lwh_* unemp year ///
	fs_1_5_inf fs_11p_inf fs_1_5_for fs_11p_for [aw = share]
gen share = 1
gen ed_group = 0
append using `tmp'
outsheet using "../Output/summaryTable1_`year'.csv", replace
restore
}
***END INTERMISSION***


keep if unemp == 1 | informal != .

egen id = group(idind)
bys id t: drop if _n>1
xtset id t

by id: replace age_group = age_group[1]
by id: replace ed_group = ed_group[1]
gen wk_group = (ed_group - 1)*4 + age_group

levelsof wk_group, local(wk_levels)
levelsof region, local(region_list)

*First: run regressions to get formal wage premia and 
* transition probabilities
foreach r of local region_list {

	preserve
	keep if region == `r'

	foreach wk_g of local wk_levels {
		
		*Overall formal wage premium
		areg lwh informal age i.t [aw = samp_weight] ///
			if n_jobs==1 & wk_group == `wk_g', ///
			abs(idind) vce(cluster idind)
		local wp_formal_`wk_g' = exp(-_b[informal])-1
		
		*Overall firm size premium
		areg lwh fs_6_10 fs_11p age i.t [aw = samp_weight] ///
			if n_jobs==1 & wk_group == `wk_g', ///
			abs(idind) vce(cluster idind)
		local wp_fs_6_10_`wk_g' = exp(_b[fs_6_10])-1
		local wp_fs_11p_`wk_g' = exp(_b[fs_11p])-1
		
		*Firm size premium conditional on sectors 
		areg lwh fs_*_X_for fs_*_X_inf age i.t [aw = samp_weight] ///
			if n_jobs==1 & wk_group == `wk_g', ///
			abs(idind) vce(cluster idind)
		local wp_fs_for_6_10_`wk_g' = exp(_b[fs_6_10_X_for]-_b[fs_1_5_X_for])-1
		local wp_fs_for_11p_`wk_g' = exp(_b[fs_11p_X_for]-_b[fs_1_5_X_for])-1
		local wp_fs_inf_6_10_`wk_g' = exp(_b[fs_6_10_X_inf])-1
		local wp_fs_inf_11p_`wk_g' = exp(_b[fs_11p_X_inf])-1
	}


	*Transitions out of unemployment 
	gen trans_U2E = (1-unemp) if L.unemp == 1
	foreach wk_g of local wk_levels {
		summ trans_U2E if wk_group == `wk_g' [aw=samp_weight]
		local trans_U2E_`wk_g' = r(mean)
	}

	*Get other moments using collapse
	gen lwh_informal = lwh if informal==1
	gen fs_6_10_inf = fs_6_10 if informal==1
	gen fs_11p_inf = fs_11p if informal==1
	gen fs_6_10_for = fs_6_10 if informal==0
	gen fs_11p_for = fs_11p if informal==0
	gen ones = 1
	collapse (mean) ed_group age_group binding informal lwh fs_6_10 fs_11p ///
		fs_6_10_inf fs_11p_inf fs_6_10_for fs_11p_for (sum) N = ones ///
		(p1) lwh_inf_p01 = lwh_informal (p5) lwh_inf_p05 = lwh_informal ///
		(p10) lwh_inf_p10 = lwh_informal (p20) lwh_inf_p20 = lwh_informal ///
		[aw = samp_weight], by(wk_group)
	gen w = exp(lwh)

	*Recover others
	gen wp_formal = .
	gen wp_fs_6_10 = .
	gen wp_fs_11p = .
	gen wp_fs_for_6_10 = .
	gen wp_fs_for_11p = .
	gen wp_fs_inf_6_10 = .
	gen wp_fs_inf_11p = .
	gen trans_U2E = .
	foreach wk_g of local wk_levels {
		replace wp_formal = `wp_formal_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_6_10 = `wp_fs_6_10_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_11p = `wp_fs_11p_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_for_6_10 = `wp_fs_for_6_10_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_for_11p = `wp_fs_for_11p_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_inf_6_10 = `wp_fs_inf_6_10_`wk_g'' ///
			if wk_group == `wk_g'
		replace wp_fs_inf_11p = `wp_fs_inf_11p_`wk_g'' ///
			if wk_group == `wk_g'
		replace trans_U2E = `trans_U2E_`wk_g'' ///
			if wk_group == `wk_g'
	}

	egen sample_size_region = total(N)
	gen group_share = N/sample_size_region
	drop N

	gen region = `r'	
	tempfile temp_`r'
	save `temp_`r''
	
	restore
}

*Add transitions out of employment


gen trans_2U = unemp if L.informal != . 
summ trans_2U
local trans_2U_overall = r(mean)

gen trans_2U_cox = 1 if F1.unemp == 1 & informal != . //Transitions into unemployment
by id: replace trans_2U_cox = 0 if _n==_N & informal != . //Right censoring

*Estimate monthly rate of transitions into unemployment
gen tenure = 1 if v427 == 1
replace tenure = v4272 + 1 if v427 == 2
replace tenure = v4275 + 13 if v427 == 3
replace tenure = 12*v4274 + 6 if v427 == 4

by id: egen mean_lwh = mean(lwh)
bys region wk_group: egen mean_lwh_region = mean(lwh)  //Control -- proxy for unobserved worker ability within groups
gen rel_mean_lwh = mean_lwh-mean_lwh_region
gen rel_mean_lwh2 = rel_mean_lwh^2
gen rel_mean_lwh3 = rel_mean_lwh^3
by region wk_group: egen mean_age_region = mean(age)
gen rel_age = age - mean_age_region

stset tenure trans_2U_cox
stcox informal i.wk_group i.region rel_age rel_mean_lwh* i.t

clear
foreach r of local region_list {
	append using `temp_`r''
}
gen N = sample_size_region * group_share
egen tot_N = total(N)
gen group_share_all = N / tot_N
drop N tot_N sample_size_region

*Put info relative to job destruction
gen rel_trans_F2U = .
local hr_inf = exp(_b[informal]) //Assumed to be the same for all types
foreach r of local region_list {
	local hr_region = exp(_b[`r'.region])
	foreach wk_g of local wk_levels {
		local hr_educ = exp(_b[`wk_g'.wk_group])
		replace rel_trans_F2U = `hr_region' * `hr_educ' ///
			if region == `r' & wk_group == `wk_g'
	}
}
gen rel_trans_I2U = rel_trans_F2U * `hr_inf'
gen total_rel_trans = ((1-informal)*rel_trans_F2U + informal*rel_trans_I2U)*group_share_all
egen total_rel_E2U = total(total_rel_trans)
gen trans_F2U = rel_trans_F2U * `trans_2U_overall' / total_rel_E2U
gen trans_I2U = rel_trans_I2U * `trans_2U_overall' / total_rel_E2U
drop rel_trans* total_rel_*


if `bootstrapId' > 0 {
	local outputFileStub "../CleanData/Bootstrap/pme_`year'_"
	local suffix "_`bootstrapId'"
}
else {
	local outputFileStub "../CleanData/pme_`year'_"
	local suffix ""
}
	
outsheet using "`outputFileStub'regions`suffix'.csv", replace comma nolabel
