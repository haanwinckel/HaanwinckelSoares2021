

clear all
set more off
set matsize 4000
set maxvar 10000
set max_memory 12g

use "census_collapsed_data.dta", clear

gen wei = .
gen formal = .
gen employed = .
sort micror year
foreach v in lincome educ_8p formal_noeduc employed_noeduc {
	by micror: gen ini_`v' = `v'[1]
}

foreach depVar in formal employed {
	foreach panel in allc noeduc { //option none also exists

		eststo clear
		replace `depVar' = `depVar'_`panel'
		replace wei = 1 / `depVar'_se_`panel'
		
		eststo r_basic: reg `depVar' educ_8p age_3* age_4* lpop ///
			i.year [aweight = wei], cluster(micror)
			
		eststo r_fe: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year [aweight = wei], abs(micror) cluster(micror)
			
		eststo r_inifor: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year i.year#c.ini_formal_noeduc ///
			[aweight = wei], abs(micror) cluster(micror)
			
		eststo r_statefe: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year##i.state i.year#c.ini_formal_noeduc ///
			[aweight = wei], abs(micror) cluster(micror)
			
		eststo r_sector: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year i.year#c.ini_formal_noeduc sh_* ///
			[aweight = wei], abs(micror) cluster(micror)
			
		eststo r_income: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year i.year#c.ini_formal_noeduc lincome ///
			[aweight = wei], abs(micror) cluster(micror)
			
		eststo r_iniall: areg `depVar' educ_8p age_3* age_4* lpop ///
			i.year i.year#c.ini_formal_noeduc ///
			i.year#c.ini_lincome i.year#c.ini_educ_8p ///
			i.year#c.ini_employed_noeduc ///
			[aweight = wei], abs(micror) cluster(micror)
		
		*Adding description
		estadd local Pop_Age_Controls "Yes" : _all
		estadd local Av_Earnings_Control "No" : _all
		estadd local Av_Earnings_Control "Yes", replace : r_income
		estadd local Microregion_FE "Yes" : _all
		estadd local Microregion_FE "No", replace : r_basic
		estadd local Time_FE "Yes" : _all
		estadd local Time_FE "No", replace : r_statefe
		estadd local State_Time_FE "No" : _all
		estadd local State_Time_FE "Yes", replace : r_statefe
		estadd local Sectoral_Shares "No" : _all
		estadd local Sectoral_Shares "Yes", replace : r_sector
		estadd local Ini_Formality_Time "Yes" : _all
		estadd local Ini_Formality_Time "No", replace : r_basic r_fe
		estadd local Other_ini_Time "No" : _all
		estadd local Other_ini_Time "Yes", replace : r_iniall
			
		esttab r_* using "../Output/Tab_`depVar'_`panel'.csv", ///
			se csv replace nonotes keep(educ*) ///
			star(* 0.1 ** 0.05 *** 0.01) compress nomtitles ///
			stats(Pop_Age_Controls Av_Earnings_Control ///
			Microregion_FE Time_FE State_Time_FE ///
			Sectoral_Shares Ini_Formality_Time Other_ini_Time ///
			N r2)
	}
}
	