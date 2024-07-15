

clear all
set more off
set matsize 4000
set maxvar 10000
set max_memory 12g

local base_folder /// Where the clean data file is stored
	"."


use "`base_folder'\census_person.dta", clear

*Keep only men, 20-50 years old, urban areas
drop if p_female == 1
drop if p_age < 20 | p_age > 50
keep if p_rural == 0

destring p_microregion, replace	
destring p_mesoregion, replace
replace p_mesoregion = 100*state + p_mesoregion
replace p_microregion = 100*p_mesoregion + p_microregion
sort newamc9100 newamc0010 year
unique p_microregion
by newamc9100 newamc0010: replace p_microregion = p_microregion[_N]
unique p_microregion

**Create microregions with constant areas over time
*Will first use a smaller dataset to create a crosswalk
*from old microregion levels to the new microregion

preserve
bys newamc9100 newamc0010 year: keep if _n == _N
*Get stable municipalities, from 1991 to 2010
encode newamc9100, gen(code1)
gen code2 = newamc0010
sort code1 code2
forvalues i=1/3 {
	bysort code1: replace code2 = code2[1]
	sort code2 code1
	bysort code2: replace code1 = code1[1]
	sort code1 code2
	by code1: egen v2 = sd(code2)
	cap assert v2 == 0 | v2 == ., fast
	if _rc == 0 {
		continue, break
	}
	drop v2
}
assert v2 == 0 | v2 == .

*Using those, replace values for microregions 
gen new_micror = p_microregion
egen stableGeo = group(code1)
unique(p_microregion), by(stableGeo) gen(num_micror)
summ stableGeo
local numGeo = r(max)
quietly {
	forvalues i = 1/`numGeo' {
		*Check if more than one microregion
		summ num_micror if stableGeo == `i'
		if r(mean) > 1 { //More than one microregion:
			levelsof p_microregion if stableGeo == `i', ///
				local(list)
			local baseMR : word 1 of `list'
			foreach mr of local list {
				replace new_micror = `baseMR' ///
					if p_microregion == `mr'
			}
		}
	}
}
keep p_microregion new_micror
duplicates drop
tempfile crosswalk
save `crosswalk'
unique p_microregion
unique new_micror

restore
merge m:1 p_microregion using `crosswalk'

egen micror = group(new_micror)
egen mesor = group(p_mesoregion)

summ micror
local num_mr = r(max)

gen p_age2 = p_age^2
gen p_age3 = p_age^3
gen p_age4 = p_age^4


bysort year micror: gen first = (_n == 1)

*Make sure all microregions are included in the regression
fvset base none micror
drop p_educ_gr1

foreach choice in allc noeduc { 
	
	if "`choice'" == "allc" {
		local controls "p_educ* p_age* i.p_race"
	}
	else {
		local controls "p_age* i.p_race"
	}
	
	foreach depVar in formal employed {
		gen `depVar'_`choice' = .
		gen `depVar'_se_`choice' = .
		foreach yr in 1991 2000 2010 {
			reg p_`depVar' `controls' i.micror if year == `yr', ///
				robust nocons
			forvalues iMicror = 1/`num_mr' {
				replace `depVar'_`choice' = _b[`iMicror'.micror] ///
					if year == `yr' & micror == `iMicror' & first == 1
				replace `depVar'_se_`choice' = _se[`iMicror'.micror] ///
					if year == `yr' & micror == `iMicror' & first == 1
			}
		}
	}
}

*Create variables to be aggregated

gen educ_8p = p_educ_gr3 + p_educ_gr4 + p_educ_gr5
label var educ_8p "Share 8+ yrs"

gen income = p_income / 1000
label var income "Av. income"

gen age_20_29 = p_age < 30
gen age_30_39 = p_age >= 30 & p_age < 40
gen age_40p = p_age >= 40

gen pop = 1
label var pop "Size workforce"

foreach sct in agr mining transf cnst util ///
	comm serv govt other {
	rename p_act_`sct' sh_`sct'
}

collapse (mean) employed_* formal_* educ_*p income age_* ///
	sh_* (sum) pop [pweight = p_weight], ///
	by(year state mesor micror)

gen lpop = log(pop)
gen lincome = log(income)

save "census_collapsed_data.dta", replace
