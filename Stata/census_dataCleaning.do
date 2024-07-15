set more off

local years 2010 1991 2000
local states     "RO AC AM RR PA AP TO MA PI CE RN PB PE AL SE BA MG ES RJ SP PR SC RS MS MT GO DF"
local cod_states "11 12 13 14 15 16 17 21 22 23 24 25 26 27 28 29 31 32 33 35 41 42 43 50 51 52 53"

local base_folder "A:\Dropbox\Pesquisa\Informality and Education in Brazil" //Folder where the clean data is saved
local md_folder "A:\Dropbox\Microdados\Censo\Microdados" //Folder with original microdata

*Create temp folder for datazoom_censo output
/*tempfile temp_f
local temp_f `"`=substr("`temp_f'",1,`=length("`temp_f'")-4')'"'
cap mkdir "`temp_f'" */
local temp_f "A:\Datazoom" //A temporary folder that will store files created by Data Zoom

foreach YR of local years {
	foreach ST of local states {
		*Won't work if file for that state/year already in place
		* or for state of TO, year 1980.
		capture confirm file "`temp_f'\person_census_`YR'_`ST'.dta"
		
		if _rc != 0 & (`YR' > 1980 | "`ST'" != "TO") {
			
			*First part:
			*Obtain data from census using Data Zoom dictionaries
			* and standardized definitions.
			
			capture confirm file ///
				"`temp_f'\CENSO`=substr("`YR'",3,2)'_`ST'_comp.dta"			
			if _rc != 0 {
				clear
				datazoom_censo, years( `YR' ) ///
					ufs( `ST' ) saving("`temp_f'") ///
					original(`md_folder') comp both
			}
			else {
				disp "Datazoom file already exists. Loading..."
				use "`temp_f'\CENSO`=substr("`YR'",3,2)'_`ST'_comp.dta", clear
			}
			rename ano year
			rename UF state
			
			*gender
			gen p_female = sexo == 0 if sexo != .
			
			*race
			rename raca p_race
			
			*Age
			rename idade p_age
			
			*Education variables
			gen p_educ_literate = alfabetizado==1
			if `YR' == 1970 {
				gen p_educ_gr1 = ( ///
					anos_estudo < 4) if anos_estudo != .
				gen p_educ_gr2 = (anos_estudo >= 4 & ///
					 anos_estudo < 8) if anos_estudo != .
				gen p_educ_gr3 = (anos_estudo >= 8 & ///
					anos_estudo < 11) if anos_estudo != .
				gen p_educ_gr4 = (anos_estudo >= 11 & ///
					anos_estudo < 15) if anos_estudo != .
				gen p_educ_gr5 = ( ///
					anos_estudo >= 15) if anos_estudo != .
			}
			else {
				quietly: tabulate anos_estudoC, generate(p_educ_gr)
			}
			
			*Employment/unemployment
			gen working_age = p_age >= 18 & p_age <= 65 if p_age != .
			if `YR' >= 2000 {
				gen p_workforce = working_age == 1 & ///
					(tomou_prov == 1 | ///
					trab_rem_sem == 1 | afast_trab_sem == 1)
				gen p_unemp = tomou_prov == 1 & ///
					(trab_rem_sem == 0 & afast_trab_sem == 0) ///
					if p_workforce == 1
				gen p_employed = trab_rem_sem == 1 | ///
					afast_trab_sem == 1
			}
			else {
				gen p_workforce = working_age == 1 & ///
					cond_ativ >= 0 & cond_ativ <= 2
				gen p_unemp_before2000 = cond_ativ != 0 ///
					if p_workforce == 1
				gen p_unemp_after2000 = 0
				gen p_employed = cond_ativ == 0
			}
			
			*keep if p_workforce == 1
			keep if freq_escola == 0 //Not going to school
			
			gen p_rural = sit_setor_C == 0			
			
			*Earnings
			gen p_income = rend_total
			
			*Formality
			if `YR' >= 2000 {
				gen p_formal = pos_ocup_sem == 1 if ///
					(pos_ocup_sem==1 | pos_ocup_sem==3)
			}
			else if `YR' == 1991 {
				*Info available for 1991, but not in Datazoom compat.
				*But due to modification in compat_censo91pess.ado,
				* line 559, this information is kept in the
				* compat file.
				gen p_formal = v0350 == 1 if v0350 == 1 | v0350 == 3
			}
			
			*Self-employment
			if `YR' >= 2000 {
				gen p_self = pos_ocup_sem == 6 if ///
					pos_ocup_sem != .
			}
			else {
				gen p_self = pos_ocup_hab == 7 if ///
					pos_ocup_hab != .
			}
			
			*Shares of workforce in each sector
			if `YR' == 1991 {
				gen ativ = ativ_hab
			}
			else if `YR' >= 2000 {
				capture drop ativ1991
				merge m:1 ativ2000 using ///
					"`base_folder'\Auxiliary Code\activityCrosswalk.dta"
				drop if _merge == 2
				gen ativ = ativ1991
			}
			else {
				gen ativ = .
			}
			*Agriculture, veg. extractive, livestock, fishery
			gen p_act_agr = ativ >= 1 & ativ <= 42
			*Mining (mineral extractive)
			gen p_act_mining = ativ >= 50 & ativ <= 59
			*Transformation industry
			gen p_act_transf = ativ >= 100 & ativ <= 300
			*Construction
			gen p_act_cnst = ativ == 340
			*Utilities
			gen p_act_util = ativ >= 351 & ativ <= 354
			*Commerce/trading
			gen p_act_comm = ativ >= 410 & ativ <= 424
			*Services (includes transp, finance, medical, educ)
			gen p_act_serv = ativ >= 451 & ativ <= 632
			*Govt
			gen p_act_govt = ativ >= 711 & ativ <= 801
			*Other
			gen p_act_other = ativ == 901 | ativ == 902 | ativ == 0
			foreach v of varlist p_act_* {
				replace `v' = . if ativ == .
			}			
			
			*migrant
			gen p_migr = t_mor_mun_80 <= 6 if ///
				p_age >= 10 & p_age != .
				
			*Survey weight
			rename peso_pess p_weight
			
			*Keep only variables of interest
			cap rename v1002 p_mesoregion
			cap rename v1003 p_microregion
			keep state year newamc9100 newamc0010 p_* id_dom ordem ///
				deflator conversor
			sort id_dom ordem	
			
			*Account for RO 1991 problems
			by id_dom ordem: gen check = _n
			drop if check > 1
			drop check
			
			save "`temp_f'\person_census_`YR'_`ST'", replace
		}
	}
}
clear
foreach YR of local years {
	foreach ST of local states {
		if `YR' > 1980 | "`ST'" != "TO" {
			append using "`temp_f'\person_census_`YR'_`ST'"
		}
	}
}

foreach v of varlist p_income {
	replace `v' = `v' / (deflator*conversor)
}
drop deflator conversor

save "`base_folder'\census_person.dta", replace

