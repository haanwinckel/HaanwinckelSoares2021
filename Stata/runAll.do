
forvalues i = 0/100 {
	do getMoments 2003 `i'
	do getMoments 2012 `i'
	do getNationalMoments `i'
}
