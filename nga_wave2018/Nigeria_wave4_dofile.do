










clear



global Nigeria_GHS_W4_raw_data 		"C:\Users\obine\Music\Documents\Smallholder lsms STATA\NGA_2018_GHSP-W4_v03_M_Stata12 (1)"
global Nigeria_GHS_W4_created_data  "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data/try_ng\nga_wave2018"




********************************************************************************
* AG FILTER *
********************************************************************************

use "${Nigeria_GHS_W4_raw_data}/sect11a_plantingw4.dta" , clear

keep hhid ag1
rename (ag1) (ag_rainy_18)
save  "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", replace



*merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

*keep if ag_rainy_18==1






********************************************************************************
* WEIGHTS *
********************************************************************************

use "${Nigeria_GHS_W4_raw_data}/secta_plantingw4.dta" , clear
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
gen rural = (sector==2)
lab var rural "1= Rural"
keep hhid zone state lga ea wt_wave4 rural
ren wt_wave4 weight
collapse (max) weight, by (hhid)
save  "${Nigeria_GHS_W4_created_data}/weight.dta", replace




************************
*Geodata Variables
************************

use "${Nigeria_GHS_W4_raw_data}\nga_plotgeovariables_y4.dta", clear
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

encode srtmslp_nga, gen( plot_slope)
ren srtm_nga  plot_elevation
ren twi_nw  plot_wetness

tab1 plot_slope plot_elevation plot_wetness, missing

collapse (sum) plot_slope plot_elevation plot_wetness, by (hhid)
sort hhid
la var plot_slope "slope of plot"
la var plot_elevation "Elevation of plot"
la var plot_wetness "Potential wetness index of plot"
save "${Nigeria_GHS_W4_created_data}\geodata_2018.dta", replace









*********************************************** 
*Used Fertilizer
***********************************************
use "${Nigeria_GHS_W4_raw_data}\secta11c2_harvestw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
//plot level variables..... no question about agency where they purchased it
 

*s11dq1a      1= if used fertilizer on plot
*s11c2q36_1   1= if used npk on plot
*s11c2q36_2   1= if used urea on plot
*s11c2q36_99  1= if used other fert on plot
*s11c2q37a    qty of npk used 
*s11c2q37b    units of npk used 
*s11c2q37a_conv conversion factor
*s11c2q38a    qty of urea used
*s11c2q38b    units of urea used 
*s11c2q38a_conv  coversion factor
*s11c2q39a    qty of other fert
*s11c2q39b    units of other fert
*s11c2q39a_conv  conversion factor
*s11dq36      1= if used org fert on plot
*s11dq37a     qty of org used
*s11dq37b     units of org used
*s11c2q37_conv   conversion


*br s11c2q38a s11c2q38b s11c2q38a_conv

*****Coversion of fertilizer's units into kilogram using 

gen fert1 = s11c2q37a*s11c2q37a_conv if s11dq1a ==1
gen fert2 = s11c2q38a*s11c2q38a_conv if s11dq1a ==1
gen fert3 = s11c2q39a*s11c2q39a_conv  if s11dq1a ==1

*br s11c2q37a s11c2q37a_conv fert1 s11c2q38a s11c2q38a_conv fert2 s11c2q39a s11c2q39a_conv fert3 s11dq1a if s11dq1a ==1
****generate the total qty*************
egen used_fert = rowtotal(fert1 fert2 fert3)
sum  used_fert , detail

replace used_fert  = 500 if used_fert > 500
tab  used_fert 
sum  used_fert , detail

***************
*organic fertilizer
***************
gen org_fert = (s11dq36==1)
tab org_fert, missing


collapse (sum) used_fert (max) org_fert , by(hhid)



merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1


sum  used_fert , detail
************winzonrizing total_qty
foreach v of varlist  used_fert  {
	_pctile `v' [aw=weight] , p(1 90) 
	gen `v'_w=`v'
	*replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 10%"
}


tab used_fert
tab used_fert_w, missing
sum used_fert used_fert_w, detail

keep hhid used_fert_w org_fert 

la var org_fert "1= if used organic fertilizer"
label var used_fert_w "quantity of inorganic fertilizer used in kg"
sort hhid
save "${Nigeria_GHS_W4_created_data}\total_qty_2018.dta", replace



****************************
*Subsidized Fertilizer
****************************

use "${Nigeria_GHS_W4_raw_data}\secta11c3_harvestw4.dta",clear   
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
 //household level variables......... distance to purchase was also asked

*inputid  1= org fert| 2-4 inorg fert
*s11c3q2  1= hhold purchased inputs
*s11c3q4a qty of purchased inputs
*s11c3q4b units of purchased inputs
*s11c3q4_conv conversion factor
*s11c3q5  cost of inputs
*s11c3q6b institute of purchased
*s11c3q7  distance to institute (km)


******conversion to kg

gen input_kg = s11c3q4a*s11c3q4_conv if inputid >=2 & inputid <=4

*br s11c3q4a s11c3q4b s11c3q4_conv input_kg

***getting the qty for inorg fertilizer

gen inorg_fert = input_kg if inputid >=2 & inputid <=4
*br input_kg inputid inorg_fert
tab inorg_fert






************
*Getting total subsidy_dummy
**********

*****I am using Agricultural input dealer and govt extension officer on the assumption that they sell subsidized fertilizer

gen subsidy_qty = inorg_fert if s11c3q6b ==1 | s11c3q6b ==3
tab subsidy_qty
tab subsidy_qty,missing
sum subsidy_qty,detail


gen subsidy_dummy = 0
replace subsidy_dummy = 1 if s11c3q6b ==1 | s11c3q6b ==3
tab subsidy_dummy, missing


collapse (sum)subsidy_qty (max) subsidy_dummy, by (hhid)


merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

************winzonrizing subsidy_qty
foreach v of varlist  subsidy_qty  {
	_pctile `v' [aw=weight] , p(1 99) 
	gen `v'_w=`v'
	*replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 1%"
}

tab subsidy_qty
tab subsidy_qty_w, missing
sum subsidy_qty subsidy_qty_w, detail







keep hhid  subsidy_qty_w subsidy_dummy



label var subsidy_qty "Quantity of Fertilizer Purchased in kg"
label var subsidy_dummy "=1 if acquired any subsidied fertilizer"
save "${Nigeria_GHS_W4_created_data}\subsidized_fert_2018.dta", replace

*********************************************** 
*Purchased Fertilizer
***********************************************

use "${Nigeria_GHS_W4_raw_data}\secta11c3_harvestw4.dta",clear   
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
 //household level variables......... distance to purchase was also asked

*inputid  1= org fert| 2-4 inorg fert
*s11c3q2  1= hhold purchased inputs
*s11c3q4a qty of purchased inputs
*s11c3q4b units of purchased inputs
*s11c3q4_conv conversion factor
*s11c3q5  cost of inputs
*s11c3q6b institute of purchased
*s11c3q7  distance to institute (km)


******conversion to kg

gen input_kg = s11c3q4a*s11c3q4_conv if inputid >=2 & inputid <=4

*br s11c3q4a s11c3q4b s11c3q4_conv input_kg

***getting the qty for inorg fertilizer

gen inorg_fert = input_kg if inputid >=2 & inputid <=4
*br input_kg inputid inorg_fert
tab inorg_fert

gen total_qty = input_kg if inputid >=2 & inputid <=4 & s11c3q6b== 6


gen cost_fert = s11c3q5 if inputid >=2 & inputid <=4
gen cost_fert_real = cost_fert if s11c3q6b== 6
*br s11c3q5 inputid cost_fert

gen tpricefert = cost_fert_real/total_qty
tab tpricefert

gen tpricefert_cens = tpricefert 
replace tpricefert_cens = 900 if tpricefert_cens > 900 & tpricefert_cens < .   //winzonrizing bottom 1%
replace tpricefert_cens = 60 if tpricefert_cens < 60
tab tpricefert_cens, missing //winzonrizing top 1%


egen medianfert_pr_ea = median(tpricefert_cens), by (ea)
egen medianfert_pr_lga = median(tpricefert_cens), by (lga)
egen medianfert_pr_state = median(tpricefert_cens), by (state)
egen medianfert_pr_zone = median(tpricefert_cens), by (zone)



egen num_fert_pr_ea = count(tpricefert_cens), by (ea)
egen num_fert_pr_lga = count(tpricefert_cens), by (lga)
egen num_fert_pr_state = count(tpricefert_cens), by (state)
egen num_fert_pr_zone = count(tpricefert_cens), by (zone)



tab medianfert_pr_ea
tab medianfert_pr_lga
tab medianfert_pr_state
tab medianfert_pr_zone



tab num_fert_pr_ea
tab num_fert_pr_lga
tab num_fert_pr_state
tab num_fert_pr_zone

gen tpricefert_cens_mrk = tpricefert_cens

replace tpricefert_cens_mrk = medianfert_pr_ea if tpricefert_cens_mrk ==. & num_fert_pr_ea >= 7

tab tpricefert_cens_mrk,missing


replace tpricefert_cens_mrk = medianfert_pr_lga if tpricefert_cens_mrk ==. & num_fert_pr_lga >= 7

tab tpricefert_cens_mrk,missing



replace tpricefert_cens_mrk = medianfert_pr_state if tpricefert_cens_mrk ==. & num_fert_pr_state >= 7

tab tpricefert_cens_mrk,missing


replace tpricefert_cens_mrk = medianfert_pr_zone if tpricefert_cens_mrk ==. & num_fert_pr_zone >= 7

tab tpricefert_cens_mrk,missing


********Distance to institute of purchased fertilizer
gen distance = s11c3q7 if inputid >=2 & inputid <=4
tab distance
replace distance = . if distance== 0
tab distance

egen medianfert_dist_ea = median(distance), by (ea)
egen medianfert_dist_lga = median(distance), by (lga)
egen medianfert_dist_state = median(distance), by (state)
egen medianfert_dist_zone = median(distance), by (zone)
egen medianfert_dist_sector = median(distance), by (sector)


egen num_fert_dist_ea = count(distance), by (ea)
egen num_fert_dist_lga = count(distance), by (lga)
egen num_fert_dist_state = count(distance), by (state)
egen num_fert_dist_zone = count(distance), by (zone)
egen num_fert_dist_sector = count(distance), by (sector)


tab medianfert_dist_ea
tab medianfert_dist_lga
tab medianfert_dist_state
tab medianfert_dist_zone



tab num_fert_dist_ea
tab num_fert_dist_lga
tab num_fert_dist_state
tab num_fert_dist_zone

gen mrk_dist = distance

replace mrk_dist = medianfert_dist_ea if mrk_dist ==. & num_fert_dist_ea >= 20

tab mrk_dist,missing


replace mrk_dist = medianfert_dist_lga if mrk_dist ==. & num_fert_dist_lga >= 20

tab mrk_dist,missing



replace mrk_dist = medianfert_dist_state if mrk_dist ==. & num_fert_dist_state >= 20

tab mrk_dist,missing


replace mrk_dist = medianfert_dist_zone if mrk_dist ==. & num_fert_dist_zone >= 20

tab mrk_dist,missing
replace mrk_dist = medianfert_dist_sector if mrk_dist ==. & num_fert_dist_sector >= 20

tab mrk_dist,missing




collapse zone lga sector ea (sum) total_qty (max)   mrk_dist  tpricefert_cens_mrk, by(hhid)


merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1


************winzonrizing total_qty
foreach v of varlist  total_qty  {
	_pctile `v' [aw=weight] , p(1 95) 
	gen `v'_w=`v'
	*replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 1%"
}


tab total_qty
tab total_qty_w, missing
sum total_qty total_qty_w, detail


************winzonrizing fertilizer distance
foreach v of varlist  mrk_dist  {
	_pctile `v' [aw=weight] , p(1 99) 
	gen `v'_w=`v'
	*replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 1%"
}

tab mrk_dist
tab mrk_dist_w, missing
sum mrk_dist mrk_dist_w, detail

gen rea_tpricefert_cens_mrk = tpricefert_cens_mrk
gen real_tpricefert_cens_mrk = rea_tpricefert_cens_mrk
tab real_tpricefert_cens_mrk
sum real_tpricefert_cens_mrk, detail


keep hhid zone lga sector ea total_qty_w mrk_dist_w real_tpricefert_cens_mrk







la var mrk_dist_w "Distance from farm to where you purchased inorg fertilizer"

label var real_tpricefert_cens_mrk "price of commercial fertilizer purchased in naira"
sort hhid
save "${Nigeria_GHS_W4_created_data}\purchased_fert_2018.dta", replace




************************************************
*Savings 
************************************************



use "${Nigeria_GHS_W4_raw_data}\sect4a1_plantingw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

*s4aq1 1= have a bank acccount
*s4aq8 1= used commmercial bank savings
*s4aq10 1=  used informal savings

ren s4aq1 formal_bank
tab formal_bank, missing
replace formal_bank =0 if formal_bank ==2 | formal_bank ==.
tab formal_bank, nolabel
tab formal_bank,missing

 ren s4aq8 formal_save 
 tab formal_save, missing
 replace formal_save =0 if formal_save ==2 | formal_save ==.
 tab formal_save, missing

 ren s4aq10 informal_save 
 tab informal_save, missing
 replace informal_save =0 if informal_save ==2 | informal_save ==.
 tab informal_save, missing

 collapse (max) formal_bank  formal_save  informal_save, by (hhid)
 la var formal_bank "=1 if respondent have an account in bank"
 la var formal_save "=1 if used formal saving group"
 la var informal_save "=1 if used informal saving group"
save "${Nigeria_GHS_W4_created_data}\savings_2018.dta", replace



*******************************************************
*Credit access 
*******************************************************

use "${Nigeria_GHS_W4_raw_data}\sect4c2_plantingw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
*s4cq2b   type of loan lenders (=<4 formal banks)
*s4cq20   <=2 if loan was approved
 
tab s4cq2b
label list S4CQ20
 gen formal_credit  =1 if s4cq20<=2 & s4cq2b <=4
 tab formal_credit,missing
 replace formal_credit =0 if formal_credit ==.
 tab formal_credit,missing
 
 
 gen informal_credit =1 if s4cq20<=2 & s4cq2b >=5
 tab informal_credit,missing
replace informal_credit =0 if informal_credit ==.
 tab informal_credit,missing


 collapse (max) formal_credit  informal_credit, by (hhid)
 la var formal_credit "=1 if borrowed from formal credit group"
 la var informal_credit "=1 if borrowed from informal credit group"
save "${Nigeria_GHS_W4_created_data}\credit_access_2018.dta", replace





******************************* 
*Extension Visit 
*******************************


use "${Nigeria_GHS_W4_raw_data}\sect11l1_plantingw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

ren s11l1q1 ext_acess 

tab ext_acess, missing
tab ext_acess, nolabel

replace ext_acess = 0 if ext_acess==2 | ext_acess==.
tab ext_acess, missing
collapse (max) ext_acess, by (hhid)
la var ext_acess "=1 if received advise from extension services"
save "${Nigeria_GHS_W4_created_data}\extension_access_2018.dta", replace





*****************************
*Community 
****************************

use "${Nigeria_GHS_W4_raw_data}\sectc2_harvestw4.dta", clear
*is_cd  222 for market infrastructure
*c2q3  distance to infrastructure in km

gen mrk_dist1 = c2q3 if is_cd==222
tab mrk_dist1 if is_cd==222, missing
egen median_lga = median(mrk_dist1), by (zone state lga)
egen median_state = median(mrk_dist1), by (zone state)
egen median_zone = median(mrk_dist1), by (zone)


replace mrk_dist1 =0 if is_cd==222 & mrk_dist1==. & c2q1==1
tab mrk_dist1 if is_cd==222, missing

replace mrk_dist1 = median_lga if mrk_dist1==. & is_cd==222
replace mrk_dist1 = median_state if mrk_dist1==. & is_cd==222
replace mrk_dist1 = median_zone if mrk_dist1==. & is_cd==222
tab mrk_dist1 if is_cd==222, missing

replace mrk_dist1= 45 if mrk_dist1>=45 & mrk_dist1<. & is_cd==222
tab mrk_dist1 if is_cd==222, missing

sort zone state ea
collapse (max) median_lga median_state median_zone mrk_dist1, by (zone state lga sector ea)
replace mrk_dist1 = median_lga if mrk_dist1 ==.
tab mrk_dist1, missing
replace mrk_dist1 = median_state if mrk_dist1 ==.
tab mrk_dist1, missing
replace mrk_dist1 = median_zone if mrk_dist1 ==.
tab mrk_dist1, missing
la var mrk_dist1 "=distance to the market"

save "${Nigeria_GHS_W4_created_data}\market_distance.dta", replace 




*********************************
*Demographics 
*********************************



use "${Nigeria_GHS_W4_raw_data}\sect1_plantingw4.dta",clear 


merge 1:1 hhid indiv using "${Nigeria_GHS_W4_raw_data}\sect2_harvestw4.dta", gen(household)

merge m:1 zone state lga sector ea using "${Nigeria_GHS_W4_created_data}\market_distance.dta", keepusing (median_lga median_state median_zone mrk_dist1)

merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
**************
*market distance
*************






*s1q2 sex
*s1q3 relationship with hhead (1= head)
*s1q6 age (in years)
sort hhid indiv 
 
gen num_mem  = 1


******** female head****

gen femhead  = 0
replace femhead = 1 if s1q2== 2 & s1q3==1
tab femhead,missing

********Age of HHead***********
ren s1q6 hh_age
gen hh_headage  = hh_age if s1q3==1

tab hh_headage

replace hh_headage = 100 if hh_headage > 100 & hh_headage < .
tab hh_headage
tab hh_headage, missing


************generating the median age**************

egen medianhh_pr_ea = median(hh_headage), by (ea)

egen medianhh_pr_lga = median(hh_headage), by (lga)

egen num_hh_pr_ea = count(hh_headage), by (ea)

egen num_hh_pr_lga = count(hh_headage), by (lga)

egen medianhh_pr_state = median(hh_headage), by (state)
egen num_hh_pr_state = count(hh_headage), by (state)

egen medianhh_pr_zone = median(hh_headage), by (zone)
egen num_hh_pr_zone = count(hh_headage), by (zone)


tab medianhh_pr_ea
tab medianhh_pr_lga
tab medianhh_pr_state
tab medianhh_pr_zone



tab num_hh_pr_ea
tab num_hh_pr_lga
tab num_hh_pr_state
tab num_hh_pr_zone



replace hh_headage = medianhh_pr_ea if hh_headage ==. & num_hh_pr_ea >= 30

tab hh_headage,missing


replace hh_headage = medianhh_pr_lga if hh_headage ==. & num_hh_pr_lga >= 30

tab hh_headage,missing



replace hh_headage = medianhh_pr_state if hh_headage ==. & num_hh_pr_state >= 30

tab hh_headage,missing


replace hh_headage = medianhh_pr_zone if hh_headage ==. & num_hh_pr_zone >= 30

tab hh_headage,missing

sum hh_headage, detail



********************Education****************************************************
*s2aq6 attend school
*s2aq9 highest level of edu completed
*s1q3 relationship with hhead (1= head)

ren  s2aq6 attend_sch 
tab attend_sch
replace attend_sch = 0 if attend_sch ==2
tab attend_sch, nolabel
*tab s1q4 if s2q7==.

replace s2aq9= 0 if attend_sch==0
tab s2aq9
tab s1q3 if _merge==1

tab s2aq9 if s1q3==1
replace s2aq9 = 16 if s2aq9==. &  s1q3==1

*** Education Dummy Variable*****

 label list S2AQ9

gen pry_edu  = 1 if s2aq9 >= 1 & s2aq9 < 16 & s1q3==1
gen finish_pry = 1 if s2aq9 >= 16 & s2aq9 < 26 & s1q3==1
gen finish_sec  = 1 if s2aq9 >= 26 & s2aq9 & s1q3==1
replace finish_sec  =0 if s2aq9==51 | s2aq9==52 & s1q3==1

replace pry_edu =0 if pry_edu==. & s1q3==1
replace finish_pry  =0 if finish_pry==. & s1q3==1
replace finish_sec =0 if finish_sec==. & s1q3==1
tab pry_edu if s1q3==1 , missing
tab finish_pry if s1q3==1 , missing 
tab finish_sec if s1q3==1 , missing

collapse (sum) num_mem (max) mrk_dist1 hh_headage femhead attend_sch  pry_edu finish_pry finish_sec, by (hhid)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1


tab mrk_dist1, missing
************winzonrizing distance to market
egen median_hhid = median(mrk_dist1), by (hhid)
replace mrk_dist1 = median_hhid if mrk_dist1==.
tab mrk_dist1, missing




foreach v of varlist  mrk_dist1  {
	_pctile `v' [aw=weight] , p(1 99) 
	gen `v'_w=`v'
	*replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 5%"
}


tab mrk_dist1
tab mrk_dist1_w, missing
sum mrk_dist1 mrk_dist1_w, detail


keep hhid mrk_dist1_w num_mem femhead hh_headage attend_sch pry_edu finish_pry finish_sec

tab attend_sch, missing
egen mid_attend= median(attend_sch)
replace attend_sch = mid_attend if attend_sch==.

tab pry_edu, missing
tab finish_pry, missing
tab finish_sec, missing

egen mid_pry_edu= median(pry_edu)
egen mid_finish_pry= median(finish_pry)
egen mid_finish_sec= median(finish_sec)

replace pry_edu = mid_pry_edu if pry_edu==.
replace finish_pry = mid_finish_pry if finish_pry==.
replace finish_sec = mid_finish_sec if finish_sec==.



la var num_mem "household size"
la var mrk_dist1_w "distance to the nearest market in km"
la var femhead  "=1 if head is female"
la var hh_headage "age of household head in years"
la var attend_sch"=1 if respondent attended school"
la var pry_edu  "=1 if household head attended pry school"
la var finish_pry "=1 if household head finished pry school"
la var finish_sec "=1 if household head finished sec school"
save "${Nigeria_GHS_W4_created_data}\demographics_2018.dta", replace

********************************* 
*Labor Age 
*********************************
use "${Nigeria_GHS_W4_raw_data}\sect1_plantingw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
ren s1q6 hh_age

gen worker = 1
replace worker = 0 if hh_age < 15 | hh_age > 65

tab worker,missing
sort hhid
collapse (sum) worker, by (hhid)
la var worker "number of members age 15 and older and less than 65"
sort hhid

save "${Nigeria_GHS_W4_created_data}\laborage_2018.dta", replace


********************************
*Safety Net
********************************

use "${Nigeria_GHS_W4_raw_data}\sect14a_harvestw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
*s14q1a__1 1= received cash 
*s14q1a__2 1= received food 
*s14q1a__3 1= received other kinds
**s14q1a__4 1= received from institutes 

gen safety_net  =1 if s14q1a__1==1 | s14q1a__2==1 | s14q1a__3==1 | s14q1a__4==1

replace safety_net =0 if safety_net==.
tab safety_net,missing
collapse (max) safety_net, by (hhid)
tab safety_net
la var safety_net "=1 if received cash transfer, cash for work, food for work or other assistance"
save "${Nigeria_GHS_W4_created_data}\safety_net_2018.dta", replace




















***********************************
*Food Prices from Community
*********************************
use "${Nigeria_GHS_W4_raw_data}\sectc2_plantingw4.dta", clear
*rice is 13, maize is 16

*br if item_cd == 20
*br if item_cd ==20 & c2q2==1
tab c2q3 if item_cd ==20 & c2q2==1
tab c2q2 if item_cd==20
tab c2q3 if item_cd ==20
tab c2q3 if item_cd ==13





gen conversion =1
tab conversion, missing
gen food_size=1 //This makes it easy for me to copy-paste existing code rather than having to write a new block
replace conversion = food_size*2.696 if c2q2 == 11
replace conversion = food_size*0.001 if  c2q2 == 2
replace conversion = food_size*0.175 if  c2q2 == 12		
replace conversion = food_size*0.23 if  c2q2 == 13
replace conversion = food_size*1.5 if  c2q2 == 20 |c2q2 == 21  |c2q2 == 30  |c2q2 == 31 	
replace conversion = food_size*0.35 if  c2q2 == 40 
replace conversion = food_size*0.70 if  c2q2 == 41
replace conversion = food_size*3.00 if  c2q2 == 51  |c2q2 == 52 
replace conversion = food_size*0.718 if  c2q2 == 70	 |c2q2 == 71  |c2q2 == 72
replace conversion = food_size*1.615 if  c2q2 == 80  |c2q2 == 81  |c2q2 == 82
replace conversion = food_size*1.135 if   c2q2 == 90  |c2q2 == 91  |c2q2 == 92
				
tab conversion, missing	



gen maize_price= c2q3* conversion if item_cd==20
tab maize_price

sum maize_price, detail


tab maize_price,missing
sum maize_price,detail
tab maize_price

replace maize_price = 500 if maize_price >500 & maize_price<.  //bottom 5%
replace maize_price = 50 if maize_price< 50       ////top 5%



egen median_pr_ea = median(maize_price), by (ea)
egen median_pr_lga = median(maize_price), by (lga)
egen median_pr_state = median(maize_price), by (state)
egen median_pr_zone = median(maize_price), by (zone)

egen num_pr_ea = count(maize_price), by (ea)
egen num_pr_lga = count(maize_price), by (lga)
egen num_pr_state = count(maize_price), by (state)
egen num_pr_zone = count(maize_price), by (zone)

tab num_pr_ea
tab num_pr_lga
tab num_pr_state
tab num_pr_zone


gen maize_price_mr = maize_price

replace maize_price_mr = median_pr_ea if maize_price_mr==. & num_pr_ea>=2
tab maize_price_mr,missing

replace maize_price_mr = median_pr_lga if maize_price_mr==. & num_pr_lga>=2
tab maize_price_mr,missing

replace maize_price_mr = median_pr_state if maize_price_mr==. & num_pr_state>=2
tab maize_price_mr,missing

replace maize_price_mr = median_pr_zone if maize_price_mr==. & num_pr_zone>=2
tab maize_price_mr,missing






****************rice

gen rice_price = c2q3* conversion if item_cd==13

sum rice_price,detail
tab rice_price

replace rice_price = 1500 if rice_price >1500 & rice_price<.   //bottom 5%
replace rice_price = 10 if rice_price< 10  //top 5%
tab rice_price,missing



egen median_rice_ea = median(rice_price), by (ea)
egen median_rice_lga = median(rice_price), by (lga)
egen median_rice_state = median(rice_price), by (state)
egen median_rice_zone = median(rice_price), by (zone)

egen num_rice_ea = count(rice_price), by (ea)
egen num_rice_lga = count(rice_price), by (lga)
egen num_rice_state = count(rice_price), by (state)
egen num_rice_zone = count(rice_price), by (zone)

tab num_rice_ea
tab num_rice_lga
tab num_rice_state
tab num_rice_zone


gen rice_price_mr = rice_price

replace rice_price_mr = median_rice_ea if rice_price_mr==. & num_rice_ea>=7
tab rice_price_mr,missing

replace rice_price_mr = median_rice_lga if rice_price_mr==. & num_rice_lga>=7
tab rice_price_mr,missing

replace rice_price_mr = median_rice_state if rice_price_mr==. & num_rice_state>=7
tab rice_price_mr,missing

replace rice_price_mr = median_rice_zone if rice_price_mr==. & num_rice_zone>=7
tab rice_price_mr,missing


sort zone state ea
collapse (max) maize_price_mr rice_price_mr  median_pr_lga median_pr_state median_pr_zone median_pr_ea , by (zone state lga sector ea)


save "${Nigeria_GHS_W4_created_data}\food_prices.dta", replace




**************
*Net Buyers and Sellers
***************
use "${Nigeria_GHS_W4_raw_data}\sect7b_plantingw4.dta", clear
merge m:1 zone state lga sector ea using "${Nigeria_GHS_W4_created_data}\food_prices.dta", keepusing (median_pr_ea median_pr_lga median_pr_state median_pr_zone maize_price_mr rice_price_mr)
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1







**************
*maize price
*************
//missing values persists even after i did this
replace maize_price_mr = median_pr_ea if maize_price_mr==.
tab maize_price_mr, missing

replace maize_price_mr = median_pr_lga if maize_price_mr==.
tab maize_price_mr, missing

replace maize_price_mr = median_pr_state if maize_price_mr==.
tab maize_price_mr, missing

replace maize_price_mr = median_pr_zone if maize_price_mr==.
tab maize_price_mr, missing


tab rice_price_mr, missing

*s7bq5a from purchases
*s7bq6a from own production

tab s7bq5a
tab s7bq6a

replace s7bq5a = 0 if s7bq5a<=0 |s7bq5a==.
tab s7bq5a,missing
replace s7bq6a = 0 if s7bq6a<=0 |s7bq6a==.
tab s7bq6a,missing

gen net_seller = 1 if s7bq6a > s7bq5a
tab net_seller,missing
replace net_seller=0 if net_seller==.
tab net_seller,missing

gen net_buyer = 1 if s7bq6a < s7bq5a
tab net_buyer,missing
replace net_buyer=0 if net_buyer==.
tab net_buyer,missing


collapse zone (max) maize_price_mr rice_price_mr net_seller net_buyer, by(hhid)
gen rea_maize_price_mr = maize_price_mr
gen real_maize_price_mr = rea_maize_price_mr
tab real_maize_price_mr
sum real_maize_price_mr, detail
gen rea_rice_price_mr = rice_price_mr
gen real_rice_price_mr = rea_rice_price_mr
tab real_rice_price_mr
sum real_rice_price_mr, detail

label var real_maize_price_mr "commercial price of maize in naira"
label var real_rice_price_mr "commercial price of rice in naira"
la var net_seller "1= if respondent is a net seller"
la var net_buyer "1= if respondent is a net buyer"
sort hhid
save "${Nigeria_GHS_W4_created_data}\food_prices_2018.dta", replace





*****************************
*Household Assests
****************************


use "${Nigeria_GHS_W4_raw_data}\sect5_plantingw4.dta",clear 
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1
sort hhid item_cd

*s5q1 qty of item
*s5q4 value of item

gen hhasset_value  = s5q4*s5q1
tab hhasset_value,missing
sum hhasset_value,detail

/*
replace hhasset_value = 1000000 if hhasset_value > 2000000 & hhasset_value <.
replace hhasset_value = 200 if hhasset_value <200
replace hhasset_value = 0 if hhasset_value ==.
*/
sum hhasset_value, detail



collapse (sum) hhasset_value, by (hhid)
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1


foreach v of varlist  hhasset_value  {
	_pctile `v' [aw=weight] , p(1 99) 
	gen `v'_w=`v'
	replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top & bottom 5%"
}


tab hhasset_value
tab hhasset_value_w, missing
sum hhasset_value hhasset_value_w, detail

gen rea_hhvalue = hhasset_value_w
gen real_hhvalue= rea_hhvalue/1000
sum hhasset_value_w real_hhvalue, detail


keep  hhid real_hhvalue






la var real_hhvalue "total value of household asset"
save "${Nigeria_GHS_W4_created_data}\household_asset_2018.dta", replace





 ********************************************************************************
* PLOT AREAS *
********************************************************************************
clear
 
 
 
 *************** Plot Size **********************

//ALT IMPORTANT NOTE: As of W4, the implied area conversions for farmer estimated units (including hectares) are markedly different from previous waves. I recommend excluding plots that do not have GPS measured areas from any area-based productivity estimates.
use "${Nigeria_GHS_W4_raw_data}/sect11a1_plantingw4.dta", clear
*merging in planting section to get cultivated status
merge 1:1 hhid plotid using "${Nigeria_GHS_W4_raw_data}/sect11b1_plantingw4.dta", nogen
*merging in harvest section to get areas for new plots
merge 1:1 hhid plotid using "${Nigeria_GHS_W4_raw_data}/secta1_harvestw4.dta", gen(plot_merge)
 merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

ren s11aq4aa area_size
ren s11aq4b area_unit
ren sa1q11 area_size2 //GPS measurement, no units in file
//ren sa1q9b area_unit2 //Not in file
ren s11aq4c area_meas_sqm
//ren sa1q9c area_meas_sqm2
gen cultivate = s11b1q27 ==1 


gen field_size= area_size if area_unit==6
replace field_size = area_size*0.0667 if area_unit==4									//reported in plots
replace field_size = area_size*0.404686 if area_unit==5		    						//reported in acres
replace field_size = area_size*0.0001 if area_unit==7									//reported in square meters

replace field_size = area_size*0.00012 if area_unit==1 & zone==1						//reported in heaps
replace field_size = area_size*0.00016 if area_unit==1 & zone==2
replace field_size = area_size*0.00011 if area_unit==1 & zone==3
replace field_size = area_size*0.00019 if area_unit==1 & zone==4
replace field_size = area_size*0.00021 if area_unit==1 & zone==5
replace field_size = area_size*0.00012 if area_unit==1 & zone==6

replace field_size = area_size*0.0027 if area_unit==2 & zone==1							//reported in ridges
replace field_size = area_size*0.004 if area_unit==2 & zone==2
replace field_size = area_size*0.00494 if area_unit==2 & zone==3
replace field_size = area_size*0.0023 if area_unit==2 & zone==4
replace field_size = area_size*0.0023 if area_unit==2 & zone==5
replace field_size = area_size*0.00001 if area_unit==2 & zone==6

replace field_size = area_size*0.00006 if area_unit==3 & zone==1						//reported in stands
replace field_size = area_size*0.00016 if area_unit==3 & zone==2
replace field_size = area_size*0.00004 if area_unit==3 & zone==3
replace field_size = area_size*0.00004 if area_unit==3 & zone==4
replace field_size = area_size*0.00013 if area_unit==3 & zone==5
replace field_size = area_size*0.00041 if area_unit==3 & zone==6



/*ALT 02.23.23*/ gen area_est = field_size
*replacing farmer reported with GPS if available
replace field_size = area_meas_sqm*0.0001 if area_meas_sqm!=.               				
gen gps_meas = (area_meas_sqm!=.)
la var gps_meas "Plot was measured with GPS, 1=Yes"
ren plotid plot_id
*Total land holding including cultivated and rented out
collapse (sum) field_size, by (hhid)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/weight.dta", gen(wgt)

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

foreach v of varlist  field_size  {
	_pctile `v' [aw=weight] , p(5 99) 
	gen `v'_w=`v'
	replace  `v'_w = r(r1) if  `v'_w < r(r1) &  `v'_w!=.
	replace  `v'_w = r(r2) if  `v'_w > r(r2) &  `v'_w!=.
	local l`v' : var lab `v'
	lab var  `v'_w  "`l`v'' - Winzorized top 5% & bottom 1%"
}

tab field_size
tab field_size_w, missing
sum field_size field_size_w, detail



sort hhid
ren field_size_w land_holding
keep hhid land_holding
label var land_holding "land holding in hectares"
save "${Nigeria_GHS_W4_created_data}\land_holding_2018.dta", replace







*******************************
*Soil Quality
*******************************

//ALT IMPORTANT NOTE: As of W4, the implied area conversions for farmer estimated units (including hectares) are markedly different from previous waves. I recommend excluding plots that do not have GPS measured areas from any area-based productivity estimates.
use "${Nigeria_GHS_W4_raw_data}/sect11a1_plantingw4.dta", clear
*merging in planting section to get cultivated status
merge 1:1 hhid plotid using "${Nigeria_GHS_W4_raw_data}/sect11b1_plantingw4.dta", nogen
*merging in harvest section to get areas for new plots
merge 1:1 hhid plotid using "${Nigeria_GHS_W4_raw_data}/secta1_harvestw4.dta", gen(plot_merge)
 
ren s11aq4aa area_size
ren s11aq4b area_unit
ren sa1q11 area_size2 //GPS measurement, no units in file
//ren sa1q9b area_unit2 //Not in file
ren s11aq4c area_meas_sqm
//ren sa1q9c area_meas_sqm2
gen cultivate = s11b1q27 ==1 


gen field_size= area_size if area_unit==6
replace field_size = area_size*0.0667 if area_unit==4									//reported in plots
replace field_size = area_size*0.404686 if area_unit==5		    						//reported in acres
replace field_size = area_size*0.0001 if area_unit==7									//reported in square meters

replace field_size = area_size*0.00012 if area_unit==1 & zone==1						//reported in heaps
replace field_size = area_size*0.00016 if area_unit==1 & zone==2
replace field_size = area_size*0.00011 if area_unit==1 & zone==3
replace field_size = area_size*0.00019 if area_unit==1 & zone==4
replace field_size = area_size*0.00021 if area_unit==1 & zone==5
replace field_size = area_size*0.00012 if area_unit==1 & zone==6

replace field_size = area_size*0.0027 if area_unit==2 & zone==1							//reported in ridges
replace field_size = area_size*0.004 if area_unit==2 & zone==2
replace field_size = area_size*0.00494 if area_unit==2 & zone==3
replace field_size = area_size*0.0023 if area_unit==2 & zone==4
replace field_size = area_size*0.0023 if area_unit==2 & zone==5
replace field_size = area_size*0.00001 if area_unit==2 & zone==6

replace field_size = area_size*0.00006 if area_unit==3 & zone==1						//reported in stands
replace field_size = area_size*0.00016 if area_unit==3 & zone==2
replace field_size = area_size*0.00004 if area_unit==3 & zone==3
replace field_size = area_size*0.00004 if area_unit==3 & zone==4
replace field_size = area_size*0.00013 if area_unit==3 & zone==5
replace field_size = area_size*0.00041 if area_unit==3 & zone==6



/*ALT 02.23.23*/ gen area_est = field_size
*replacing farmer reported with GPS if available
replace field_size = area_meas_sqm*0.0001 if area_meas_sqm!=.               				
gen gps_meas = (area_meas_sqm!=.)
la var gps_meas "Plot was measured with GPS, 1=Yes"

keep zone state lga sector ea hhid plotid field_size

merge 1:1 hhid plotid using "${Nigeria_GHS_W4_raw_data}\sect11b1_plantingw4.dta"
merge m:1 hhid using "${Nigeria_GHS_W4_created_data}/ag_rainy_18.dta", gen(filter)

keep if ag_rainy_18==1

ren s11b1q45 soil_quality
tab soil_quality, missing
order field_size soil_quality hhid 
sort hhid


egen max_fieldsize = max(field_size), by (hhid)
replace max_fieldsize= . if max_fieldsize!= max_fieldsize
order field_size soil_quality hhid max_fieldsize
sort hhid
keep if field_size== max_fieldsize
sort hhid plotid field_size

duplicates report hhid

duplicates tag hhid, generate(dup)
tab dup
list field_size soil_quality dup


list hhid plotid field_size soil_quality dup if dup>0

egen soil_qty_rev = min(soil_quality) 
gen soil_qty_rev2 = soil_quality

replace soil_qty_rev2 = soil_qty_rev if dup>0

list hhid plotid  field_size soil_quality soil_qty_rev soil_qty_rev2 dup if dup>0



egen med_soil_ea = median(soil_qty_rev2), by (ea)
egen med_soil_lga = median(soil_qty_rev2), by (lga)
egen med_soil_state = median(soil_qty_rev2), by (state)
egen med_soil_zone = median(soil_qty_rev2), by (zone)

replace soil_qty_rev2= med_soil_ea if soil_qty_rev2==.
tab soil_qty_rev2, missing
replace soil_qty_rev2= med_soil_lga if soil_qty_rev2==.
tab soil_qty_rev2, missing
replace soil_qty_rev2= med_soil_state if soil_qty_rev2==.
tab soil_qty_rev2, missing
replace soil_qty_rev2= med_soil_zone if soil_qty_rev2==.
tab soil_qty_rev2, missing

replace soil_qty_rev2= 2 if soil_qty_rev2==1.5
tab soil_qty_rev2, missing

la define soil 1 "Good" 2 "fair" 3 "poor"

*la value soil soil_qty_rev2

collapse (mean) soil_qty_rev2 , by (hhid)
la var soil_qty_rev2 "1=Good 2= fair 3=Bad "
save "${Nigeria_GHS_W4_created_data}\soil_quality_2018.dta", replace





























************************* Merging Agricultural Datasets ********************

use "${Nigeria_GHS_W4_created_data}\purchased_fert_2018.dta", replace


*******All observations Merged*****

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\subsidized_fert_2018.dta"
drop _merge
sort hhid


merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\total_qty_2018.dta"
drop _merge
sort hhid

merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\weight.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\savings_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\credit_access_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\extension_access_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\demographics_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\laborage_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\safety_net_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\food_prices_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\soil_quality_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\geodata_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\household_asset_2018.dta"
drop _merge
sort hhid
merge 1:1 hhid using "${Nigeria_GHS_W4_created_data}\land_holding_2018.dta"
drop _merge

gen year = 2018
sort hhid
save "${Nigeria_GHS_W4_created_data}/Nigeria_wave4_completedata_2018.dta", replace


tabstat total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding [aweight = weight], statistics( mean median sd min max ) columns(statistics)








misstable summarize subsidy_dummy femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer 
proportion subsidy_dummy femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2















*****************Appending all Nigeria Datasets*****************
use "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\nga_wave2010\Nigeria_wave1_complete_data.dta",clear
append using "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\nga_wave2012\Nigeria_wave2_complete_data.dta.dta"
append using "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\nga_wave2015\Nigeria_wave3_completedata_2015.dta"
append using "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\nga_wave2018\Nigeria_wave4_completedata_2018.dta"

order year



save "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\complete_files", replace







tabstat total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding [aweight = weight], statistics( mean median sd min max ) columns(statistics)








misstable summarize subsidy_dummy femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer 
proportion subsidy_dummy femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2







