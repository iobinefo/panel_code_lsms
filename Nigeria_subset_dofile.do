
use "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\complete_files.dta", clear


gen dummy = 1

collapse (sum) dummy, by (hhid)
keep if dummy==4
sort hhid

save "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\subset_complete_files", replace


merge 1:m hhid using "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\complete_files.dta"

drop if _merge==2




gen year_2010 = (year==2010)
gen year_2012 = (year==2012)
gen year_2015 = (year==2015)
gen year_2018 = (year==2018)

gen commercial_dummy = (total_qty_w>0)





local time_avg "total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  subsidy_dummy femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2"

foreach x in `time_avg' {

	bysort hhid : egen TAvg_`x' = mean(`x')

}

** OLS **
reg total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 i.zone i.year


** OLS with HH fixed effects
xtreg total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 i.zone i.year, fe i(hhid) cluster(hhid)


** CRE-TOBIT 
tobit total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 i.zone i.year, ll(0)

margins, predict(ystar(0,.)) dydx(*) post




** Double Hurdle **

capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	summarize bsdEy_dxj
	
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot



****************************
*Subsidy Analysis
****************************

	
//% of HHs that bought commercial fertilizer by each survey wave
bysort year : tabstat subsidy_dummy [w=weight], stat(mean sem)

// By HH, sum the binary variable of commerical fert market particapation for all waves
bysort hhid : egen sum_4waves_com_fer_bin = sum(subsidy_dummy) 



****************************
*Commercial Analysis
****************************

	
//% of HHs that bought commercial fertilizer by each survey wave
bysort year : tabstat commercial_dummy [w=weight], stat(mean sem)

// By HH, sum the binary variable of commerical fert market particapation for all waves
bysort hhid : egen sum_4waves_com_fer_bin = sum(commercial_dummy) 






keep if year==2010

gen group = 1 if sum_4waves_com_fer_bin == 4 // HHs that purchased comm fert in all waves

replace group= 2 if sum_4waves_com_fer_bin>0 & sum_4waves_com_fer_bin<4 //HHs that purchased comm fert in some waves

replace group = 3 if sum_4waves_com_fer_bin==0 // HHs that never purchased comm fert in all waves

//% of HHs in each comm fert particiaption group 
tab group

global outreg "C:\Users\obine\Music\Documents\Smallholder lsms STATA\analyzed_data\try_ng\"


* Difference between always vs (sometimes and never)
preserve
gen always1others0 = 1 if group==1 
replace always1others0 = 0 if always1others0==.

reg always1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone

reg always1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(ea)

reg always1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(lga)
outreg2 using "$outreg/particiaption_group.doc", title (Table 3: Factors affecting if households bought commercial fertilizer in every survey wave vs others) ctitle(OLS) se dec(3) label replace

restore
            





* Difference between never vs (always and sometimes)
preserve

gen never1others0 = 1 if group == 3
replace never1others0 = 0 if never1others0==.

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(ea)

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(lga)
outreg2 using "$outreg/particiaption_group.doc", title (Table 4: Factors affecting if household never bought commercial fertilizer in any survey wave vs others) ctitle(OLS) se dec(3) label replace

restore








































                                   *********************************************** 
								   *Crowding out estimate for each asset quintile
								   ***********************************************
	
	xtile asset_quintiles=real_hhvalue, nq(5)

	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	
	summarize bsdEy_dxj // Average crowding out estimate for reference. Use SE from bootstrap 

	tabulate asset_quintiles, summarize (bsdEy_dxj) // Crowding out estimate for each quintile
	restore


** Bootstraped SE for each quintiles

*Asset q1
capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	keep if asset_quintiles==1
	
	summarize bsdEy_dxj
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est_asset_q1 = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot

*Asset q2
capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	keep if asset_quintiles==2
	
	summarize bsdEy_dxj
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est_asset_q2 = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot

*Asset q3
capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	keep if asset_quintiles==3
	
	summarize bsdEy_dxj
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est_asset_q3 = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot


*Asset q4
capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	keep if asset_quintiles==4
	
	summarize bsdEy_dxj
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est_asset_q4 = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot

*Asset q5
capture program drop APEboot
program define APEboot, rclass
	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	keep if asset_quintiles==5
	
	summarize bsdEy_dxj
	return scalar ape_xj =r(mean)
	matrix ape_xj=r(ape_xj)
	restore
end
bootstrap crowd_out_est_asset_q5 = r(ape_xj), reps(250) cluster(hhid) idcluster(newid): APEboot


//////////////////////////////Mean asset by quantiles////////////////////////
sum real_hhvalue if asset_quintiles==1, detail
sum real_hhvalue if asset_quintiles==2, detail
sum real_hhvalue if asset_quintiles==3, detail
sum real_hhvalue if asset_quintiles==4, detail
sum real_hhvalue if asset_quintiles==5, detail



sum real_hhvalue[w=weight] if asset_quintiles==1
sum real_hhvalue[w=weight] if asset_quintiles==2
sum real_hhvalue[w=weight] if asset_quintiles==3
sum real_hhvalue[w=weight] if asset_quintiles==4
sum real_hhvalue[w=weight] if asset_quintiles==5





















                                   *********************************************** 
								   *Crowding out estimate for each year
								   ***********************************************


	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	
	summarize bsdEy_dxj // Average crowding out estimate for reference. Use SE from bootstrap 

	tabulate year, summarize (bsdEy_dxj) // Crowding out estimate for each quintile
	restore

	
	
	



                                   *********************************************** 
								   *Crowding out estimate for land_holding
								   ***********************************************
	
	xtile land_quintiles=land_holding, nq(5)

	preserve
	craggit commercial_dummy subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding  femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018, second(total_qty_w subsidy_qty_w mrk_dist_w real_tpricefert_cens_mrk num_mem hh_headage real_hhvalue worker real_maize_price_mr real_rice_price_mr land_holding femhead informal_save formal_credit informal_credit ext_acess attend_sch pry_edu finish_pry finish_sec safety_net net_seller net_buyer soil_qty_rev2 TAvg_total_qty_w TAvg_subsidy_qty_w TAvg_mrk_dist_w TAvg_real_tpricefert_cens_mrk TAvg_num_mem TAvg_hh_headage TAvg_real_hhvalue TAvg_worker TAvg_real_maize_price_mr TAvg_real_rice_price_mr TAvg_land_holding TAvg_femhead TAvg_informal_save TAvg_formal_credit TAvg_informal_credit TAvg_ext_acess TAvg_attend_sch TAvg_pry_edu TAvg_finish_pry TAvg_finish_sec TAvg_safety_net TAvg_net_seller TAvg_net_buyer TAvg_soil_qty_rev2 year_2010 year_2012 year_2015 year_2018) cluster(hhid)

	predict bsx1g, eq(Tier1)
	predict bsx2b, eq(Tier2)
	predict  bssigma , eq(sigma)
	generate bsIMR = normalden(bsx2b/bssigma)/normal(bsx2b/bssigma)
	generate bsdEy_dxj = 												///
				[Tier1]_b[subsidy_qty_w]*normalden(bsx1g)*(bsx2b+bssigma*bsIMR) ///
				+[Tier2]_b[subsidy_qty_w]*normal(bsx1g)*(1-bsIMR*(bsx2b/bssigma+bsIMR))
	
	
	summarize bsdEy_dxj // Average crowding out estimate for reference. Use SE from bootstrap 

	tabulate land_quintiles, summarize (bsdEy_dxj) // Crowding out estimate for each quintile
	restore	

	
//////////////////////////////Mean land_holding by quantiles////////////////////////
sum land_holding[w=weight] if land_quintiles==1
sum land_holding[w=weight] if land_quintiles==2
sum land_holding[w=weight] if land_quintiles==3
sum land_holding[w=weight] if land_quintiles==4
sum land_holding[w=weight] if land_quintiles==5



///////////Average commercial fertilizer bought by land_quintiles//////////////
sum total_qty_w if land_quintiles==1
sum total_qty_w if land_quintiles==2
sum total_qty_w if land_quintiles==3
sum total_qty_w if land_quintiles==4
sum total_qty_w if land_quintiles==5

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////Mean asset////////////////////

*************North Central South
gen north_land = land_holding if zone== 2 | zone==3
gen north_fert = total_qty_w if zone== 2 | zone==3

gen central_land = land_holding if zone==1
gen central_fert = total_qty_w if zone== 1

gen south_land = land_holding if zone== 4 | zone==5 | zone==6
gen south_fert = total_qty_w if zone== 4 | zone==5 | zone==6


sum north_land [w=weight] if land_holding <2
sum north_land [w=weight] if land_holding >=2 & land_holding <=5
sum north_land [w=weight] if  land_holding >5

sum north_fert [w=weight] if land_holding <2
sum north_fert [w=weight] if land_holding >=2 & land_holding <=5
sum north_fert [w=weight] if  land_holding >5



sum central_land [w=weight] if land_holding <2
sum central_land [w=weight] if land_holding >=2 & land_holding <=5
sum central_land [w=weight] if  land_holding >5

sum central_fert [w=weight] if land_holding <2
sum central_fert [w=weight] if land_holding >=2 & land_holding <=5
sum central_fert [w=weight] if  land_holding >5

sum south_land [w=weight] if land_holding <2
sum south_land [w=weight] if land_holding >=2 & land_holding <=5
sum south_land [w=weight] if  land_holding >5

sum south_fert [w=weight] if land_holding <2
sum south_fert [w=weight] if land_holding >=2 & land_holding <=5
sum south_fert [w=weight] if  land_holding >5








*************Urban Rural
gen rural_land = land_holding if sector== 2 
gen rural_fert = total_qty_w if sector== 2

gen urban_land = land_holding if sector==1
gen urban_fert = total_qty_w if sector== 1


sum  rural_land [w=weight]  if land_holding <2
sum rural_land [w=weight] if land_holding >=2 & land_holding <=5
sum rural_land [w=weight] if  land_holding >5

sum rural_fert [w=weight] if land_holding <2
sum rural_fert [w=weight] if land_holding >=2 & land_holding <=5
sum rural_fert [w=weight] if  land_holding >5




sum urban_land [w=weight] if land_holding <2
sum urban_land [w=weight] if land_holding >=2 & land_holding <=5
sum urban_land [w=weight] if  land_holding >5

sum urban_fert [w=weight] if land_holding <2
sum urban_fert [w=weight] if land_holding >=2 & land_holding <=5
sum urban_fert [w=weight] if  land_holding >5	

*********************National- Total

sum land_holding [w=weight] if land_holding <2
sum land_holding [w=weight] if land_holding >=2 & land_holding <=5
sum land_holding [w=weight] if  land_holding >5

sum total_qty_w [w=weight] if land_holding <2
sum total_qty_w [w=weight] if land_holding >=2 & land_holding <=5
sum total_qty_w [w=weight] if  land_holding >5



//////% of sample that used fertilizer in the group


gen use_fert = (total_qty_w>0)
tab use_fert
/*
tab use_fert if zone== 2 | zone==3 & land_holding <2
tab use_fert if zone== 2 | zone==3 & land_holding >=2 & land_holding <=5
tab use_fert if zone== 2 | zone==3 & land_holding >5

tab use_fert if  zone== 4 | zone==5 | zone==6 & land_holding <2
tab use_fert if  zone== 4 | zone==5 | zone==6 & land_holding >=2 & land_holding <=5
tab use_fert if  zone== 4 | zone==5 | zone==6 & land_holding >5

*/





gen land_group = 1 if land_holding <2
replace land_group =2 if land_holding >=2 & land_holding <=5
replace land_group =3 if land_holding >5



gen region = 1 if zone== 2 | zone==3
replace region =2 if zone== 1
replace region =3 if zone== 4 | zone==5 | zone==6

tab land_group  region, column


gen rural_used =1 if sector== 2 
replace rural_used =2 if sector ==1

tab land_group rural_used, column


tab land_group use_fert, column















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////checking for 2 or more waves////////////////



	
//% of HHs that bought commercial fertilizer by each survey wave
bysort year : tabstat commercial_dummy [w=weight], stat(mean sem)

// By HH, sum the binary variable of commerical fert market particapation for all waves
bysort hhid : egen sum_4waves_com_fer_bin = sum(commercial_dummy) 






keep if year==2010

gen group = 1 if sum_4waves_com_fer_bin == 4 // HHs that purchased comm fert in all waves

*replace group= 2 if sum_4waves_com_fer_bin>0 & sum_4waves_com_fer_bin<4 //HHs that purchased comm fert in some waves

replace group = 3 if sum_4waves_com_fer_bin==0 // HHs that never purchased comm fert in all waves

replace group = 4 if sum_4waves_com_fer_bin==3  // HHs that purchased comm fert in 3 waves

replace group = 5 if sum_4waves_com_fer_bin ==2   //HHs that purchased comm fert in 2 waves



//% of HHs in each comm fert particiaption group 
tab group
	
*br total_qty_w  subsidy_qty_w real_hhvalue land_holding real_tpricefert_cens_mrk num_mem femhead attend_sch if group ==3
	
	
	
	
	


gen use_two = 1 if group ==5 | group ==4 | group ==1
*gen use_two = 1 if group ==2
replace use_two = 0 if use_two ==.
tab use_two
	
	

* Difference between always vs (sometimes and never)
preserve

reg use_two subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone

reg use_two subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(ea)

reg use_two subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(lga)
outreg2 using "$outreg/particiaption_group.doc", title (Table 3: Factors affecting if households bought commercial fertilizer in every survey wave vs others) ctitle(OLS) se dec(3) label replace

restore


* Difference between never vs (always and sometimes)
preserve

gen never1others0 = 1 if group == 3
replace never1others0 = 0 if never1others0==.

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(ea)

reg never1others0 subsidy_dummy real_tpricefert_cens_mrk real_maize_price_mr real_rice_price_mr land_holding num_mem worker hh_headage attend_sch pry_edu finish_pry finish_sec femhead informal_save formal_credit informal_credit net_seller net_buyer mrk_dist_w ext_acess safety_net soil_qty_rev2 real_hhvalue i.zone, cluster(lga)
outreg2 using "$outreg/particiaption_group.doc", title (Table 4: Factors affecting if household never bought commercial fertilizer in any survey wave vs others) ctitle(OLS) se dec(3) label replace

restore

