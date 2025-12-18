********************************************************************************
* 论文案例1: Taxation and Innovation in the Twentieth Century
* Akcigit, Grigsby, Nicholas, and Stantcheva（2022，QJE）
********************************************************************************

*----- 图像设置 -----*
clear all
set scheme s2color

*----- 基础设置 -----*
global main "/Users/lmm/Documents/Binscatter论文复现/"
global data "$main/02_data_raw"
global temp "$main/03_temp"
global output "$main/04_output"

cd "$main"

use "$data/CCFF_2024_AER--AGNS.dta", clear

********************************************************************************
* 方法比对
********************************************************************************

*-------------------- A1. 残差化因变量 lnpat --------------------*
areg lnpat top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 i.statenum ///
    [aw=pop1940], a(year)
predict resid_lnpat, res
egen mean_lnpat = mean(lnpat)
replace resid_lnpat = resid_lnpat + mean_lnpat  // "均值校准"：将残差平移回原始量级附近

*-------------------- A2. 残差化核心解释变量 mtr90_lag3 --------------------*
areg mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 i.statenum ///
    [aw=pop1940], a(year)
predict resid_mtr90_lag3, res
egen mean_mtr90_lag3 = mean(mtr90_lag3)
replace resid_mtr90_lag3 = resid_mtr90_lag3 + mean_mtr90_lag3   // 同样做均值校准，使横轴变量回到原量级附近

*-------------------- A3. 残差化 binscatter：叠加线性对照线 --------------------*
binsreg resid_lnpat resid_mtr90_lag3 [aw=pop1940], ///
    nbins(50) polyreg(1) ///
    dotsplotopt(mcolor(dkorange)) ///
    polyregplotopt(lcolor(black)) ///
    savedata($temp/tmpAGNS1binscatter) replace ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(5.7(0.1)6.2, nogrid) xlabel(-0.525(0.05)-0.325)

graph export "$output/AGNS_binscatter.png", as(png) replace width(800) height(600)


*-------------------- B1. 生成协变量调整 binsreg 的 savedata --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], ///
    nbins(50) absorb(statenum year) at(mean) polyreg(1) ///
    savedata($temp/tmpAGNS1binsreg) replace
preserve

*-------------------- B2. 合并两套 savedata，并标记来源 --------------------*
use "$temp/tmpAGNS1binscatter", clear
gen binsreg = 0
append using "$temp/tmpAGNS1binsreg"
replace binsreg = 1 if missing(binsreg)  // 1：协变量调整 binsreg 数据（append 进来时 binsreg 变量缺失）

*-------------------- B3. 人工增加两条观测，用于正确的表示原来的比例尺（虚线） --------------------*
local new = _N + 2
set obs `new'    // 扩充两行观测作为线段的端点

replace poly_x = -0.8 if _n == `new'-1
replace poly_x =  0   if _n == `new'   // 使横轴范围与协变量调整图的横轴范围一致

replace binsreg = 2 if _n >= `new'-1   
gen poly_fit2 = .

reg poly_fit poly_x if binsreg==0 // 用残差化方法（binsreg==0）的 poly_fit 与 poly_x 拟合线性关系

replace poly_fit2 = _b[_cons] + _b[poly_x]*poly_x if binsreg==2

*-------------------- B4. 统一绘图：点（残差化）+ 实线（残差化线性）+ 虚线延伸 --------------------*
tw ///
    (scatter dots_fit dots_x if binsreg==0, color(dkorange)) ///
    (line poly_fit poly_x if binsreg==0, color(black) lwidth(medthick)) ///
    (line poly_fit2 poly_x if binsreg==2, color(black) lpattern(dot) lwidth(medthick)), ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(5(.5)8) legend(off)

graph export "$output/AGNS_covariateAdjustments_binscatter.png", replace as(png) width(800) height(600)

restore


*-------------------- C1. 协变量调整 binsreg（不含线性对照线） --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], ///
    nbins(50) ///
    absorb(statenum year) replace at(mean) ///
    dotsplotopt(mcolor(black)) ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0)
    // 移除 polyreg(1)，仅展示协变量调整后的分箱点/估计形态

graph export "$output/AGNS_covariateAdjustments_binsreg_woutRegLine.png", replace as(png) width(800) height(600)


********************************************************************************
* 分箱数的选择
********************************************************************************
use $data/CCFF_2024_AER--AGNS, clear

* 5 bins -- too few for consistent estimation, interpret as deciles
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], absorb(statenum year) nbins(5) at(mean) polyreg(1) dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) plotregion(lcolor(black)) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9) 
graph export "$output/AGNS_nbins5.png", replace as(png) width(800) height(600)

* IMSE optimal
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], absorb(statenum year) vce(cluster fiveyrblockbystate year) randcut(1) at(mean) polyreg(1) dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) plotregion(lcolor(black)) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9) 
graph export "$output/AGNS_nbinsOptimal.png", replace as(png) width(800) height(600)

* 50 bins -- undersmoothed
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], absorb(statenum year) nbins(50) at(mean)  polyreg(1) dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) plotregion(lcolor(black)) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9) 
*xlabel(-0.9(0.3)0) ylabel(5(.5)8)
graph export "$output/AGNS_nbins50.png", replace as(png) width(800) height(600)


********************************************************************************
* 箱数对置信带的影响
********************************************************************************

binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], randcut(1) absorb(statenum year) nbins(5) at(mean) cb(0 0) dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) plotregion(lcolor(black)) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/AGNS_fixedJconfidenceBandwDots.png", replace as(png) width(800) height(600)

use "$data/CCFF_2024_AER--AGNS.dta", clear
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], vce(cluster fiveyrblockbystate year) absorb(statenum year) randcut(1) cb(T) replace at(mean) dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) plotregion(lcolor(black)) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/AGNS_confidenceBandwDots.png", as(png) replace 

********************************************************************************
* 置信带的推断
********************************************************************************
use $data/CCFF_2024_AER--AGNS, clear

binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 [aw=pop1940], vce(cluster fiveyrblockbystate year) absorb(statenum year) randcut(1) cb(T) at(mean) nodraw savedata($temp/tmp3) replace nsims(5000) simsgrid(200)

use $temp/tmp3, clear

qui summ CB_r if CB_x < -0.8
gen horline = `r(max)'

gen CB_avg = 0.4*CB_r+0.6*CB_l
reg CB_avg CB_x
predict inline, xb 

tw (rarea CB_l CB_r CB_x, fcolor(edkblue) fintensity(20) lwidth(none)) (line inline CB_x, lcolor(red)), graphregion(color(white) margin(large)) plotregion(lcolor(black)) ylabel(, nogrid nogextend) legend(off) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率")
graph export "$output/AGNS_confidenceBandwLine.png", replace as(png) width(800) height(600)

tw (rarea CB_l CB_r CB_x, fcolor(edkblue) fintensity(20) lwidth(none)) (scatter dots_fit dots_x, mcolor(black)) (line horline CB_x, lcolor(black) lpattern(dash)), graphregion(color(white) margin(large)) plotregion(lcolor(black)) ylabel(, nogrid nogextend) legend(off) ytitle("专利数量") xtitle("前 90% 收入者的边际净税率")
graph export "$output/AGNS_confidenceBandwHorLine.png", replace as(png) width(800) height(600)






