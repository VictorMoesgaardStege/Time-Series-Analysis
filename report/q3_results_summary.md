# Question 3 Results Summary

## 3.1 Time Series
**Figure:** `report/figures/q3_31_time_series.pdf`
**Ph range:** 14 – 98 W
**Tdelta range:** 12.6 – 25.4 °C
**Gv range:** -2.2 – 876 W/m²
**N total:** 231 hourly observations

## 3.2 Train/Test Split
**Training set:** 167 obs (rows 1–167, thour up to 185), up to 2013-02-06
**Test set:** 64 obs, thour 168–231, from 2013-02-06 01:00:00

## 3.3 Exploratory Analysis
**Figure A (scatter):** `report/figures/q3_33_scatter.pdf`
**Figure B (ACF/PACF/CCF):** `report/figures/q3_33_acf_ccf.pdf`
**ACF Ph:** slow exponential decay — strong AR dynamics needed
**PACF Ph:** cuts off sharply after lag 1–2 — suggests low-order AR
**Raw CCF Tdelta→Ph:** positive at lag 0, persists several lags
**Raw CCF Gv→Ph:** strongly negative at lag 0 (solar gain replaces heating)
**Note:** raw CCF is biased when input is autocorrelated — see prewhitened version in 3.4

## 3.4 Impulse Response (Prewhitening)
**Figure:** `report/figures/q3_34_impulse_response.pdf`
**AR order used to prewhiten Tdelta:** 2
**AR order used to prewhiten Gv:** 1
**Tdelta IR peak lag:** 0 h, CCF = 0.5
**Gv IR peak lag:** 0 h, CCF = 0.837

## 3.5 Linear Regression (no AR)
**Figure:** `report/figures/q3_35_lm_diagnostics.pdf`
**Table:** `report/tables/q3_35_lm_coef.tex`
**R²:** 0.8981 | **Adj R²:** 0.8968
**Residual std error:** 5.182 W
**omega1 (Tdelta):** 3.3357 (p=1.91e-57)
**omega2 (Gv):** -0.1114 (p=1.05e-64)
**Residual ACF:** significant autocorrelation remains → AR terms needed
**Residual CCF:** structure visible vs both inputs → transfer function model needed

## 3.6 ARX(1)
**Figure:** `report/figures/q3_36_arx1_diagnostics.pdf`
**Table:** `report/tables/q3_36_arx1_coef.tex`
**R²:** 0.9749 | **Adj R²:** 0.9744
**Residual std error:** 2.581 W
**phi1 (Ph.l1):** 0.4046 (p=2.02e-51)
**omega1 (Tdelta):** 2.096 (p=1.31e-55)
**omega2 (Gv):** -0.0851 (p=2.3e-81)
**Improvement over OLS:** residual autocorrelation substantially reduced

## 3.7 Model Order Selection (AIC / BIC)
**Figure:** `report/figures/q3_37_aic_bic.pdf`
**Table:** `report/tables/q3_37_ic.tex`
**Best order by AIC:** 6 (AIC = 688.73)
**Best order by BIC:** 6 (BIC = 751.09)
**AIC values (p=1..10):** 796.6, 750.7, 734.7, 709.9, 702.4, 688.7, 690.3, 695.1, 699.4, 698.9
**BIC values (p=1..10):** 812.2, 775.6, 769, 753.6, 755.4, 751.1, 762, 776.2, 789.9, 798.6

## 3.8 Test-Set RMSE vs. Model Order
**Figure:** `report/figures/q3_38_rmse.pdf`
**Table:** `report/tables/q3_38_rmse.tex`
**Best order by RMSE:** 3 (RMSE = 2.8304 W)
**RMSE values (p=1..10):** 3.663, 3.195, 2.83, 2.958, 3.046, 3.314, 3.283, 3.274, 3.282, 3.301

## 3.9 Multi-Step Simulation
**Figure:** `report/figures/q3_39_multistep_simulation.pdf`
**Table (coefficients):** `report/tables/q3_39_arx6_coef.tex`
**Model used:** ARX(6)
**Simulation RMSE — training period:** 1.769 W
**Simulation RMSE — test period:** 3.404 W
**Coefficients:**
  - (Intercept): 4.9112
  - Ph.l1: 0.2599
  - Ph.l2: -0.0922
  - Ph.l3: -0.0551
  - Ph.l4: 0.0384
  - Ph.l5: 0.0786
  - Ph.l6: 0.0115
  - Tdelta.l0: 0.1061
  - Tdelta.l1: 1.6519
  - Tdelta.l2: -0.7116
  - Tdelta.l3: 2.6783
  - Tdelta.l4: -2.4155
  - Tdelta.l5: 1.5028
  - Gv.l0: -0.0997
  - Gv.l1: -0.0013
  - Gv.l2: -0.0141
  - Gv.l3: 0.0024
  - Gv.l4: -0.0124
  - Gv.l5: 0.0078

## 3.10 Summary
| Criterion | Best order |
|-----------|------------|
| AIC       | 6 |
| BIC       | 6 |
| Test RMSE | 3 |
| Chosen    | 6 (BIC) |

