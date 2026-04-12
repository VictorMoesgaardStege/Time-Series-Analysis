# Question 3 Results Summary

## 3.1 Time Series
**Figure:** `report/figures/q3_31_time_series.pdf`
**Ph range:** 14 – 98 W
**Tdelta range:** 12.6 – 25.4 °C
**Gv range:** -2.2 – 876 W/m²
**N total:** 231 hourly observations

## 3.2 Train/Test Split
**Training set:** 149 obs, thour 1–167, up to 2013-02-05 06:00:00
**Test set:** 82 obs, thour 168–231, from 2013-02-05 07:00:00

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
**AR order used to prewhiten Gv:** 2
**Tdelta IR peak lag:** 0 h, CCF = 0.517
**Gv IR peak lag:** 0 h, CCF = 0.823

## 3.5 Linear Regression (no AR)
**Figure:** `report/figures/q3_35_lm_diagnostics.pdf`
**Table:** `report/tables/q3_35_lm_coef.tex`
**R²:** 0.9029 | **Adj R²:** 0.9015
**Residual std error:** 4.913 W
**omega1 (Tdelta):** 3.4212 (p=5.63e-54)
**omega2 (Gv):** -0.1131 (p=1.69e-57)
**Residual ACF:** significant autocorrelation remains → AR terms needed
**Residual CCF:** structure visible vs both inputs → transfer function model needed

## 3.6 ARX(1)
**Figure:** `report/figures/q3_36_arx1_diagnostics.pdf`
**Table:** `report/tables/q3_36_arx1_coef.tex`
**R²:** 0.9766 | **Adj R²:** 0.9761
**Residual std error:** 2.419 W
**phi1 (Ph.l1):** 0.4149 (p=1.14e-46)
**omega1 (Tdelta):** 2.0841 (p=3.5e-49)
**omega2 (Gv):** -0.0845 (p=1.21e-70)
**Improvement over OLS:** residual autocorrelation substantially reduced

## 3.7 Model Order Selection (AIC / BIC)
**Figure:** `report/figures/q3_37_aic_bic.pdf`
**Table:** `report/tables/q3_37_ic.tex`
**Best order by AIC:** 5 (AIC = 585.83)
**Best order by BIC:** 3 (BIC = 638.57)
**AIC values (p=1..10):** 643.5, 609.3, 599.5, 594.2, 585.8, 591.6, 594.8, 598.9, 601.4, 605.6
**BIC values (p=1..10):** 664.5, 639.3, 638.6, 642.3, 642.9, 657.6, 669.9, 683, 694.5, 707.8

## 3.8 Test-Set RMSE vs. Model Order
**Figure:** `report/figures/q3_38_rmse.pdf`
**Table:** `report/tables/q3_38_rmse.tex`
**Best order by RMSE:** 4 (RMSE = 3.0341 W)
**RMSE values (p=1..10):** 3.278, 3.057, 3.053, 3.034, 3.165, 3.161, 3.162, 3.177, 3.198, 3.21

## 3.9 Multi-Step Simulation
**Figure:** `report/figures/q3_39_multistep_simulation.pdf`
**Table (coefficients):** `report/tables/q3_39_arx3_coef.tex`
**Model used:** ARX(3)
**Simulation RMSE — training period:** 1.652 W
**Simulation RMSE — test period:** 3.208 W
**Coefficients:**
  - (Intercept): 4.1155
  - Ph.l1: 0.1117
  - Ph.l2: 0.1306
  - Ph.l3: 0.164
  - Tdelta.l0: -0.1185
  - Tdelta.l1: 2.3108
  - Tdelta.l2: -0.8871
  - Tdelta.l3: 0.8741
  - Gv.l0: -0.0973
  - Gv.l1: -0.0164
  - Gv.l2: 0.0072
  - Gv.l3: 0.0154

## 3.10 Summary
| Criterion | Best order |
|-----------|------------|
| AIC       | 5 |
| BIC       | 3 |
| Test RMSE | 4 |
| Chosen    | 3 (BIC) |

