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
**Note:** raw CCF is only exploratory here; in 3.4 the impulse responses are estimated jointly by FIR regression rather than pre-whitened CCF

## 3.4 Impulse Response (Joint FIR Regression)
**Figure:** `report/figures/q3_34_impulse_response.png`
**Table:** `report/tables/q3_34_fir_coef.tex`
**CSV output:** `output/tables/q3_34_impulse_response.csv`
**Method:** joint FIR regression with lags 0–10 for both Tdelta and Gv
**Largest |Tdelta impulse-response coefficient| at lag:** 3 h
**Largest |Gv impulse-response coefficient| at lag:** 0 h
**Comment:** pre-whitened CCF is a standard textbook identification idea for single-input transfer-function models, but with two potentially correlated inputs it is less reliable. Therefore the impulse responses are estimated jointly here through finite distributed lags.

## 3.5 Linear Regression (no AR)
**Figure:** `report/figures/q3_35_lm_diagnostics.pdf`
**Table:** `report/tables/q3_35_lm_coef.tex`
**R²:** 0.9945 | **Adj R²:** 0.9944
**Residual std error:** 5.44 W
**omega1 (Tdelta):** 3.8948 (p=2.79e-186)
**omega2 (Gv):** -0.1099 (p=1.66e-61)
**Residual ACF:** significant autocorrelation remains → AR terms needed
**Residual CCF:** structure visible vs both inputs → transfer function model needed

## 3.6 ARX(1)
**Figure:** `report/figures/q3_36_arx1_diagnostics.pdf`
**Table:** `report/tables/q3_36_arx1_coef.tex`
**R²:** 0.9987 | **Adj R²:** 0.9986
**Residual std error:** 2.697 W
**phi1 (Ph.l1):** 0.4187 (p=4.61e-52)
**omega1 (Tdelta):** 2.3121 (p=2.58e-73)
**omega2 (Gv):** -0.0836 (p=1.17e-78)
**Improvement over OLS:** residual autocorrelation substantially reduced

## 3.7 Model Order Selection (AIC / BIC)
**Figure:** `report/figures/q3_37_aic_bic.pdf`
**Table:** `report/tables/q3_37_ic.tex`
**Best order by AIC:** 6 (AIC = 703.26)
**Best order by BIC:** 5 (BIC = 759.28)
**AIC values (p=1..10):** 810.3, 760.1, 741.5, 722, 709.4, 703.3, 707.3, 710.5, 712.8, 707.6
**BIC values (p=1..10):** 822.7, 782, 772.7, 762.5, 759.3, 762.5, 775.9, 788.5, 800.1, 804.3

## 3.8 Test-Set RMSE vs. Model Order
**Figure:** `report/figures/q3_38_rmse.pdf`
**Table:** `report/tables/q3_38_rmse.tex`
**Best order by RMSE:** 3 (RMSE = 2.8755 W)
**RMSE values (p=1..10):** 3.752, 3.245, 2.875, 3.039, 3.126, 3.332, 3.323, 3.269, 3.26, 3.302

## 3.9 Multi-Step Simulation
**Figure:** `report/figures/q3_39_multistep_simulation.pdf`
**Table (coefficients):** `report/tables/q3_39_arx5_coef.tex`
**Model used:** ARX(5)
**Simulation RMSE — training period:** 2.076 W
**Simulation RMSE — test period:** 3.732 W
**Coefficients:**
  - Ph.l1: 0.2898
  - Ph.l2: 0.0443
  - Ph.l3: 0.0322
  - Ph.l4: 0.1381
  - Ph.l5: 0.0214
  - Tdelta.l0: 0.3702
  - Tdelta.l1: 1.483
  - Tdelta.l2: -0.5818
  - Tdelta.l3: 2.0902
  - Tdelta.l4: -1.4383
  - Gv.l0: -0.0989
  - Gv.l1: -2e-04
  - Gv.l2: 8e-04
  - Gv.l3: 0.0124
  - Gv.l4: 0.0069

## 3.10 Summary
| Criterion | Best order |
|-----------|------------|
| AIC       | 6 |
| BIC       | 5 |
| Test RMSE | 3 |
| Chosen    | 5 (BIC) |

