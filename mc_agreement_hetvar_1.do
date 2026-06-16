/**********************************************************************
REPEATED-MEASURES AGREEMENT ANALYSIS
HETEROGENEOUS RESIDUAL VARIANCE BY DEVICE
(assumes variance differs by device but not occasion)
======================================================================

MODEL
----------------------------------------------------------------------

    value_ijk =
        μ
        + β_k          fixed effect of device k
        + γ_j          fixed effect of occasion j
        + u_i          subject random intercept       ~ N(0, σ²_u)
        + v_ik         subject-by-device slope        ~ N(0, σ²_v)
        + w_ij         subject-by-occasion slope      ~ N(0, σ²_w)
        + e_ijk        residual                       ~ N(0, σ²_k)

Residual variance is allowed to differ by device (σ²_1, σ²_2).
Occasion residual variance is assumed homogeneous (no evidence
for heterogeneity by occasion, and a more complex model is harder
to interpret clinically).

======================================================================
VARIANCE COMPONENT LAYOUT IN e(b)
----------------------------------------------------------------------

After -mixed- with residuals(independent, by(device)), e(b) contains:

  Fixed effects:  _cons  device  occasion
  Random effects (log-SD parameterisation, in order of || specs):
    lns1_1_1:_cons   log-SD of u_i   (subject intercept)
    lns2_1_1:_cons   log-SD of v_ik  (subject-by-device)
    lns3_1_1:_cons   log-SD of w_ij  (subject-by-occasion)
  Residual log-SDs:
    lns?_?_1:_cons   log-SD of e for device 1   (last-1)
    lns?_?_2:_cons   log-SD of e for device 2   (last)

We read these positionally: first lns = σ_u, last two = σ_1, σ_2.
All three intermediate slopes are present but not used directly.

======================================================================
FORMULAS
----------------------------------------------------------------------

DEVICE AGREEMENT

  A difference between devices for subject i on occasion j is:
    D_ij = (μ + β_2 + u_i + v_i2 + w_ij + e_ij2)
          -(μ + β_1 + u_i + v_i1 + w_ij + e_ij1)
         = bias + (v_i2 - v_i1) + (e_ij2 - e_ij1)

  where:
    bias            = β_2 - β_1  (fixed device effect)
    Var(v_i2-v_i1)  = 2σ²_v      (subject-by-device slope cancels on
                                   average but contributes to spread)
    Var(e_ij2-e_ij1)= σ²_1+σ²_2  (independent residuals)

  SD of differences:
    SD_D = sqrt(2σ²_v + σ²_1 + σ²_2)

  Limits of agreement:
    bias ± 1.96 * SD_D

  Lin's CCC:
    cov(Y_k1, Y_k2) = σ²_u         (shared subject intercept only;
                                     subject-by-device slopes are
                                     independent across devices)
    Var(Y_1) = σ²_u + σ²_v + σ²_w + σ²_1
    Var(Y_2) = σ²_u + σ²_v + σ²_w + σ²_2

    CCC = 2*σ²_u / (Var(Y_1) + Var(Y_2) + bias²)
        = 2*σ²_u /
          (2σ²_u + 2σ²_v + 2σ²_w + σ²_1 + σ²_2 + bias²)

  CV (device):
    Mean measurement noise SD = sqrt((σ²_1 + σ²_2) / 2)
    CV = 100 * sqrt((σ²_1 + σ²_2) / 2) / μ

OCCASION REPEATABILITY

  A difference between occasions for subject i on device k is:
    D_ik = (γ_2 - γ_1) + (w_i2 - w_i1) + (e_i2k - e_i1k)

  where:
    bias_occ            = γ_2 - γ_1
    Var(w_i2 - w_i1)    = 2σ²_w
    Var(e_i2k - e_i1k)  = 2 * mean(σ²_1, σ²_2)
                          [residuals are device-specific but the
                           occasion contrast averages across devices]

  SD of occasion differences:
    SD_occ = sqrt(2σ²_w + σ²_1 + σ²_2)

  ICC:
    Between-subject variance seen by a single device on a single
    occasion = σ²_u + σ²_v + σ²_w
    Within-subject (occasion) variance = σ²_w + mean(σ²_1, σ²_2)

    ICC = σ²_u / (σ²_u + σ²_v + σ²_w + mean(σ²_1, σ²_2))

    Rationale: ICC is the expected correlation between two measurements
    on the same subject on different occasions on the same device.
    The shared components are σ²_u only; σ²_v is shared within device
    but not part of the within-occasion total variance denominator;
    the denominator is total variance for a single observation.

  CV (occasion):
    CV = 100 * sqrt(σ²_w + mean(σ²_1, σ²_2)) / μ

======================================================================
OUTPUTS
----------------------------------------------------------------------
    Lin's CCC, bias, LoA, RC, CV       (device agreement)
    ICC, bias, LoA, RC, CV             (occasion repeatability)
    Variance components with 95% CIs   (from bootstrap)
**********************************************************************/

/**********************************************************************
0. LOAD DATA
**********************************************************************/

use ///
    "~\OneDrive - University College London\Stata code\Bland Altman PEFR data (repro)\Bland_pefr_long_long_v1.dta", ///
    clear

rename test   occasion
rename method device
rename pefr_  value


/**********************************************************************
1. MAIN PROGRAM
**********************************************************************/

capture program drop mc_agreement_hetvar

program define mc_agreement_hetvar, rclass

    version 18.0
    preserve

    // ----------------------------------------------------------------
    // Encode string variables if needed
    // ----------------------------------------------------------------

    capture confirm numeric variable device
    if _rc != 0 {
        encode device, gen(device_n)
        local dvar device_n
    }
    else {
        local dvar device
    }

    capture confirm numeric variable occasion
    if _rc != 0 {
        encode occasion, gen(occasion_n)
        local ovar occasion_n
    }
    else {
        local ovar occasion
    }

    // ----------------------------------------------------------------
    // Fit model: heterogeneous residual variance by device
    // ----------------------------------------------------------------

    quietly mixed value ///
        i.`dvar' ///
        i.`ovar' ///
        || id: ///
        || id: `dvar', nocons ///
        || id: `ovar', nocons ///
        reml nolog ///
        residuals(independent, by(`dvar'))

    // ----------------------------------------------------------------
    // Fixed effects
    // ----------------------------------------------------------------

    quietly levelsof `dvar', local(dlevels)
    local d2 : word 2 of `dlevels'
    capture scalar bias_device = _b[2.`dvar']
    if _rc != 0 {
        scalar bias_device = _b[`d2'.`dvar']
    }

    quietly levelsof `ovar', local(olevels)
    local o2 : word 2 of `olevels'
    capture scalar bias_occasion = _b[2.`ovar']
    if _rc != 0 {
        scalar bias_occasion = _b[`o2'.`ovar']
    }

    // ----------------------------------------------------------------
    // Grand mean from fixed-effects prediction only
    // ----------------------------------------------------------------

    predict double xb, xb
    quietly summarize xb
    scalar mu = r(mean)

    // ----------------------------------------------------------------
    // Extract variance components from e(b)
    //
    // Confirmed layout from matrix list e(b) — 10 columns total:
    //
    //   col 1-5  : fixed effects (device, occasion, _cons)
    //   col 6    : lns1_1_1:_cons   log-SD of u_i  (subject intercept)
    //   col 7    : lns1_2_1:_cons   log-SD of v_ik (subject-by-device)
    //   col 8    : lns1_3_1:_cons   log-SD of w_ij (subject-by-occasion)
    //   col 9    : lnsig_e:_cons    log-SD of e for device 1
    //   col 10   : r_lns2ose:_cons  log-SD of e for device 2
    //
    // The three random-effect log-SDs all share the prefix "lns1_"
    // and variable "_cons". The two residual log-SDs use different
    // equation name patterns ("lnsig_e" and "r_lns2ose").
    //
    // Strategy: collect ALL columns with variable "_cons" that are
    // not fixed effects (i.e. equation name != "value"). Split into:
    //   - random-effect log-SDs : equation name matches "^lns1_"
    //   - residual log-SDs      : the remaining two (any other prefix)
    // ----------------------------------------------------------------

    tempname bmat
    matrix `bmat' = e(b)

    local ncols   = colsof(`bmat')
    local eqnames : coleq    `bmat'
    local vnnames : colnames `bmat'

    local re_cols  ""   // random-effect log-SDs (lns1_*)
    local res_cols ""   // residual log-SDs      (lnsig_e, r_lns2ose, etc.)

    forvalues c = 1/`ncols' {
        local eq  : word `c' of `eqnames'
        local var : word `c' of `vnnames'
        if "`var'" == "_cons" & "`eq'" != "value" {
            if regexm("`eq'", "^lns1_") {
                local re_cols "`re_cols' `c'"
            }
            else if regexm("`eq'", "^lns") | regexm("`eq'", "^r_lns") {
                local res_cols "`res_cols' `c'"
            }
        }
    }

    local nre  = wordcount("`re_cols'")
    local nres = wordcount("`res_cols'")

    if `nre' < 3 | `nres' < 2 {
        di as error "Could not locate all variance parameters."
        di as error "Random-effect log-SDs found (need 3): `nre'"
        di as error "Residual log-SDs found (need 2): `nres'"
        di as error "Run -matrix list e(b)- to inspect."
        exit 499
    }

    // Random effects: subject intercept, by-device slope, by-occasion slope
    local pos_u : word 1 of `re_cols'
    local pos_v : word 2 of `re_cols'
    local pos_w : word 3 of `re_cols'

    // Residuals: device 1 (lnsig_e), device 2 (r_lns2ose)
    local pos_r1 : word 1 of `res_cols'
    local pos_r2 : word 2 of `res_cols'

    // Print for transparency
    foreach p in pos_u pos_v pos_w pos_r1 pos_r2 {
        local col ``p''
        local eq  : word `col' of `eqnames'
        di as text "  [`p'] [`eq':_cons] = " ///
            %8.5f `bmat'[1,`col'] "  (SD = " %7.4f exp(`bmat'[1,`col']) ")"
    }

    scalar sigma2_u  = exp(`bmat'[1, `pos_u'])^2
    scalar sigma2_v  = exp(`bmat'[1, `pos_v'])^2
    scalar sigma2_w  = exp(`bmat'[1, `pos_w'])^2
    scalar sigma2_r1 = exp(`bmat'[1, `pos_r1'])^2
    scalar sigma2_r2 = exp(`bmat'[1, `pos_r2'])^2

    // ----------------------------------------------------------------
    // Device agreement statistics
    // ----------------------------------------------------------------

    scalar sd_D_device      = sqrt(2*sigma2_v + sigma2_r1 + sigma2_r2)
    scalar var_Y1           = sigma2_u + sigma2_v + sigma2_w + sigma2_r1
    scalar var_Y2           = sigma2_u + sigma2_v + sigma2_w + sigma2_r2

    scalar ccc_device       = (2 * sigma2_u) / ///
                              (var_Y1 + var_Y2 + bias_device^2)
    scalar loa_lower_device = bias_device - 1.96 * sd_D_device
    scalar loa_upper_device = bias_device + 1.96 * sd_D_device
    scalar rc_device        = 1.96 * sd_D_device
    scalar cv_device        = 100 * sqrt((sigma2_r1 + sigma2_r2) / 2) / mu

    // ----------------------------------------------------------------
    // Occasion repeatability statistics
    // ----------------------------------------------------------------

    scalar mean_sigma2_r    = (sigma2_r1 + sigma2_r2) / 2
    scalar sd_D_occasion    = sqrt(2*sigma2_w + sigma2_r1 + sigma2_r2)

    scalar icc_occasion       = sigma2_u / ///
                                (sigma2_u + sigma2_v + sigma2_w + mean_sigma2_r)
    scalar loa_lower_occasion = bias_occasion - 1.96 * sd_D_occasion
    scalar loa_upper_occasion = bias_occasion + 1.96 * sd_D_occasion
    scalar rc_occasion        = 1.96 * sd_D_occasion
    scalar cv_occasion        = 100 * sqrt(sigma2_w + mean_sigma2_r) / mu

    // ----------------------------------------------------------------
    // Return
    // ----------------------------------------------------------------

    return scalar sigma2_u   = sigma2_u
    return scalar sigma2_v   = sigma2_v
    return scalar sigma2_w   = sigma2_w
    return scalar sigma2_r1  = sigma2_r1
    return scalar sigma2_r2  = sigma2_r2
    return scalar mu         = mu

    return scalar ccc_device       = ccc_device
    return scalar bias_device      = bias_device
    return scalar loa_lower_device = loa_lower_device
    return scalar loa_upper_device = loa_upper_device
    return scalar rc_device        = rc_device
    return scalar cv_device        = cv_device

    return scalar icc_occasion       = icc_occasion
    return scalar bias_occasion      = bias_occasion
    return scalar loa_lower_occasion = loa_lower_occasion
    return scalar loa_upper_occasion = loa_upper_occasion
    return scalar rc_occasion        = rc_occasion
    return scalar cv_occasion        = cv_occasion

    restore

end

/**********************************************************************
2. RUN ONCE — inspect variance components and check ordering
**********************************************************************/

mc_agreement_hetvar

di ""
di "==============================================================="
di "VARIANCE COMPONENTS"
di "==============================================================="
di "  sigma2_u  (subject intercept)     " %10.4f r(sigma2_u)
di "  sigma2_v  (subject-by-device)     " %10.4f r(sigma2_v)
di "  sigma2_w  (subject-by-occasion)   " %10.4f r(sigma2_w)
di "  sigma2_r1 (residual, device 1)    " %10.4f r(sigma2_r1)
di "  sigma2_r2 (residual, device 2)    " %10.4f r(sigma2_r2)
di "  mu        (grand mean)            " %10.4f r(mu)
di "==============================================================="

/**********************************************************************
3. BOOTSTRAP
**********************************************************************/

bootstrap ///
    Dev_CCC=r(ccc_device)            ///
    Dev_Bias=r(bias_device)          ///
    Dev_LoA_LL=r(loa_lower_device)   ///
    Dev_LoA_UL=r(loa_upper_device)   ///
    Dev_RC=r(rc_device)              ///
    Dev_CV=r(cv_device)              ///
    Occ_ICC=r(icc_occasion)          ///
    Occ_Bias=r(bias_occasion)        ///
    Occ_LoA_LL=r(loa_lower_occasion) ///
    Occ_LoA_UL=r(loa_upper_occasion) ///
    Occ_RC=r(rc_occasion)            ///
    Occ_CV=r(cv_occasion)            ///
    VC_u=r(sigma2_u)                 ///
    VC_v=r(sigma2_v)                 ///
    VC_w=r(sigma2_w)                 ///
    VC_r1=r(sigma2_r1)               ///
    VC_r2=r(sigma2_r2),              ///
    reps(100)                        ///
    cluster(id)                      ///
    seed(42):                        ///
    mc_agreement_hetvar

estat bootstrap, all
