/**********************************************************************
REPEATED-MEASURES AGREEMENT ANALYSIS
RANDOM-SLOPE HIERARCHICAL MODEL
**********************************************************************/

use "~\OneDrive - University College London\Stata code\Bland Altman PEFR data (repro)\Bland_pefr_long_long_v1.dta", clear

rename test   occasion
rename method device
rename pefr_  value

/**********************************************************************
1. DEFINE PROGRAM
**********************************************************************/

capture program drop mc_agreement

program define mc_agreement, rclass

    version 18.0
    preserve

    // Encode string variables if needed
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

    // Fit random-slope mixed model
    quietly mixed value ///
        i.`dvar' ///
        i.`ovar' ///
        || id: ///
        || id:`dvar' ///
        || id:`ovar', ///
        reml nolog

    // Fixed effects
    scalar intercept     = _b[_cons]
    levelsof `dvar', local(dlevels)
    local d2 : word 2 of `dlevels'
    capture scalar bias_device = _b[2.`dvar']
    if _rc != 0 scalar bias_device = _b[`d2'.`dvar']

    levelsof `ovar', local(olevels)
    local o2 : word 2 of `olevels'
    capture scalar bias_occasion = _b[2.`ovar']
    if _rc != 0 scalar bias_occasion = _b[`o2'.`ovar']

    // Variance components
    predict double fitted_value, fitted
    predict double residual,     residuals

    quietly summarize residual
    scalar sigma2_residual = r(Var)

    bysort id: egen double id_mean = mean(fitted_value)
    quietly summarize id_mean
    scalar sigma2_subject = r(Var)

    scalar sigma2_total = sigma2_subject + sigma2_residual

    quietly summarize fitted_value
    scalar mu = r(mean)

    // Device agreement metrics
    scalar ccc_device       = (2 * sigma2_subject) / (2 * sigma2_total + bias_device^2)
    scalar sd_device        = sqrt(2 * sigma2_residual)
    scalar loa_lower_device = bias_device - 1.96 * sd_device
    scalar loa_upper_device = bias_device + 1.96 * sd_device
    scalar rc_device        = 1.96 * sqrt(2 * sigma2_residual)
    scalar cv_device        = 100 * sqrt(sigma2_residual / (mu^2))

    // Occasion repeatability metrics
    scalar icc_occasion       = sigma2_subject / (sigma2_subject + sigma2_residual)
    scalar sd_occasion        = sqrt(2 * sigma2_residual)
    scalar loa_lower_occasion = bias_occasion - 1.96 * sd_occasion
    scalar loa_upper_occasion = bias_occasion + 1.96 * sd_occasion
    scalar rc_occasion        = 1.96 * sqrt(2 * sigma2_residual)
    scalar cv_occasion        = 100 * sqrt(sigma2_residual / (mu^2))

    // Return results
    return scalar ccc_device         = ccc_device
    return scalar bias_device        = bias_device
    return scalar loa_lower_device   = loa_lower_device
    return scalar loa_upper_device   = loa_upper_device
    return scalar rc_device          = rc_device
    return scalar cv_device          = cv_device

    return scalar icc_occasion       = icc_occasion
    return scalar bias_occasion      = bias_occasion
    return scalar loa_lower_occasion = loa_lower_occasion
    return scalar loa_upper_occasion = loa_upper_occasion
    return scalar rc_occasion        = rc_occasion
    return scalar cv_occasion        = cv_occasion

    restore

end

/**********************************************************************
2. RUN & BOOTSTRAP
**********************************************************************/

bootstrap ///
    Dev_CCC=r(ccc_device)         ///
    Dev_Bias=r(bias_device)       ///
    Dev_LoA_LL=r(loa_lower_device) ///
    Dev_LoA_UL=r(loa_upper_device) ///
    Dev_RC=r(rc_device)           ///
    Dev_CV=r(cv_device)           ///
    Occ_ICC=r(icc_occasion)       ///
    Occ_Bias=r(bias_occasion)     ///
    Occ_LoA_LL=r(loa_lower_occasion) ///
    Occ_LoA_UL=r(loa_upper_occasion) ///
    Occ_RC=r(rc_occasion)         ///
    Occ_CV=r(cv_occasion),        ///
    reps(100)                     ///
    cluster(id)                   ///
    seed(42):                     ///
    mc_agreement

estat bootstrap, all