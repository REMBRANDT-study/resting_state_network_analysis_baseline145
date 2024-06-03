################################################################################
# calc_HC_p_MC uses comparison to theoretic null to determine p value for HC
# In:  hc - HC statistic
#      n_test - number of tests
#      n_sim - number of simulations to run
#      alpha - absolute cutoff for p values considered for HC
#      k1 - fraction or number of p values to keep
#      emp - flag for using variance of empirical distribution
# Out: p - p value
################################################################################

calc_HC_p_MC = function(hc, n_test, n_sim, alpha, k1, emp) {

    alpha_cutoff = FALSE
    if(missing(k1)) {
        k1 = n_test
        if(missing(alpha)) {
            alpha = 1
        } else {
            alpha_cutoff = TRUE
        }
    } else {
        if(missing(alpha)) { alpha = 1 }
        if(k1 <= 1) { k1 = floor(k1 * n_test) }
    }
    if(missing(n_sim)) { n_sim = 10000 }
    if(missing(emp)) { emp = FALSE }
    hc_vals = rep(NA, n_sim)
    if(alpha_cutoff) {
        for(i in 1:n_sim) {
            hc_vals[i] = calc_HC(p=runif(n_test), alpha=alpha, emp=emp) }
    } else {
        for(i in 1:n_sim) {
            hc_vals[i] = calc_HC(p=runif(n_test), k1=k1, emp=emp) }
    }
    p = rep(NA, length(hc))
    for(i in 1:length(hc)) {
        p[i] = 1 - sum(hc[i] > hc_vals) / n_sim }
    return(p)

}
