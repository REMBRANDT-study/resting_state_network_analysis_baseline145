
################################################################################
# calc_HC calculates the higher criticism statistic on a set of p values
# In:  p - p-values to be highly criticized
#      alpha - absolute cutoff for p values considered for HC
#      k1 - fraction or number of p values to keep
#      emp - flag for using variance of empirical distribution
#      plot - flag to plot summary items
# Out: p - maximum p value under cutoff
#      hc - maximum HC statistic
################################################################################

calc_HC = function(p, alpha, k1, emp, plot, debug) {

    if(missing(debug)) { debug = FALSE }
    # Sort p values in ascending order
    n = length(p)
    p.sorted = sort(p)
    # Determine cutoff
    if(missing(k1)) {
        print_msg('WARNING! k1 missing from calc_HC')
        k1 = length(p)
        if(missing(alpha)) {
            alpha = 1
        } else {
            k1 = max(which(p.sorted < alpha))
        }
    } else {
        if(missing(alpha)) { alpha = 1 }
        if(k1 <= 1) { k1 = floor(k1 * length(p)) }
    }
    p.sorted = p.sorted[1:k1]

    if(debug) { print_msg(sprintf('alpha %0.2f k1 %i emp %s', alpha, k1, emp)) }
    # Assume variance of theoretical null
    if(missing(emp)) { emp = FALSE }
    if(missing(plot)) { plot = FALSE }

    if(min(p.sorted) < alpha) {

        # Typical case
        i.vals = 1:k1

        # Calculate higher criticism
        # page 966 of Donoho2004 (unnumbered equation for HC^*_n)
        hc1 = rep(NA, k1)
        for (i in i.vals) {
            hc1[i] = sqrt(n) * (i / n - p.sorted[i]) /
                sqrt(p.sorted[i] * (1 - p.sorted[i])) }
        # eq (1) of Donoho2018
        hc2 = rep(NA, k1)
        for (i in i.vals) {
            hc2[i] = sqrt(n) * (i / n - p.sorted[i]) /
                sqrt(i / n * (1 - i / n)) }
        # select theoretical or empirical results
        if(emp) { hc = hc1 } else { hc = hc2 }

        # Remove points below Bonferroni correction
        hc[p.sorted < (0.05 / n)] = 0
        # Summarize results
        p.max = p.sorted[which.max(hc)]

        # Plot HC
        if(plot) {
            dev.new()
            plot(p.sorted, hc,
                pch='.', col='blue', cex=2,
                xlim=c(0,alpha),
                xlab='Fraction of Tests (ordered by p-value)',
                ylab='HC Statistic')
            lines(c(0,alpha), c(2,2), col='red', lwd=3) }

    } else {

        # No p values below cutoff
        hc = 0
        p.max = 0

    }

    return(max(hc))

}
