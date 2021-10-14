import numpy as np
from scipy import integrate
from scipy.stats import binom
import statistics
from .utils import match_arg

def stats2pvalue(i, Tp, M, formula = "exact", alternative = "right_tail"):
    available_alternatives = ["left_tail", "right_tail", "two_tail"]
    alternative = match_arg(alternative, available_alternatives)
    T0 = Tp[i]
    B = len(Tp) - 1
    if alternative == "right_tail":
        b = sum(Tp > T0)
    elif alternative == "left_tail":
        b = sum(Tp < T0)
    else:
        b = 2 * min(sum(Tp <= T0), sum(Tp > T0)) - 1
    return get_p(b, B, M, formula)

def get_p(b, B, M, formula):
    available_formulae = ["exact", "upper_bound", "estimate"]
    type = match_arg(formula, available_formulae)
    if type == "estimate":
        return b / B
    elif type == "upper_bound":
        return (b + 1) / (B + 1)
    else:
        return phipson_smyth_pvalue(b, B, M)

def phipson_smyth_pvalue(b, B, M):
    if M <= 10000:
        pt = np.linspace(1, M + 1, M + 1) / (M + 1)
        return statistics.mean(binom.cdf(k = b, n = B, p = pt))
  
    corr = integrate.quadrature(binom.cdf, 0, 0.5 / (M + 1), args = {k: b, n: B})
    return (b + 1) / (B + 1) - corr[0]

def combine_pvalues(p, combine_with = "tippett"):
    available_combine_withs = ["fisher", "tippett"]
    combine_with = match_arg(combine_with, available_combine_withs)
    if combine_with == "tippett":
        return 1 - min(p)
    else:
        return -2 * sum(np.log(p))

def run_permutation_scheme(formula, alternative, stats, B, perm_data, stat_data, M, combine_with, **kwargs):
    formula = match_arg(formula, c("exact", "upper_bound", "estimate"))
    alternative = match_arg(alternative, c("left_tail", "right_tail", "two_tail"))
    nstats = len(stats)
    npc = nstats > 1
    if npc == True:
        altern = "right_tail"
    else:
        altern = alternative
    if len(alternative) == 1:
        alternative = np.repeat(alternative, nstats)

    if npc == False:
        Tp = [get_permuted_statistic(i, 
                                     perm_data = perm_data, 
                                     stat_data = stat_data, 
                                     stat_fun = stats[0], 
                                     **kwargs) for i in np.arange(B+1)]
    else:
        Tp = [[get_permuted_statistic(i, 
                                     perm_data = perm_data, 
                                     stat_data = stat_data, 
                                     stat_fun = s, 
                                     **kwargs) for i in np.arange(B+1)] for s in stats]
        Tp = [[stats2pvalue(i, Tp = tval, M = M, formula = "upper_bound", alternative = aval) for tval, aval in zip(Tp, alternative)] for i in np.arange(B+1)+1]
        Tp = [combine_pvalues(tvals, combine_with = combine_with) for tvals in Tp]

    return {
        observed: Tp[0],
        pvalue: stats2pvalue(0, Tp, M, formula = formula, alternative = altern),
        null_distribution: Tp[1:],
        permutations: perm_data
    }

def get_permuted_statistic(i, perm_data, stat_data, stat_fun, **kwargs):
    stat_fun(stat_data, perm_data[, i + 1], **kwargs)
}
