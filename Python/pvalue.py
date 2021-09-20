import numpy
from scipy import integrate
from scipy.stats import binom
import statistics
import utils

def stats2pvalue(i, Tp, M, formula = "exact", alternative = "right_tail"):
    available_alternatives = ["left_tail", "right_tail", "two_tail"]
    alternative = utils.match_arg(alternative, available_alternatives)
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
    type = utils.match_arg(formula, available_formulae)
    if type == "estimate":
        return b / B
    elif type == "upper_bound":
        return (b + 1) / (B + 1)
    else:
        return phipson_smyth_pvalue(b, B, M)

def phipson_smyth_pvalue(b, B, M):
    if M <= 10000:
        pt = numpy.linspace(1, M + 1, M + 1) / (M + 1)
        return statistics.mean(binom.cdf(k = b, n = B, p = pt))
  
    corr = integrate.quadrature(binom.cdf, 0, 0.5 / (M + 1), args = {k: b, n: B})
    return (b + 1) / (B + 1) - corr[0]

def combine_pvalues(p, combine_with = "tippett"):
    available_combine_withs = ["fisher", "tippett"]
    combine_with = utils.match_arg(combine_with, available_combine_withs)
    if combine_with == "tippett":
        return 1 - min(p)
    else:
        return -2 * sum(numpy.log(p))
