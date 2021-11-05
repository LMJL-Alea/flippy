import sys
import numpy as np
from triarray import TriMatrix
from .utils import set_doc

two_sample_stats_doc = """
Test Statistics for the Two-Sample Problem

This is a collection of functions that provide test statistics to be used
into the permutation scheme for performing two-sample testing. These test
statistics can be divided into two categories: traditional statistics that
use empirical moments and inter-point statistics that only rely on pairwise
dissimilarities between data points.

Parameters:
    data (list or triarray.TriMatrix): Either a list of the `n1 + n2` 
concatenated observations with the original `n1` observations from the first 
sample on top and the original `n2` observations from the second sample below. 
Or a dissimilarity matrix stored as a `triarray.TriMatrix` object for all 
inter-point statistics whose function name should end with `_ip()`.
    indices1 (range or list of integers): An integer vector specifying the 
indices in `data` that are considered to belong to the first sample.
    alpha (double): A scalar value in (0, 2] specifying the power to which the 
dissimilarities should be elevated in the computation of the `stat_energy_ip()` 
statistic. Default is `1`.
    standardize (bool): A boolean specifying whether the distance between 
medoids in the `stat_dom_ip()` function should be normalized by the pooled 
corresponding variances. Default is `True`.
    **kwargs: Extra parameters specific to some statistics.

Value:
A real scalar giving the value of test statistic for the permutation specified 
by the integer vector `indices`.

Traditional Test Statistics:
    `stat_hotelling()` implements Hotelling's \eqn{T^2} statistic for 
multivariate data with \eqn{p < n}.
    `stat_student()` or `stat_t()` implements Student's statistic (originally 
assuming equal variances and thus using the pooled empirical variance estimator).
    `stat_welch()` implements Student-Welch statistic which is essentially a 
modification of Student's statistic accounting for unequal variances.
    `stat_fisher()` or `stat_f()` implements Fisher's variance ratio statistic.
    `stat_mean()` implements a statistic that computes the difference between 
the means.
    `stat_bs()` implements the statistic proposed by Bai & Saranadasa (1996) for 
high-dimensional multivariate data.

Inter-Point Test Statistics:
    `stat_student_ip()` or `stat_t_ip()` implements a Student-like test 
statistic based on inter-point distances only as described in Lovato et al. 
(2020).
    `stat_fisher_ip()` or `stat_f_ip()` implements a Fisher-like test statistic 
based on inter-point distances only as described in Lovato et al. (2020).
    `stat_bg_ip()` implements the statistic proposed by Biswas & Ghosh (2014).
    `stat_energy_ip()` implements the class of energy-based statistics as 
described in Székely & Rizzo (2013);
    `stat_cq_ip()` implements the statistic proposed by Chen & Qin (2010).
    `stat_mod_ip()` implements a statistic that computes the mean of inter-point 
distances.
    `stat_dom_ip()` implements a statistic that computes the distance between 
the medoids of the two samples, possibly standardized by the pooled 
corresponding variances.

References:
    Bai, Z., & Saranadasa, H. (1996). Effect of high dimension: by an example of
a two sample problem. Statistica Sinica, 311-329.
    Lovato, I., Pini, A., Stamm, A., & Vantini, S. (2020). Model-free two-sample
test for network-valued data. Computational Statistics & Data Analysis, 144,
106896.
    Biswas, M., & Ghosh, A. K. (2014). A nonparametric two-sample test applicable
to high dimensional data. Journal of Multivariate Analysis, 123, 160-171.
    Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of
statistics based on distances. Journal of statistical planning and inference,
143(8), 1249-1272.
    Chen, S. X., & Qin, Y. L. (2010). A two-sample test for high-dimensional data
with applications to gene-set testing. The Annals of Statistics, 38(2),
808-835.

Examples:
import numpy as np
from scipy.spatial.distance import pdist
from triarray import TriMatrix

n = 10
x = np.random.rand(n, 1)
y = np.random.rand(n, 1)
dists = pdist(np.vstack([x, y]))
D = TriMatrix(dists, upper = True, diag_val = 0)

x = list(x)
y = list(y)

stat_welch(x + y, range(n))
stat_t(x + y, range(n))
stat_f(x + y, range(n))
stat_mean(x + y, range(n))
stat_hotelling(x + y, range(n))
stat_bs(x + y, range(n))

stat_t_ip(D, range(n))
stat_f_ip(D, range(n))
stat_bg_ip(D, range(n))
stat_energy_ip(D, range(n))
stat_cq_ip(D, range(n))
stat_mod_ip(D, range(n))
stat_dom_ip(D, range(n))
"""

@set_doc(two_sample_stats_doc)
def stat_welch(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = [data[i] for i in l["idx1"]]
    x2 = [data[i] for i in l["idx2"]]
    m1 = np.mean(x1)
    m2 = np.mean(x2)
    v1 = np.var(x1, ddof = 1)
    v2 = np.var(x2, ddof = 1)
    stderr = np.sqrt(v1 / n1 + v2 / n2)
    return (m1 - m2) / stderr

@set_doc(two_sample_stats_doc)
def stat_student(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = [data[i] for i in l["idx1"]]
    x2 = [data[i] for i in l["idx2"]]
    m1 = np.mean(x1)
    m2 = np.mean(x2)
    v1 = np.var(x1)
    v2 = np.var(x2)
    sig2 = (n1 * v1 + n2 * v2) / (n1 + n2 - 2)
    stderr = np.sqrt(sig2 * (1 / n1 + 1 / n2))
    return (m2 - m1) / stderr
 
@set_doc(two_sample_stats_doc)
def stat_t(data, indices1, **kwargs):
    return stat_student(data, indices1, **kwargs)
 
@set_doc(two_sample_stats_doc)
def stat_fisher(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = [data[i] for i in l["idx1"]]
    x2 = [data[i] for i in l["idx2"]]
    v1 = np.var(x1, ddof = 1)
    v2 = np.var(x2, ddof = 1)
    return v2 / v1

@set_doc(two_sample_stats_doc)
def stat_f(data, indices1, **kwargs):
    return stat_fisher(data, indices1, **kwargs)

@set_doc(two_sample_stats_doc)
def stat_mean(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = [data[i] for i in l["idx1"]]
    x2 = [data[i] for i in l["idx2"]]
    m1 = np.mean(x1)
    m2 = np.mean(x2)
    return m1 - m2

@set_doc(two_sample_stats_doc)
def stat_hotelling(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = np.array([data[i] for i in l["idx1"]])
    x2 = np.array([data[i] for i in l["idx2"]])
    m1 = np.mean(x1, 0)
    m2 = np.mean(x2, 0)
    v1 = np.cov(x1.transpose(), ddof = 1)
    v2 = np.cov(x2.transpose(), ddof = 1)
    vp = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    vinv = np.linalg.inv(vp)
    d = m1 - m2
    return d.transpose() @ vinv @ d

@set_doc(two_sample_stats_doc)
def stat_bs(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    x1 = np.array([data[i] for i in l["idx1"]])
    x2 = np.array([data[i] for i in l["idx2"]])
    m1 = np.mean(x1, 0)
    m2 = np.mean(x2, 0)
    v1 = np.cov(x1.transpose(), ddof = 1)
    v2 = np.cov(x2.transpose(), ddof = 1)
    vn = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2)
    d = m1 - m2
    return d.transpose() @ d - (n1 + n2) / (n1 * n2) * sum(np.diag(vn))

@set_doc(two_sample_stats_doc)
def stat_student_ip(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D11 = [[data.get_item(l["idx1"][i], l["idx1"][j]) 
        for j in range(i + 1, n1)] for i in range(n1)]
    M11 = (np.array([item for sublist in D11 for item in sublist])**2).mean()
    D22 = [[data.get_item(l["idx2"][i], l["idx2"][j]) 
        for j in range(i + 1, n2)] for i in range(n2)]
    M22 = (np.array([item for sublist in D22 for item in sublist])**2).mean()
    D12 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx1"]]
    M12 = (np.array(D12)**2).mean()
    v1 = M11 / 2
    v2 = M22 / 2
    numerator = M12 - v1 - v2
    vn = (v1 / n1 + v2 / n2)
    if vn < np.sqrt(np.finfo(float).eps):
        return numerator
    return numerator / vn

@set_doc(two_sample_stats_doc)
def stat_t_ip(data, indices1, **kwargs):
    return stat_student_ip(data, indices1, **kwargs)

@set_doc(two_sample_stats_doc)
def stat_fisher_ip(data, indices1, **kwarg):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D11 = [[data.get_item(l["idx1"][i], l["idx1"][j]) 
        for j in range(i + 1, n1)] for i in range(n1)]
    M11 = (np.array([item for sublist in D11 for item in sublist])**2).mean()
    D22 = [[data.get_item(l["idx2"][i], l["idx2"][j]) 
        for j in range(i + 1, n2)] for i in range(n2)]
    M22 = (np.array([item for sublist in D22 for item in sublist])**2).mean()
    v1 = M11 / 2
    v2 = M22 / 2
    vmin = min(v1, v2)
    vmax = max(v1, v2)
    if vmin < np.sqrt(np.finfo(float).eps):
        return vmax
    return vmax / vmin

@set_doc(two_sample_stats_doc)
def stat_f_ip(data, indices1, **kwargs):
    return stat_fisher_ip(data, indices1, **kwargs)

@set_doc(two_sample_stats_doc)
def stat_bg_ip(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D11 = [[data.get_item(l["idx1"][i], l["idx1"][j]) 
        for j in range(i + 1, n1)] for i in range(n1)]
    M11 = np.array([item for sublist in D11 for item in sublist]).mean()
    D22 = [[data.get_item(l["idx2"][i], l["idx2"][j]) 
        for j in range(i + 1, n2)] for i in range(n2)]
    M22 = np.array([item for sublist in D22 for item in sublist]).mean()
    D12 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx1"]]
    M12 = np.array(D12).mean()
    return (M11 - M12)**2 + (M22 - M12)**2

@set_doc(two_sample_stats_doc)
def stat_energy_ip(data, indices1, alpha = 1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D11 = [[data.get_item(i, j) for j in l["idx1"]] for i in l["idx1"]]
    M11 = (np.array(D11)**alpha).mean()
    D22 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx2"]]
    M22 = (np.array(D22)**alpha).mean()
    D12 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx1"]]
    M12 = (np.array(D12)**alpha).mean()
    return M12 - (M11 + M22) / 2

@set_doc(two_sample_stats_doc)
def stat_cq_ip(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D11 = [[data.get_item(i, j) for j in l["idx1"] if j != i] for i in l["idx1"]]
    M11 = (np.array(D11)).mean()
    D22 = [[data.get_item(i, j) for j in l["idx2"] if j != i] for i in l["idx2"]]
    M22 = (np.array(D22)).mean()
    D12 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx1"]]
    M12 = (np.array(D12)).mean()
    return M11 + M22 - 2 * M12
  
@set_doc(two_sample_stats_doc)
def stat_mod_ip(data, indices1, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    D12 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx1"]]
    M12 = (np.array(D12)).mean()
    return M12

@set_doc(two_sample_stats_doc)
def stat_dom_ip(data, indices1, standardize = True, **kwargs):
    l = two_sample_prep(data, indices1)
    n1 = len(l["idx1"])
    n2 = len(l["idx2"])
    
    D11 = [[data.get_item(i, j) for j in l["idx1"]] for i in l["idx1"]]
    ssd1_vec = (np.array(D11)**2).sum(axis = 0)
    km1 = l["idx1"][ssd1_vec.argmin()]
    
    D22 = [[data.get_item(i, j) for j in l["idx2"]] for i in l["idx2"]]
    ssd2_vec = (np.array(D22)**2).sum(axis = 0)
    km2 = l["idx2"][ssd2_vec.argmin()]
    
    stat = data.get_item(km1, km2)
    
    if not standardize:
        return stat
    
    ssd1 = min(ssd1_vec)
    ssd2 = min(ssd2_vec)
    vp = (ssd1 + ssd2) / (n1 + n2 - 2)
    return stat / np.sqrt(vp)

def two_sample_prep(data, indices1):
    if isinstance(data, list):
        n = len(data)
    elif isinstance(data, TriMatrix):
        n = data.size
    else:
        sys.exit("The `data` input should be of class either list or dist.")
    indices1 = list(indices1)
    indices1.sort()
    indices2 = list(set(range(n)) - set(indices1))
    indices2.sort()
    return {"idx1": indices1, "idx2": indices2}
