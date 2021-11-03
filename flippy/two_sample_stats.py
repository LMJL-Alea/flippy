import sys
import numpy as np
from triarray import TriMatrix
from .utils import set_doc

two_sample_stats_doc = """Test Statistics for the Two-Sample Problem

Description:
This is a collection of functions that provide test statistics to be used
into the permutation scheme for performing two-sample testing. These test
statistics can be divided into two categories: traditional statistics that
use empirical moments and inter-point statistics that only rely on pairwise
dissimilarities between data points.

Arguments:
- data: Either a list of the `n1 + n2` concatenated observations with the
original `n1` observations from the first sample on top and the original
`n2` observations from the second sample below. Or a dissimilarity matrix
stored as a \code{\link[stats]{dist}} object for all inter-point statistics
whose function name should end with `_ip()`.
- indices1: An integer vector specifying the indices in `data` that are
considered to belong to the first sample.
- alpha: A scalar value specifying the power to which the dissimilarities
should be elevated in the computation of the inter-point energy statistic.
Default is `1L`.
- standardize: A boolean specifying whether the distance between medoids
in the \code{\link{stat_dom_ip}} function should be normalized by the
pooled corresponding variances. Default is `TRUE`.
- **kwargs: Extra parameters specific to some statistics.

Value:
A real scalar giving the value of test statistic for the permutation
specified by the integer vector `indices`.

Traditional Test Statistics:
- \code{\link{stat_hotelling}} implements Hotelling's \eqn{T^2} statistic for
' multivariate data with \eqn{p < n}.
- \code{\link{stat_student}} or \code{\link{stat_t}} implements Student's
statistic (originally assuming equal variances and thus using the pooled
empirical variance estimator). See \code{\link[stats]{t.test}} for details.
- \code{\link{stat_welch}} implements Student-Welch statistic which is
essentially a modification of Student's statistic accounting for unequal
variances. See \code{\link[stats]{t.test}} for details.
- \code{\link{stat_fisher}} or \code{\link{stat_f}} implements Fisher's
variance ratio statistic. See \code{\link[stats]{var.test}} for details.
- \code{\link{stat_mean}} implements a statistic that computes the difference
between the means.
- \code{\link{stat_bs}} implements the statistic proposed by Bai & Saranadasa
' (1996) for high-dimensional multivariate data.

Inter-Point Test Statistics:
- \code{\link{stat_student_ip}} or \code{\link{stat_t_ip}} implements a
Student-like test statistic based on inter-point distances only as described
in Lovato et al. (2020).
- \code{\link{stat_fisher_ip}} or \code{\link{stat_f_ip}} implements a
Fisher-like test statistic based on inter-point distances only as described
in Lovato et al. (2020).
- \code{\link{stat_bg_ip}} implements the statistic proposed by Biswas &
Ghosh (2014).
- \code{\link{stat_energy_ip}} implements the class of energy-based
statistics as described in Székely & Rizzo (2013);
- \code{\link{stat_cq_ip}} implements the statistic proposed by Chen & Qin
(2010).
- \code{\link{stat_mod_ip}} implements a statistic that computes the mean of
inter-point distances.
- \code{\link{stat_dom_ip}} implements a statistic that computes the distance
between the medoids of the two samples, possibly standardized by the pooled
corresponding variances.

References:
- Bai, Z., & Saranadasa, H. (1996). Effect of high dimension: by an example of
a two sample problem. Statistica Sinica, 311-329.
- Lovato, I., Pini, A., Stamm, A., & Vantini, S. (2020). Model-free two-sample
test for network-valued data. Computational Statistics & Data Analysis, 144,
106896.
- Biswas, M., & Ghosh, A. K. (2014). A nonparametric two-sample test applicable
to high dimensional data. Journal of Multivariate Analysis, 123, 160-171.
- Székely, G. J., & Rizzo, M. L. (2013). Energy statistics: A class of
statistics based on distances. Journal of statistical planning and inference,
143(8), 1249-1272.
- Chen, S. X., & Qin, Y. L. (2010). A two-sample test for high-dimensional data
with applications to gene-set testing. The Annals of Statistics, 38(2),
808-835.

Examples:
n = 10
mx = 0
sigma = 1
delta = 10
my = mx + delta
x = rnorm(n = n, mean = mx, sd = sigma)
y = rnorm(n = n, mean = my, sd = sigma)
D = dist(c(x, y))

x = as.list(x)
y = as.list(y)

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

# #' @rdname two-sample-stats
# #' @export
# stat_student_ip <- function(data, indices1, ...) {
#   l <- two_sample_prep(data, indices1)
#   stat_student_impl(data, l$idx1, l$idx2)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_t_ip <- stat_student_ip
# 
# #' @rdname two-sample-stats
# #' @export
# stat_fisher_ip <- function(data, indices1, ...) {
#   l <- two_sample_prep(data, indices1)
#   stat_fisher_impl(data, l$idx1, l$idx2)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_f_ip <- stat_fisher_ip
# 
# #' @rdname two-sample-stats
# #' @export
# stat_bg_ip <- function(data, indices1, ...) {
#   l <- two_sample_prep(data, indices1)
#   stat_bg_impl(data, l$idx1, l$idx2)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_energy_ip <- function(data, indices1, alpha = 1L, ...) {
#   l <- two_sample_prep(data, indices1)
#   stat_energy_impl(data, l$idx1, l$idx2, alpha)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_cq_ip <- function(data, indices1, ...) {
#   l <- two_sample_prep(data, indices1)
#   stat_cq_impl(data, l$idx1, l$idx2)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_mod_ip <- function(data, indices1, ...) {
#   l <- two_sample_prep(data, indices1)
#   l <- purrr::cross(l)
#   dist_values <- purrr::map_dbl(l, ~ getElement(
#     distObject = data,
#     rowIndex = .x$idx1,
#     colIndex = .x$idx2
#   ))
#   mean(dist_values)
# }
# 
# #' @rdname two-sample-stats
# #' @export
# stat_dom_ip <- function(data, indices1, standardize = TRUE, ...) {
#   l <- two_sample_prep(data, indices1)
#   n1 <- length(l$idx1)
#   n2 <- length(l$idx2)
# 
#   ssd1_vec <- purrr::map_dbl(l$idx1, function(idx) {
#     sum(purrr::map_dbl(l$idx1, ~ getElement(
#       distObject = data,
#       rowIndex = idx,
#       colIndex = .x
#     ))^2)
#   })
#   km1 <- l$idx1[which.min(ssd1_vec)]
# 
#   ssd2_vec <- purrr::map_dbl(l$idx2, function(idx) {
#     sum(purrr::map_dbl(l$idx2, ~ getElement(
#       distObject = data,
#       rowIndex = idx,
#       colIndex = .x
#     ))^2)
#   })
#   km2 <- l$idx2[which.min(ssd2_vec)]
# 
#   stat <- getElement(data, km1, km2)
# 
#   if (!standardize)
#     return(stat)
# 
#   ssd1 <- min(ssd1_vec)
#   ssd2 <- min(ssd2_vec)
#   pooled_variance <- (ssd1 + ssd2) / (n1 + n2 - 2)
#   stat / sqrt(pooled_variance)
# }

def two_sample_prep(data, indices1):
    if isinstance(data, list):
        n = len(data)
    elif isinstance(data, TriMatrix):
        n = data.size
    else:
        sys.exit("The `data` input should be of class either list or dist.")
    indices1 = list(indices1)
    indices2 = list(set(range(n)) - set(indices1))
    return {"idx1": indices1, "idx2": indices2}
