"""
Helper functions for pyferret stats external functions.
"""
import math
import numpy
import scipy.stats
import scipy.special
import pyferret

# The number of supported distributions
NUM_DISTRIBS = 22

def getdistrib(distribname=None, distribparams=None):
    """
    Creates and returns scipy.stats "frozen" probability distribution
    object.  Converts the "standard" parameters (including ordering)
    for a distribution in order to appropriately call the constructor
    for the scipy.stats frozen distribution object.

    If distribparams is None (or not given), this instead returns a tuple
    of (param_name, param_descript) string pairs.  The value param_name
    is the parameter name and the value param_descript is a description
    of the parameter.

    If distribname is None (or not given), this instead returns a tuple
    of (dist_name, dist_descript) string pairs.  The value dist_name is
    the abbreviated name of the probability distribution and the value
    dist_descript is a full name of the distribution with parameters
    names.

    Arguments:
       distribname - name of the distribution
       distribparams - tuple/list/array of input parameters

    Returns:
       if distribname is None (or not given), a tuple of (dist_name,
           dist_descript) string pairs; otherwise,
       if distribparams is None (or not given), a tuple of (param_name,
           param_descript) string pairs; otherwise,
       the scipy.stats "frozen" distribution object described by
           distribname and distribparams

    Raises:
       ValueError if the distribution name is not recognized by this routine,
                  if the incorrect number of parameters are given, or
                  if the distribution parameters are invalid
    """
    if distribname == None:
        return ( ( "beta", "Beta(ALPHA, BETA)", ),
                 ( "binom", "Binomial(N, P)", ),
                 ( "cauchy", "Cauchy(M, GAMMA)", ),
                 ( "chi", "Chi(DF)", ),
                 ( "chi2", "Chi-Square(DF)", ),
                 ( "expon", "Exponential(LAMBDA)", ),
                 ( "exponweib", "Exponentiated-Weibull(K, LAMBDA, ALPHA)", ),
                 ( "f", "F or Fisher(DFN, DFD)", ),
                 ( "gamma", "Gamma(ALPHA, THETA)" ),
                 ( "geom", "Geometric or Shifted-Geometric(P)", ),
                 ( "hypergeom", "Hypergeometric(NGOOD, NTOTAL, NDRAWN)", ),
                 ( "invgamma", "Inverse-Gamma(ALPHA, BETA)", ),
                 ( "laplace", "Laplace(MU, B)", ),
                 ( "lognorm", "Log-Normal(MU, SIGMA)", ),
                 ( "nbinom", "Negative-Binomial(N, P)", ),
                 ( "norm", "Normal(MU, SIGMA)", ),
                 ( "pareto", "Pareto(XM, ALPHA)", ),
                 ( "poisson", "Poisson(MU)", ),
                 ( "randint", "Random-Integer or Discrete-Uniform(MIN, MAX)", ),
                 ( "t", "Students-T(DF)", ),
                 ( "uniform", "Uniform(MIN, MAX)", ),
                 ( "weibull_min", "Weibull(K, LAMBDA)", ),
               )

    lcdistname = str(distribname).lower()
    distrib = None
    if lcdistname == "beta":
        if distribparams == None:
            return ( ( "ALPHA", "first shape", ),
                     ( "BETA", "second shape", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Beta distribution")
        alpha = float(distribparams[0])
        beta = float(distribparams[1])
        if (alpha <= 0.0) or (beta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Beta distribution")
        distrib = scipy.stats.beta(alpha, beta)
    elif (lcdistname == "binom") or (lcdistname == "binomial"):
        if distribparams == None:
            return ( ( "N", "number of trials", ),
                     ( "P", "success probability in each trial", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Binomial distribution")
        nflt = float(distribparams[0])
        prob = float(distribparams[1])
        if (nflt < 0.0) or (prob < 0.0) or (prob > 1.0):
            raise ValueError("Invalid parameter(s) for the Binomial distribution")
        distrib = scipy.stats.binom(nflt, prob)
    elif lcdistname == "cauchy":
        if distribparams == None:
            return ( ( "M", "location (median)", ),
                     ( "GAMMA", "scale (half-width at half-maximum)", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Cauchy distribution")
        m = float(distribparams[0])
        gamma = float(distribparams[1])
        if gamma <= 0.0:
            raise ValueError("Invalid parameter for the Cauchy distribution")
        distrib = scipy.stats.cauchy(loc=m, scale=gamma)
    elif lcdistname == "chi":
        if distribparams == None:
            return ( ( "DF", "degrees of freedom", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Chi distribution")
        degfree = float(distribparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Chi distribution")
        distrib = scipy.stats.chi(degfree)
    elif (lcdistname == "chi2") or (lcdistname == "chi-square"):
        if distribparams == None:
            return ( ( "DF", "degrees of freedom", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Chi-Square distribution")
        degfree = float(distribparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Chi-Square distribution")
        distrib = scipy.stats.chi2(degfree)
    elif (lcdistname == "expon") or (lcdistname == "exponential"):
        if distribparams == None:
            return ( ( "LAMBDA", "rate (inverse scale)", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Exponential distribution")
        lambdaflt = float(distribparams[0])
        if lambdaflt <= 0.0:
            raise ValueError("Invalid parameter for the Exponential distribution")
        distrib = scipy.stats.expon(scale=(1.0/lambdaflt))
    elif (lcdistname == "exponweib") or (lcdistname == "exponentiated-weibull"):
        if distribparams == None:
            return ( ( "K", "Weibull shape", ),
                     ( "LAMBDA", "scale", ),
                     ( "ALPHA", "power shape", ), )
        if len(distribparams) != 3:
            raise ValueError("Three parameters expected for the Exponentiated-Weibull distribution")
        k =  float(distribparams[0])
        lambdaflt = float(distribparams[1])
        alpha = float(distribparams[2])
        if (k <= 0.0) or (lambdaflt <= 0.0) or (alpha <= 0):
            raise ValueError("Invalid parameter(s) for the Exponentiated-Weibull distribution")
        distrib = scipy.stats.exponweib(alpha, k, scale=lambdaflt)
    elif (lcdistname == "f") or (lcdistname == "fisher"):
        if distribparams == None:
            return ( ( "DFN", "numerator degrees of freedom", ),
                     ( "DFD", "denominator degrees of freedom", ), )
        if len(distribparams) != 2:
           raise ValueError("Two parameters expected for the F distribution")
        dfnum = float(distribparams[0])
        dfdenom = float(distribparams[1])
        if (dfnum <= 0.0) or (dfdenom <= 0.0):
           raise ValueError("Invalid parameter(s) for the F distribution")
        distrib = scipy.stats.f(dfnum, dfdenom)
    elif lcdistname == "gamma":
        if distribparams == None:
            return ( ( "ALPHA", "shape", ),
                     ( "THETA", "scale", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Gamma distribution")
        alpha = float(distribparams[0])
        theta = float(distribparams[1])
        if (alpha <= 0.0) or (theta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Gamma distribution")
        distrib = scipy.stats.gamma(alpha, scale=theta)
    elif (lcdistname == "geom") or (lcdistname == "geometric") or (lcdistname == "shifted-geometric"):
        if distribparams == None:
            return ( ( "P", "success probability", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Shifted-Geometric distribution")
        prob = float(distribparams[0])
        if (prob < 0.0) or (prob > 1.0):
            raise ValueError("Invalid parameter for the Shifted-Geometric distribution")
        distrib = scipy.stats.geom(prob)
    elif (lcdistname == "hypergeom") or (lcdistname == "hypergeometric"):
        if distribparams == None:
            return ( ( "NTOTAL", "total number of items", ),
                     ( "NGOOD", "total number of 'success' items", ),
                     ( "NDRAWN", "number of items selected", ), )
        if len(distribparams) != 3:
           raise ValueError("Three parameters expected for the Hypergeometric distribution")
        numtotal = float(distribparams[0])
        numgood = float(distribparams[1])
        numdrawn = float(distribparams[2])
        if (numtotal <= 0.0) or (numgood < 0.0) or (numdrawn < 0.0):
           raise ValueError("Invalid parameter(s) for the Hypergeometric distribution")
        distrib = scipy.stats.hypergeom(numtotal, numgood, numdrawn)
    elif (lcdistname == "invgamma") or (lcdistname == "inverse-gamma"):
        if distribparams == None:
            return ( ( "ALPHA", "shape", ),
                     ( "BETA", "scale", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Inverse-Gamma distribution")
        alpha = float(distribparams[0])
        beta = float(distribparams[1])
        if (alpha <= 0.0) or (beta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Inverse-Gamma distribution")
        distrib = scipy.stats.invgamma(alpha, scale=beta)
    elif lcdistname == "laplace":
        if distribparams == None:
            return ( ( "MU", "location (mean)", ),
                     ( "B", "scale", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Laplace distribution")
        mu = float(distribparams[0])
        b = float(distribparams[1])
        if b <= 0.0:
            raise ValueError("Invalid parameter for the Laplace distribution")
        distrib = scipy.stats.laplace(loc=mu, scale=b)
    elif (lcdistname == "lognorm") or (lcdistname == "log-normal"):
        if distribparams == None:
            return ( ( "MU", "log-scale (mean of the natural log of the distribution)", ),
                     ( "SIGMA", "shape (std. dev. of the natural log of the distribution)", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Log-Normal distribution")
        mu = math.exp(float(distribparams[0]))
        sigma = float(distribparams[1])
        if sigma <= 0.0:
            raise ValueError("Invalid parameter for the Log-Normal distribution")
        distrib = scipy.stats.lognorm(sigma, scale=mu)
    elif (lcdistname == "nbinom") or (lcdistname == "negative-binomial"):
        if distribparams == None:
            return ( ( "N", "number of successes to stop", ),
                     ( "P", "success probability in each trial", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Negative-Binomial distribution")
        numsuccess = float(distribparams[0])
        prob = float(distribparams[1])
        if (numsuccess < 1.0) or (prob <= 0.0) or (prob > 1.0):
            raise ValueError("Invalid parameter(s) for the Negative-Binomial distribution")
        distrib = scipy.stats.nbinom(numsuccess, prob)
    elif (lcdistname == "norm") or (lcdistname == "normal"):
        if distribparams == None:
            return ( ( "MU", "mean value", ),
                     ( "SIGMA", "standard deviation", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Normal distribution")
        mu = float(distribparams[0])
        sigma = float(distribparams[1])
        if sigma <= 0.0:
            raise ValueError("Invalid parameter for the Normal distribution")
        distrib = scipy.stats.norm(loc=mu, scale=sigma)
    elif lcdistname == "pareto":
        if distribparams == None:
            return ( ( "XM", "scale (minimum abscissa value)", ),
                     ( "ALPHA", "shape", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Pareto distribution")
        xm =  float(distribparams[0])
        alpha = float(distribparams[1])
        if (xm <= 0.0) or (alpha <= 0.0):
            raise ValueError("Invalid parameter(s) for the Pareto distribution")
        distrib = scipy.stats.pareto(alpha, scale=xm)
    elif lcdistname == "poisson":
        if distribparams == None:
            return ( ( "MU", "expected number of occurences", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Poisson distribution")
        mu = float(distribparams[0])
        if mu <= 0.0:
            raise ValueError("Invalid parameter for the Poisson distribution")
        distrib = scipy.stats.poisson(mu)
    elif (lcdistname == "randint") or (lcdistname == "random-integer") or (lcdistname == "discrete-uniform"):
        if distribparams == None:
            return ( ( "MIN", "minimum integer", ),
                     ( "MAX", "maximum integer (included)", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Random-Integer distribution")
        min = int(distribparams[0])
        max = int(distribparams[1])
        # randint takes int values, thus float values are truncated
        # this could lead to unexpected behavior (eg, one might expect
        # (0.9,10.1) to be treated as [1,11) but instead it becomes [0,10)
        minflt = float(distribparams[0])
        maxflt = float(distribparams[1])
        if (min >= max) or (min != minflt) or (max != maxflt):
            raise ValueError("Invalid parameters for the Random-Integer distribution")
        distrib = scipy.stats.randint(min, max + 1)
    elif (lcdistname == "t") or (lcdistname == "students-t"):
        if distribparams == None:
            return ( ( "DF", "degrees of freedom", ), )
        if len(distribparams) != 1:
            raise ValueError("One parameter expected for the Students-T distribution")
        degfree = float(distribparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Students-T distribution")
        distrib = scipy.stats.t(degfree)
    elif lcdistname == "uniform":
        if distribparams == None:
            return ( ( "MIN", "minimum", ),
                     ( "MAX", "maximum", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Uniform distribution")
        min = float(distribparams[0])
        max = float(distribparams[1])
        if min >= max:
            raise ValueError("Invalid parameters for the Uniform distribution")
        distrib = scipy.stats.uniform(loc=min, scale=(max - min))
    elif (lcdistname == "weibull_min") or (lcdistname == "weibull"):
        if distribparams == None:
            return ( ( "K", "shape", ),
                     ( "LAMBDA", "scale", ), )
        if len(distribparams) != 2:
            raise ValueError("Two parameters expected for the Weibull distribution")
        k =  float(distribparams[0])
        lambdaflt = float(distribparams[1])
        if (k <= 0.0) or (lambdaflt <= 0.0):
            raise ValueError("Invalid parameter(s) for the Weibull distribution")
        distrib = scipy.stats.weibull_min(k, scale=lambdaflt)
    else:
        raise ValueError("Unknown probability function %s" % str(distribname))
    if distrib == None:
        raise ValueError("Unexpected problem obtaining the probability distribution object")
    return distrib


def getfitparams(values, distribname, estparams):
    """
    Returns a tuple of "standard" parameters (including ordering) for a
    continuous probability distribution type named in distribname that
    best fits the distribution of data given in values (a 1-D array of
    data with no missing values).  Initial estimates for these "standard"
    parameters are given in estparams.
    """
    lcdistname = str(distribname).lower()
    fitparams = None
    if lcdistname == "beta":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Beta distribution")
        alpha = float(estparams[0])
        beta = float(estparams[1])
        if (alpha <= 0.0) or (beta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Beta distribution")
        fitparams = scipy.stats.beta.fit(values, a=alpha, b=beta, loc=0.0, scale=1.0)
        estparams = ( alpha, beta, 0.0, 1.0 )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif lcdistname == "cauchy":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Cauchy distribution")
        m = float(estparams[0])
        gamma = float(estparams[1])
        if gamma <= 0.0:
            raise ValueError("Invalid parameter for the Cauchy distribution")
        fitparams = scipy.stats.cauchy.fit(values, loc=m, scale=gamma)
        estparams = ( m, gamma, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif lcdistname == "chi":
        if len(estparams) != 1:
            raise ValueError("One parameter expected for the Chi distribution")
        degfree = float(estparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Chi distribution")
        fitparams = scipy.stats.chi.fit(values, df=degfree, loc=0.0, scale=1.0)
        estparams = ( degfree, 0.0, 1.0 )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif (lcdistname == "chi2") or (lcdistname == "chi-square"):
        if len(estparams) != 1:
            raise ValueError("One parameter expected for the Chi-Square distribution")
        degfree = float(estparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Chi-Square distribution")
        fitparams = scipy.stats.chi2.fit(values, df=degfree, loc=0.0, scale=1.0)
        estparams = ( degfree, 0.0, 1.0 )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif (lcdistname == "expon") or (lcdistname == "exponential"):
        if len(estparams) != 1:
            raise ValueError("One parameter expected for the Exponential distribution")
        lambdaflt = float(estparams[0])
        if lambdaflt <= 0.0:
            raise ValueError("Invalid parameter for the Exponential distribution")
        fitparams = scipy.stats.expon.fit(values, loc=0.0, scale=(1.0/lambdaflt))
        estparams = ( 0.0, 1.0/lambdaflt, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( 1.0 / fitparams[1], fitparams[0], )
    elif (lcdistname == "exponweib") or (lcdistname == "exponentiated-weibull"):
        if len(estparams) != 3:
            raise ValueError("Three parameters expected for the Exponentiated-Weibull distribution")
        k =  float(estparams[0])
        lambdaflt = float(estparams[1])
        alpha = float(estparams[2])
        if (k <= 0.0) or (lambdaflt <= 0.0) or (alpha <= 0):
            raise ValueError("Invalid parameter(s) for the Exponentiated-Weibull distribution")
        fitparams = scipy.stats.exponweib.fit(values, a=alpha, c=k, loc=0.0, scale=lambdaflt)
        estparams = ( alpha, k, 0.0, lambdaflt, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[1], fitparams[3], fitparams[0], fitparams[2], )
    elif (lcdistname == "f") or (lcdistname == "fisher"):
        if len(estparams) != 2:
           raise ValueError("Two parameters expected for the F distribution")
        dfnum = float(estparams[0])
        dfdenom = float(estparams[1])
        if (dfnum <= 0.0) or (dfdenom <= 0.0):
           raise ValueError("Invalid parameter(s) for the F distribution")
        fitparams = scipy.stats.f.fit(values, dfn=dfnum, dfd=dfdenom, loc=0.0, scale=1.0)
        estparams = ( dfnum, dfdenom, 0.0, 1.0, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif lcdistname == "gamma":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Gamma distribution")
        alpha = float(estparams[0])
        theta = float(estparams[1])
        if (alpha <= 0.0) or (theta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Gamma distribution")
        fitparams = scipy.stats.gamma.fit(values, a=alpha, loc=0.0, scale=theta)
        estparams = ( alpha, 0.0, theta, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[0], fitparams[2], fitparams[1], )
    elif (lcdistname == "invgamma") or (lcdistname == "inverse-gamma"):
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Inverse-Gamma distribution")
        alpha = float(estparams[0])
        beta = float(estparams[1])
        if (alpha <= 0.0) or (beta <= 0.0):
            raise ValueError("Invalid parameter(s) for the Inverse-Gamma distribution")
        fitparams = scipy.stats.invgamma.fit(values, a=alpha, loc=0.0, scale=beta)
        estparams = ( alpha, 0.0, beta, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[0], fitparams[2], fitparams[1], )
    elif lcdistname == "laplace":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Laplace distribution")
        mu = float(estparams[0])
        b = float(estparams[1])
        if b <= 0.0:
            raise ValueError("Invalid parameter for the Laplace distribution")
        fitparams = scipy.stats.laplace.fit(values, loc=mu, scale=b)
        estparams = ( mu, b, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif (lcdistname == "lognorm") or (lcdistname == "log-normal"):
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Log-Normal distribution")
        mu = math.exp(float(estparams[0]))
        sigma = float(estparams[1])
        if sigma <= 0.0:
            raise ValueError("Invalid parameter for the Log-Normal distribution")
        fitparams = scipy.stats.lognorm.fit(values, s=sigma, loc=0.0, scale=mu)
        estparams = ( sigma, 0.0, mu, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( math.log(fitparams[2]), fitparams[0], fitparams[1], )
    elif (lcdistname == "norm") or (lcdistname == "normal"):
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Normal distribution")
        mu = float(estparams[0])
        sigma = float(estparams[1])
        if sigma <= 0.0:
            raise ValueError("Invalid parameter for the Normal distribution")
        fitparams = scipy.stats.norm.fit(values, loc=mu, scale=sigma)
        estparams = ( mu, sigma, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif lcdistname == "pareto":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Pareto distribution")
        xm =  float(estparams[0])
        alpha = float(estparams[1])
        if (xm <= 0.0) or (alpha <= 0.0):
            raise ValueError("Invalid parameter(s) for the Pareto distribution")
        fitparams = scipy.stats.pareto.fit(values, b=alpha, loc=0.0, scale=xm)
        estparams = ( alpha, 0.0, xm, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[2], fitparams[0], fitparams[1], )
    elif (lcdistname == "t") or (lcdistname == "students-t"):
        if len(estparams) != 1:
            raise ValueError("One parameter expected for the Students-T distribution")
        degfree = float(estparams[0])
        if degfree <= 0.0:
            raise ValueError("Invalid parameter for the Students-T distribution")
        fitparams = scipy.stats.t.fit(values, df=degfree, loc=0.0, scale=1.0)
        estparams = ( degfree, 0.0, 1.0, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
    elif lcdistname == "uniform":
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Uniform distribution")
        min = float(estparams[0])
        max = float(estparams[1])
        if min >= max:
            raise ValueError("Invalid parameters for the Uniform distribution")
        fitparams = scipy.stats.uniform.fit(values, loc=min, scale=(max - min))
        estparams = ( min, (max - min), )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[0], fitparams[0] + fitparams[1], )
    elif (lcdistname == "weibull_min") or (lcdistname == "weibull"):
        if len(estparams) != 2:
            raise ValueError("Two parameters expected for the Weibull distribution")
        k =  float(estparams[0])
        lambdaflt = float(estparams[1])
        if (k <= 0.0) or (lambdaflt <= 0.0):
            raise ValueError("Invalid parameter(s) for the Weibull distribution")
        fitparams = scipy.stats.weibull_min.fit(values, c=k, loc=0.0, scale=lambdaflt)
        estparams = ( k, 0.0, lambdaflt, )
        if numpy.allclose(fitparams, estparams):
            raise ValueError("Fit parameters identical to estimate parameters")
        fitparams = ( fitparams[0], fitparams[2], fitparams[1], )
    else:
        raise ValueError("Unknown probability function %s" % str(distribname))
    if fitparams == None:
        raise ValueError("Unexpected problem parameterizing a distribution to fit given values")
    return fitparams


def getinitdict(distribname, funcname):
    """
    Returns a dictionary appropriate for the return value of ferret_init
    in a Ferret stats_<disribname>_<funcname> PyEF

    Arguments:
       distribname - name of the probability distribution
       funcname - name of the scipy.stats function
    """
    # generate a long function name from the scipy.stats function name
    if ( funcname == "cdf" ):
        funcaction = "calculate"
        funcreturn = "cumulative density function values"
    elif ( funcname == "isf" ):
        funcaction = "calculate"
        funcreturn = "inversion survival function values"
    elif ( funcname == "pdf" ):
        funcaction = "calculate"
        funcreturn = "probability distribution function values"
    elif ( funcname == "pmf" ):
        funcaction = "calculate"
        funcreturn = "probability mass function values"
    elif ( funcname == "ppf" ):
        funcaction = "calculate"
        funcreturn = "percent point function values"
    elif ( funcname == "sf" ):
        funcaction = "calculate"
        funcreturn = "survival function values"
    elif ( funcname == "rvs" ):
        funcaction = "assign"
        funcreturn = "random variates"
    else:
        raise ValueError("Unsupported scipy.stats function name '%s'" % funcname)
    # Get the distribution parameters information
    paramdescripts = getdistrib(distribname, None)
    numargs = len(paramdescripts) + 1
    if ( numargs == 2 ):
        # info for distributions with one parameter
        descript = "Returns (X=PTS,Y=%s) array of %s for %s prob. distrib." % \
                   (paramdescripts[0][0], funcreturn, distribname)
        axes = ( pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_DOES_NOT_EXIST,
                 pyferret.AXIS_DOES_NOT_EXIST, )
        argnames = ( "PTS", paramdescripts[0][0], )
        argdescripts = ( "Point(s) at which to %s the %s" % (funcaction, funcreturn),
                         "Parameter(s) defining the %s" % paramdescripts[0][1], )
        argtypes = ( pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, )
        influences = ( ( False, False, False, False, ),
                       ( False, False, False, False, ), )
    elif (numargs == 3):
        # info for distributions with two parameters
        descript = "Returns (X=PTS,Y=%s,Z=%s) array of %s for %s prob. distrib." % \
                   (paramdescripts[0][0], paramdescripts[1][0], funcreturn, distribname)
        axes = ( pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_DOES_NOT_EXIST, )
        argnames = ( "PTS", paramdescripts[0][0], paramdescripts[1][0], )
        argdescripts = ( "Point(s) at which to %s the %s" % (funcaction, funcreturn),
                         "Parameter(s) defining the %s" % paramdescripts[0][1],
                         "Parameter(s) defining the %s" % paramdescripts[1][1], )
        argtypes = ( pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, )
        influences = ( ( False, False, False, False, ),
                       ( False, False, False, False, ),
                       ( False, False, False, False, ), )
    elif (numargs == 4):
        # info for distributions with three parameters
        descript = "Returns (X=PTS,Y=%s,Z=%s,T=%s) array of %s for %s prob. distrib." % \
                   (paramdescripts[0][0], paramdescripts[1][0], paramdescripts[2][0],
                    funcreturn, distribname)
        axes = ( pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM,
                 pyferret.AXIS_CUSTOM, )
        argnames = ( "PTS", paramdescripts[0][0], paramdescripts[1][0], paramdescripts[2][0], )
        argdescripts = ( "Point(s) at which to %s the %s" % (funcaction, funcreturn),
                         "Parameter(s) defining the %s" % paramdescripts[0][1],
                         "Parameter(s) defining the %s" % paramdescripts[1][1],
                         "Parameter(s) defining the %s" % paramdescripts[2][1], )
        argtypes = ( pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, pyferret.FLOAT_ARG, )
        influences = ( ( False, False, False, False, ),
                       ( False, False, False, False, ),
                       ( False, False, False, False, ),
                       ( False, False, False, False, ), )
    else:
        raise ValueError("Unexpected number of arguments: %d" % numargs)
    # Create and return the dictionary
    return { "numargs": numargs,
             "descript": descript,
             "axes": axes,
             "argnames": argnames,
             "argdescripts": argdescripts,
             "argtypes": argtypes,
             "influences": influences, }


def getcustomaxisvals(id, distribname):
    """
    Returns a 4-tuple of custom axis values appropriate for the return value
    of ferret_custom_axis in a Ferret stats_<disribname>_<funcname> PyEF

    Arguments:
       distribname - name of the probability distribution
    """
    # Get the distribution parameters information
    paramdescripts = getdistrib(distribname, None)
    numargs = len(paramdescripts) + 1
    namelist = [ "PTS", None, None, None ]
    for k in xrange(1, numargs):
        namelist[k] = paramdescripts[k-1][0]
    argvals = ( pyferret.ARG1, pyferret.ARG2, pyferret.ARG3, pyferret.ARG4 )
    axisvals = ( pyferret.X_AXIS, pyferret.Y_AXIS, pyferret.Z_AXIS, pyferret.T_AXIS )
    customaxisvals = [ None, None, None, None ]
    for k in xrange(numargs):
        arglen = 1
        for axis in axisvals:
            axis_info = pyferret.get_axis_info(id, argvals[k], axis)
            num = axis_info.get("size", -1)
            if num > 0:
                arglen *= num
        # if all axes have undefined lengths, assume it is a single value
        customaxisvals[k] = ( 1, arglen, 1, namelist[k], False )
    return customaxisvals


def getdistribfunc(distrib, funcname):
    """
    Returns the distrib.funcname function for recognized funcnames
    """
    if ( funcname == "cdf" ):
        return distrib.cdf
    elif ( funcname == "isf" ):
        return distrib.isf
    elif ( funcname == "pdf" ):
        return distrib.pdf
    elif ( funcname == "pmf" ):
        return distrib.pmf
    elif ( funcname == "ppf" ):
        return distrib.ppf
    elif ( funcname == "sf" ):
        return distrib.sf
    elif ( funcname == "rvs" ):
        return distrib.rvs
    else:
        raise ValueError("Unsupported scipy.stats function name '%s'" % funcname)


def assignresultsarray(distribname, funcname, result, resbdf, inputs, inpbdfs):
    """
    Assigns result with the funcname function values for the distribname
    probability distributions defined by parameters in inputs[1:]
    using the abscissa or template values given in inputs[0].
    """
    ptvals = inputs[0].reshape(-1, order='F')
    badmask = ( numpy.fabs(ptvals - inpbdfs[0]) < 1.0E-5 )
    badmask = numpy.logical_or(badmask, numpy.isnan(ptvals))
    goodmask = numpy.logical_not(badmask)
    numparams = len(inputs) - 1
    if numparams == 1:
        par1vals = inputs[1].reshape(-1, order='F')
        # check that result is the required shape
        expshape = ( len(ptvals), len(par1vals), 1, 1 )
        if result.shape != expshape:
            raise ValueError("Results array size mismatch; expected: %s; found %s" % \
                             (str(expshape), str(result.shape)))
        for j in xrange(expshape[1]):
            try:
                distrib = getdistrib(distribname, ( par1vals[j], ))
                distribfunc = getdistribfunc(distrib, funcname)
                if funcname == "rvs":
                    # goodmask is one-D
                    result[goodmask, j, 0, 0] = getdistribfunc(distrib, funcname)(len(goodmask))
                else:
                    result[goodmask, j, 0, 0] = getdistribfunc(distrib, funcname)(ptvals[goodmask])
                result[badmask, j, 0, 0] = resbdf
            except ValueError, msg:
                # print msg
                result[:, j, 0, 0] = resbdf
    elif numparams == 2:
        par1vals = inputs[1].reshape(-1, order='F')
        par2vals = inputs[2].reshape(-1, order='F')
        # check that result is the required shape
        expshape = ( len(ptvals), len(par1vals), len(par2vals), 1 )
        if result.shape != expshape:
            raise ValueError("Results array size mismatch; expected: %s; found %s" % \
                             (str(expshape), str(result.shape)))
        for j in xrange(expshape[1]):
            for k in xrange(expshape[2]):
                try:
                    distrib = getdistrib(distribname, ( par1vals[j], par2vals[k], ))
                    distribfunc = getdistribfunc(distrib, funcname)
                    if funcname == "rvs":
                        # goodmask is one-D
                        result[goodmask, j, k, 0] = getdistribfunc(distrib, funcname)(len(goodmask))
                    else:
                        result[goodmask, j, k, 0] = getdistribfunc(distrib, funcname)(ptvals[goodmask])
                    result[badmask, j, k, 0] = resbdf
                except ValueError, msg:
                    # print msg
                    result[:, j, k, 0] = resbdf
    elif numparams == 3:
        par1vals = inputs[1].reshape(-1, order='F')
        par2vals = inputs[2].reshape(-1, order='F')
        par3vals = inputs[3].reshape(-1, order='F')
        # check that result is the required shape
        expshape = ( len(ptvals), len(par1vals), len(par2vals), len(par3vals) )
        if result.shape != expshape:
            raise ValueError("Results array size mismatch; expected: %s; found %s" % \
                             (str(expshape), str(result.shape)))
        for j in xrange(expshape[1]):
            for k in xrange(expshape[2]):
                for q in xrange(expshape[3]):
                    try:
                        distrib = getdistrib(distribname, ( par1vals[j], par2vals[k], par3vals[q], ))
                        distribfunc = getdistribfunc(distrib, funcname)
                        if funcname == "rvs":
                            # goodmask is one-D
                            result[goodmask, j, k, q] = getdistribfunc(distrib, funcname)(len(goodmask))
                        else:
                            result[goodmask, j, k, q] = getdistribfunc(distrib, funcname)(ptvals[goodmask])
                        result[badmask, j, k, q] = resbdf
                    except ValueError, msg:
                        # print msg
                        result[:, j, k, q] = resbdf
    else:
        raise ValueError("Unexpected number of parameters: %d" % numparams)


#
# The rest of this is just for testing this module at the command line
#
if __name__ == "__main__":
    # Test the distribution scipy name and parameters given to getdistrib
    # give the expected distribution.  (Primarily that the parameters
    # are interpreted and assigned correctly.)  Testing of the long names
    # is performed by the stats_helper.py script.


    # Number of distributions
    distdescripts = getdistrib(None, None)
    if len(distdescripts) != NUM_DISTRIBS:
        raise ValueError("number of distribution description pairs: expected %d; found %d" % \
                         (NUM_DISTRIBS, len(distdescripts)))


    # Beta distribution
    distname = "beta"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d" % \
                         (distname, len(descript)))
    alpha = 1.5
    beta = 2.75
    distparms = [ alpha, beta ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( alpha / (alpha + beta),
                      alpha * beta / ((alpha + beta)**2 * (alpha + beta + 1.0)),
                      2.0 * (beta - alpha) / (2.0 + alpha + beta) *
                          math.sqrt((1.0 + alpha + beta) / (alpha * beta)),
                      6.0 * (alpha**3 + alpha**2 * (1.0 - 2.0 * beta) + \
                          beta**2 * (1.0 + beta) - 2.0 * alpha * beta * (2.0 + beta)) / \
                          (alpha * beta * (alpha + beta + 2.0) * (alpha + beta + 3.0)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc and scale to the expected params
    distparms.append(0.0)
    distparms.append(1.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, alpha, beta, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Binomial distribution
    distname = "binom"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    ntrials = 20.0
    prob = 0.25
    distparms = [ ntrials, prob ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( ntrials * prob,
                      ntrials * prob * (1.0 - prob),
                      (1.0 - 2.0 * prob) / math.sqrt(ntrials * prob * (1.0 - prob)),
                      (1.0 - 6.0 * prob * (1.0 - prob)) / (ntrials * prob * (1.0 - prob)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    # no binom.fit function

    del descript, ntrials, prob, distparms, distf, foundstats, expectedstats
    print "%s: PASS" % distname


    # Cauchy distribution
    distname = "cauchy"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    m = 5.0
    gamma = 2.0
    distparms = [ m, gamma ]
    distf = getdistrib(distname, distparms)
    # mean, variance, skew, kurtosis undefined; instead check some pdf values
    xvals = numpy.arange(0.0, 10.1, 0.5)
    foundpdfs = distf.pdf(xvals)
    expectedpdfs = (gamma / numpy.pi) / ((xvals - m)**2 + gamma**2)
    if not numpy.allclose(foundpdfs, expectedpdfs):
        print "%s: FAIL" % distname
        raise ValueError("pdfs(0.0:10.1:0.5) of %s(%#.1f,%#.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedpdfs), str(foundpdfs)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, m, gamma, distparms, distf, xvals, foundpdfs, expectedpdfs, sample, fitparms
    print "%s: PASS" % distname


    # Chi distribution
    distname = "chi"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))
    degfree = 10
    distparms = [ degfree ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    mean = math.sqrt(2.0) * scipy.special.gamma(0.5 * (degfree + 1.0)) / \
                            scipy.special.gamma(0.5 * degfree)
    variance = degfree - mean**2
    stdev = math.sqrt(variance)
    skew = mean * (1.0 - 2.0 * variance) / stdev**3
    expectedstats = ( mean,
                      variance,
                      skew,
                      2.0 * (1.0 - mean * stdev * skew - variance) / variance,
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%d.0): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc and scale to the expected params
    distparms.append(0.0)
    distparms.append(1.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.4, atol=0.4):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, degfree, distparms, distf, foundstats, mean, variance, stdev, skew, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Chi-squared distribution
    distname = "chi2"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))
    degfreestr = "10"
    distparms = [ degfreestr ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    degfree = float(degfreestr)
    expectedstats = ( degfree,
                      2.0 * degfree,
                      math.sqrt(8.0 / degfree),
                      12.0 / degfree,
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%s): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc and scale to the expected params
    distparms = [ degfree, 0.0, 1.0 ]
    if not numpy.allclose(fitparms, distparms, rtol=0.4, atol=0.4):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, degfreestr, distparms, distf, foundstats, degfree, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Exponential distribution
    distname = "expon"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))

    lambdaflt = 11.0
    distparms = [ lambdaflt ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( 1.0 / lambdaflt, 1.0 / lambdaflt**2, 2.0, 6.0 )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, lambdaflt, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname

    # Exponentiated Weibull distribution
    distname = "exponweib"
    descript = getdistrib(distname, None)
    if len(descript) != 3:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 3; found %d:" % \
                         (distname, len(descript)))
    k = 3.0
    lambdaflt = 5.0
    alpha = 2.5
    distparms = [ k, lambdaflt, alpha ]
    distf = getdistrib(distname, distparms)
    # don't know the formula for the mean, variance, skew, kurtosis
    # instead check some cdf values
    xvals = numpy.arange(0.0, 10.1, 0.5)
    foundcdfs = distf.cdf(xvals)
    expectedcdfs = numpy.power(1.0 - numpy.exp(-1.0 * numpy.power(xvals / lambdaflt, k)), alpha)
    if not numpy.allclose(foundcdfs, expectedcdfs):
        print "%s: FAIL" % distname
        raise ValueError("cdfs(0.0:10.1:0.5) of %s(%#.1f,%#.1f%#.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], distparms[2], str(expectedcdfs), str(foundcdfs)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, k, lambdaflt, alpha, distparms, distf, xvals, foundcdfs, expectedcdfs, sample, fitparms
    print "%s: PASS" % distname


    # F distribution
    distname = "f"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    dofn = 7.0
    dofd = 11.0   # needs to be larger than 8.0 for kurtosis formula
    distparms = [ dofn, dofd ]
    distf = getdistrib(distname, distparms)
    # foundstats = distf.stats("mvsk")
    foundstats = distf.stats("mv")
    expectedstats = ( dofd / (dofd - 2.0),
                      2.0 * dofd**2 * (dofn + dofd - 2.0) / \
                          (dofn * (dofd - 2.0)**2 * (dofd - 4.0)),
                      # ((2.0 * dofn + dofd - 2.0) / (dofd - 6.0)) * \
                      #     math.sqrt(8.0 * (dofd - 4.0) / (dofn * (dofn + dofd - 2.0))),
                      # 12.0 * (20.0 * dofd - 8.0 * dofd**2 + dofd**3 + 44.0 * dofn - 32.0 * dofn * dofd + \
                      #     5.0 * dofd**2 * dofn - 22.0 * dofn**2 + 5.0 * dofd * dofn**2 - 16.0) / \
                      #     (dofn * (dofd - 6.0) * (dofd - 8.0) * (dofn + dofd - 2)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        # raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
        print "%s: FAIL" % distname
        raise ValueError("(mean, var) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))
    # since skew and kurtosis is not coming out as expected, check some pdf values
    xvals = numpy.arange(0.0, 10.1, 0.5)
    foundpdfs = distf.pdf(xvals)
    factor = scipy.special.gamma(0.5 * (dofn + dofd)) / \
             (scipy.special.gamma(0.5 * dofn) * scipy.special.gamma(0.5 *dofd))
    factor *= math.pow(dofn, 0.5 * dofn) * math.pow(dofd, 0.5 * dofd)
    expectedpdfs = factor * numpy.power(xvals, 0.5 * dofn - 1.0) / \
                   numpy.power(dofd + dofn * xvals, 0.5 * (dofn + dofd))
    if not numpy.allclose(foundpdfs, expectedpdfs):
        print "%s: FAIL" % distname
        raise ValueError("pdfs(0.0:10.1:0.5) of %s(%#.1f,%#.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedpdfs), str(foundpdfs)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc and scale to the expected params
    distparms.append(0.0)
    distparms.append(1.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.2, atol=0.4):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, dofn, dofd, distparms, distf, foundstats, expectedstats
    del xvals, foundpdfs, factor, expectedpdfs, sample, fitparms
    print "%s: PASS" % distname


    # Gamma distribution
    distname = "gamma"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    alpha = 5.0
    theta = 3.0
    distparms = [ alpha, theta ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( alpha * theta, alpha * theta**2, 2.0 / math.sqrt(alpha), 6.0 / alpha )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, alpha, theta, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Geometric distribution
    distname = "geom"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))
    prob = 0.25
    distparms = [ prob ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( 1.0 / prob,
                     (1.0 - prob) / prob**2,
                     (2.0 - prob) / math.sqrt(1.0 - prob),
                     6.0 + prob**2 / (1.0 - prob),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    # no geom.fit function

    del descript, prob, distparms, distf, foundstats, expectedstats
    print "%s: PASS" % distname


    # Hypergeometric distribution
    distname = "hypergeom"
    descript = getdistrib(distname, None)
    if len(descript) != 3:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 3; found %d:" % \
                         (distname, len(descript)))
    numtotal = 29.0
    numgood = 13.0
    numdrawn = 17.0
    distparms = [ numtotal, numgood, numdrawn ]
    distf = getdistrib(distname, distparms)
    # foundstats = distf.stats("mvsk")
    foundstats = distf.stats("mvs")
    expectedstats = ( numdrawn * numgood / numtotal,
                      numdrawn * numgood * (numtotal - numdrawn) * (numtotal - numgood) / \
                          (numtotal**2 * (numtotal - 1.0)),
                      math.sqrt(numtotal - 1.0) * (numtotal - 2.0 * numdrawn) * (numtotal - 2.0 * numgood) / \
                          (math.sqrt(numdrawn * numgood * (numtotal - numdrawn) * (numtotal - numgood)) * \
                              (numtotal - 2.0)),
                      # (numtotal**2 * (numtotal - 1.0) / \
                      #         (numdrawn * (numtotal - 2.0) * (numtotal - 3.0) * (numtotal - numdrawn))) * \
                      #     ((numtotal * (numtotal + 1.0) - 6.0 * numtotal * (numtotal - numdrawn)) / \
                      #      (numgood * (numtotal - numgood)) + \
                      #      3.0 * numdrawn * (numtotal - numdrawn) * (numtotal + 6.0) / numtotal**2 - 6.0),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], distparms[2], str(expectedstats), str(foundstats)))

    # no hypergeom.fit function

    del descript, numtotal, numgood, numdrawn, distparms, distf, foundstats, expectedstats
    print "%s: PASS" % distname


    # Inverse-Gamma distribution
    distname = "invgamma"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    alpha = 7.0  # must be > 4 for the kurtosis formula
    beta = 3.0
    distparms = [ alpha, beta ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( beta / (alpha - 1.0),
                      beta**2 / ((alpha - 1.0)**2 * (alpha - 2.0)),
                      4.0 * math.sqrt(alpha - 2.0) / (alpha - 3.0),
                      (30.0 * alpha - 66.0)/ ((alpha - 3.0) * (alpha - 4.0)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    # get rid of any major outliers - seems to give invgamma.fit problems
    sample = sample[ numpy.logical_and((sample > 0.1), (sample < 2.0)) ]
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.2, atol=0.4):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, alpha, beta, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Laplace distribution
    distname = "laplace"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    mu = 5.0
    b = 3.0
    distparms = [ mu, b ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( mu, 2.0 * b**2, 0.0, 3.0 )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, mu, b, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Log-normal distribution
    distname = "lognorm"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    mu = 0.8
    sigma = 0.5
    distparms = [ mu, sigma ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( math.exp(mu + 0.5 * sigma**2),
                      math.exp(2.0 * mu + sigma**2) * (math.exp(sigma**2) - 1.0),
                      (2.0 + math.exp(sigma**2)) * math.sqrt(math.exp(sigma**2) - 1.0),
                      math.exp(4.0 * sigma**2) + 2.0 * math.exp(3.0 * sigma**2) + 3.0 * math.exp(2.0 * sigma**2) - 6,
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, mu, sigma, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Negative-binomial distribution
    distname = "nbinom"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    numsuccess = 5.0
    prob = 0.25
    distparms = [ numsuccess, prob ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( numsuccess * (1.0 - prob) / prob,
                      numsuccess * (1.0 - prob) / prob**2,
                      (2.0 - prob) / math.sqrt(numsuccess * (1.0 - prob)),
                      (prob**2 - 6.0 * prob + 6.0) / (numsuccess * (1.0 - prob)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    # no nbinom.fit function

    del descript, numsuccess, prob, distparms, distf, foundstats, expectedstats
    print "%s: PASS" % distname


    # Normal distribution
    distname = "norm"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    mu = 5.0
    sigma = 3.0
    distparms = numpy.array([ mu, sigma ], dtype=numpy.float32)
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( mu, sigma**2, 0.0, 0.0 )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, mu, sigma, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Pareto distribution
    distname = "pareto"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    xm = 3.0
    alpha = 5.0  # must be larger than 4 for kurtosis formula
    distparms = [ xm, alpha ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( alpha * xm / (alpha - 1.0),
                      xm**2 * alpha / ((alpha - 1.0)**2 * (alpha - 2.0)),
                      2.0 * ((alpha + 1.0) / (alpha - 3.0)) * math.sqrt((alpha - 2.0) / alpha),
                      6.0 * (alpha**3 + alpha**2 - 6.0 * alpha - 2.0) / \
                          (alpha * (alpha - 3.0) * (alpha - 4.0)),
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, xm, alpha, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname


    # Poisson distribution
    distname = "poisson"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))
    mu = 7.0
    distparms = [ mu ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( mu, mu, 1.0 / math.sqrt(mu), 1.0 / mu )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    # no poisson.fit function

    del descript, mu, distparms, distf, foundstats, expectedstats
    print "%s: PASS" % distname


    # Random Integer (Discrete Uniform) distribution
    distname = "randint"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    a = -5.0
    b = 13.0
    distparms = [ a, b ]
    distf = getdistrib(distname, distparms)
    # foundstats = distf.stats("mvsk")
    foundstats = distf.stats("mvs")
    n = b - a + 1.0
    # expectedstats = ( 0.5 * (a + b), (n**2 - 1.0) / 12.0, 0.0, -6.0 * (n**2 + 1) / (5.0 * (n**2 - 1)) )
    expectedstats = ( 0.5 * (a + b), (n**2 - 1.0) / 12.0, 0.0, )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        # raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
        raise ValueError("(mean, var, skew) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))
    xvals = numpy.arange(a - 1.0, b + 1.1, 1.0)
    expectedpmfs = numpy.ones((n+2,), dtype=float) / n
    expectedpmfs[0] = 0.0
    expectedpmfs[n+1] = 0.0
    foundpmfs = distf.pmf(xvals)
    if not numpy.allclose(foundpmfs, expectedpmfs):
        print "%s: FAIL" % distname
        raise ValueError("pmfs(%.1f:%.1f:1.0) of %s(%.1f, %.1f): expected %s; found %s" % \
              (a - 1.0, b + 1.1, distname, distparms[0], distparms[1], str(expectedpmfs), str(foundpmfs)))

    # no randint.fit function

    del descript, a, b, distparms, distf, foundstats, n, expectedstats, xvals, expectedpmfs, foundpmfs
    print "%s: PASS" % distname


    # Student's-t distribution
    distname = "t"
    descript = getdistrib(distname, None)
    if len(descript) != 1:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 1; found %d:" % \
                         (distname, len(descript)))
    degfree = 11.0
    distparms = [ degfree ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( 0.0, degfree / (degfree - 2.0), 0.0, 6.0 / (degfree - 4.0) )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f): expected %s; found %s" % \
                          (distname, distparms[0], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc and scale to the expected params
    distparms.append(0.0)
    distparms.append(1.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, degfree, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname

    # Uniform distribution
    distname = "uniform"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    a = -5.0
    b = 13.0
    distparms = [ a, b ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    expectedstats = ( 0.5 * (a + b), (b - a)**2 / 12.0, 0.0, -6.0 / 5.0 )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, a, b, distparms, distf, foundstats, expectedstats, sample, fitparms
    print "%s: PASS" % distname

    # Weibull distribution
    distname = "weibull_min"
    descript = getdistrib(distname, None)
    if len(descript) != 2:
        print "%s: FAIL" % distname
        raise ValueError("number of parameter description pairs for %s: expected 2; found %d:" % \
                         (distname, len(descript)))
    k = 3.0
    lambdaflt = 5.0
    distparms = [ k, lambdaflt ]
    distf = getdistrib(distname, distparms)
    foundstats = distf.stats("mvsk")
    gam1 = scipy.special.gamma(1.0 + 1.0 / k)
    gam2 = scipy.special.gamma(1.0 + 2.0 / k)
    gam3 = scipy.special.gamma(1.0 + 3.0 / k)
    gam4 = scipy.special.gamma(1.0 + 4.0 / k)
    mu = lambdaflt * gam1
    sigma = lambdaflt * math.sqrt(gam2 - gam1**2)
    expectedstats = ( mu,
                      sigma**2,
                      (lambdaflt**3 * gam3  - 3.0 * mu * sigma**2 - mu**3) / sigma**3,
                      (gam4 - 4.0 * gam1 * gam3 - 3.0 * gam2**2 + 12.0 * gam1**2 * gam2 - 6.0 * gam1**4) / \
                      (gam2 - gam1**2)**2,
                    )
    if not numpy.allclose(foundstats, expectedstats):
        print "%s: FAIL" % distname
        raise ValueError("(mean, var, skew, kurtosis) of %s(%.1f, %.1f): expected %s; found %s" % \
                          (distname, distparms[0], distparms[1], str(expectedstats), str(foundstats)))

    sample = distf.rvs(25000)
    fitparms = getfitparams(sample, distname, distparms)
    # append the default loc to the expected params
    distparms.append(0.0)
    if not numpy.allclose(fitparms, distparms, rtol=0.1, atol=0.2):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s: expected %s; found %s" % \
                                      (distname, str(distparms), str(fitparms)))
    if numpy.allclose(fitparms, distparms):
        print "%s: FAIL" % distname
        raise ValueError("fitparams of %s identical to input params" % distname)

    del descript, k, lambdaflt, distparms, distf, foundstats
    del gam1, gam2, gam3, gam4, mu, sigma, expectedstats, sample, fitparms
    print "%s: PASS" % distname

    # All successful
    print "Success"
