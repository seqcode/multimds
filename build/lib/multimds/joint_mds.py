"""
Jointly perform multi-dimensional Scaling (MDS) on two datasets
"""
# original author: Nelle Varoquaux <nelle.varoquaux@gmail.com>
# modified by: Lila Rieber <lur159@psu.edu>
# License: BSD

import numpy as np
import sys
import warnings

from sklearn.base import BaseEstimator
from sklearn.metrics import euclidean_distances
from sklearn.utils import check_random_state, check_array, check_symmetric
from sklearn.externals.joblib import Parallel
from sklearn.externals.joblib import delayed
from sklearn.isotonic import IsotonicRegression

def squared_dist(x1, x2):
     """Computes squared Euclidean distance between coordinate x1 and coordinate x2"""
     return sum([(i1 - i2)**2 for i1, i2 in zip(x1, x2)]) 

def ssd(X1, X2):
     """Computes sum of squared distances between coordinates X1 and coordinates X2"""
     return sum([squared_dist(x1, x2) for x1, x2 in zip(X1, X2)])

def moore_penrose(V):
     """Computes Moore-Penrose inverse of matrix V"""
     n = len(V)
     return np.linalg.inv(V + np.ones((n,n))) - n**-2 * np.ones((n,n))

def initialize(dissimilarities, random_state, init, n_samples, n_components):
    random_state = check_random_state(random_state)
    sim_flat = ((1 - np.tri(n_samples)) * dissimilarities).ravel()
    sim_flat_w = sim_flat[sim_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = random_state.rand(n_samples * n_components)
        X = X.reshape((n_samples, n_components))
    else:
        n_components = init.shape[1]
        if n_samples != init.shape[0]:
            raise ValueError("init matrix should be of shape (%d, %d)" %
                             (n_samples, n_components))
        X = init

    return X, sim_flat, sim_flat_w

def nonmetric_disparities(dis, sim_flat, n_samples):
    dis_flat = dis.ravel()
    # dissimilarities with 0 are considered as missing values
    dis_flat_w = dis_flat[sim_flat != 0]

    # Compute the disparities using a monotonic regression
    disparities_flat = ir.fit_transform(sim_flat_w, dis_flat_w)
    disparities = dis_flat.copy()
    disparities[sim_flat != 0] = disparities_flat
    disparities = disparities.reshape((n_samples, n_samples))
    disparities *= np.sqrt((n_samples * (n_samples - 1) / 2) /
                                   (disparities ** 2).sum())

    return disparities

def guttman(X1, X2, disparities, inv_V, V2, dis):
    # avoid division by 0
    dis[dis == 0] = 1e-5

    # B: error between distance matrix and embedding
    ratio = disparities / dis
    B = - ratio
    B[np.arange(len(B)), np.arange(len(B))] += ratio.sum(axis=1)
    
    return np.dot(inv_V, (np.dot(B, X1) + np.dot(V2, X2)))

def _smacof_single(dissimilarities1, dissimilarities2, p, weights1=None, weights2=None, metric=True, n_components=2, 
                   init1=None, init2=None, max_iter=300, verbose=0, eps=1e-3, 
                   random_state1=None, random_state2=None):
    """
    Computes multidimensional scaling using SMACOF algorithm

    Parameters
    ----------
    dissimilarities : ndarray, shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : boolean, optional, default: True
        Compute metric or nonmetric SMACOF algorithm.

    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : ndarray, shape (n_samples, n_components), optional, default: None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    random_state : integer or numpy.RandomState, optional
        The generator used to initialize the centers. If an integer is
        given, it fixes the seed. Defaults to the global numpy random
        number generator.

    Returns
    -------
    X : ndarray, shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).

    n_iter : int
        The number of iterations corresponding to the best stress.
    """
    dissimilarities1 = check_symmetric(dissimilarities1, raise_exception=True)
    dissimilarities2 = check_symmetric(dissimilarities2, raise_exception=True)

    if dissimilarities1.shape != dissimilarities2.shape:
         print("Error. Distance matrices have different shapes.")
         sys.exit("Error. Distance matrices have different shapes.")

    n_samples = dissimilarities1.shape[0]

    X1, sim_flat1, sim_flat_w1 = initialize(dissimilarities1, random_state1, 
    init1, n_samples, n_components)
    X2, sim_flat2, sim_flat_w2 = initialize(dissimilarities2, random_state2, 
    init2, n_samples, n_components)  

    #Default: equal weights
    if weights1 is None:
        weights1 = np.ones((n_samples, n_samples))
    if weights2 is None:
        weights2 = np.ones(n_samples)

    # Disparity-specific weights (V in Borg)
    V1 = np.zeros((n_samples,n_samples))
    for i in range(n_samples):
        diagonal = 0
        for j in range(n_samples):
            V1[i,j] = -weights1[i,j]
            diagonal += weights1[i,j]
        V1[i,i] = diagonal

    # Locus-specific weights
    V2 = np.zeros((n_samples,n_samples))
    for i, weight in enumerate(weights2):
         V2[i,i] = weight * p * n_samples

    inv_V = moore_penrose(V1+V2)

    old_stress = None
    ir = IsotonicRegression()
    for it in range(max_iter):
        # Compute distance and monotonic regression
        dis1 = euclidean_distances(X1)
        dis2 = euclidean_distances(X2)

        if metric:
            disparities1 = dissimilarities1
            disparities2 = dissimilarities2
        else:
            disparities1 = nonmetric_disparities1(dis1, sim_flat1, n_samples)
            disparities2 = nonmetric_disparities2(dis2, sim_flat2, n_samples)

        # Compute stress
        stress = ((dis1.ravel() - disparities1.ravel()) ** 2).sum() + ((dis2.ravel() - disparities2.ravel()) ** 2).sum() + n_samples * p * ssd(X1, X2)   #multiply by n_samples to make ssd term comparable in magnitude to embedding error terms

        # Update X1 using the Guttman transform
        X1 = guttman(X1, X2, disparities1, inv_V, V2, dis1)

        # Update X2 using the Guttman transform
        X2 = guttman(X2, X1, disparities2, inv_V, V2, dis2)
 
        # Test stress
        dis1 = np.sqrt((X1 ** 2).sum(axis=1)).sum()
        dis2 = np.sqrt((X2 ** 2).sum(axis=1)).sum()
        dis = np.mean((dis1, dis2))
        if verbose >= 2:
            print('it: %d, stress %s' % (it, stress))
        if old_stress is not None:
            if np.abs(old_stress - stress / dis) < eps:
                if verbose:
                    print('breaking at iteration %d with stress %s' % (it,
                                                                       stress))
                break
        old_stress = stress / dis

    return X1, X2, stress, it + 1


def smacof(dissimilarities1, dissimilarities2, p, weights1, weights2, metric=True, n_components=2, init1=None, init2=None, 
           n_init=8, n_jobs=1, max_iter=300, verbose=0, eps=1e-3, random_state1=None, random_state2=None,
           return_n_iter=False):
    """
    Computes multidimensional scaling using the SMACOF algorithm.

    The SMACOF (Scaling by MAjorizing a COmplicated Function) algorithm is a
    multidimensional scaling algorithm which minimizes an objective function
    (the *stress*) using a majorization technique. Stress majorization, also
    known as the Guttman Transform, guarantees a monotone convergence of
    stress, and is more powerful than traditional techniques such as gradient
    descent.

    The SMACOF algorithm for metric MDS can summarized by the following steps:

    1. Set an initial start configuration, randomly or not.
    2. Compute the stress
    3. Compute the Guttman Transform
    4. Iterate 2 and 3 until convergence.

    The nonmetric algorithm adds a monotonic regression step before computing
    the stress.

    Parameters
    ----------
    dissimilarities : ndarray, shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : boolean, optional, default: True
        Compute metric or nonmetric SMACOF algorithm.

    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : ndarray, shape (n_samples, n_components), optional, default: None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    n_init : int, optional, default: 8
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress. If ``init`` is
        provided, this option is overridden and a single run is performed.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For ``n_jobs`` below -1,
        (``n_cpus + 1 + n_jobs``) are used. Thus for ``n_jobs = -2``, all CPUs
        but one are used.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize the centers. If an integer is given,
        it fixes the seed. Defaults to the global numpy random number
        generator.

    return_n_iter : bool, optional, default: False
        Whether or not to return the number of iterations.

    Returns
    -------
    X : ndarray, shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).

    n_iter : int
        The number of iterations corresponding to the best stress. Returned
        only if ``return_n_iter`` is set to ``True``.

    Notes
    -----
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)
    """

    if p < 0:
        sys.exit('Error. Penalty must be non-negative.')

    dissimilarities1 = check_array(dissimilarities1)
    dissimilarities2 = check_array(dissimilarities2)
    random_state1 = check_random_state(random_state1)
    random_state2 = check_random_state(random_state2)

    if hasattr(init1, '__array__'):
        init1 = np.asarray(init1).copy()
        if not n_init == 1:
            warnings.warn(
                'Explicit initial positions passed: '
                'performing only one init of the MDS instead of {}'.format(n_init))
            n_init = 1

    if hasattr(init2, '__array__'):
        init2 = np.asarray(init2).copy()
        if not n_init == 1:
            warnings.warn(
                'Explicit initial positions passed: '
                'performing only one init of the MDS instead of {}'.format(n_init))
            n_init = 1

    best_pos1, best_pos2, best_stress = None, None, None

    if n_jobs == 1:
        for it in range(n_init):
            pos1, pos2, stress, n_iter_ = _smacof_single(
                dissimilarities1, dissimilarities2, p, metric=metric,
                n_components=n_components, init1=init1,
                init2=init2, max_iter=max_iter, 
                verbose=verbose, eps=eps, random_state1=random_state1,
                random_state2=random_state2)
            if best_stress is None or stress < best_stress:
                best_stress = stress
                best_pos1 = pos1.copy()
                best_pos2 = pos2.copy()
                best_iter = n_iter_
    else:
        seeds1 = random_state1.randint(np.iinfo(np.int32).max, size=n_init)
        seeds2 = random_state2.randint(np.iinfo(np.int32).max, size=n_init)
        results = Parallel(n_jobs=n_jobs, verbose=max(verbose - 1, 0))(
            delayed(_smacof_single)(
                dissimilarities1, dissimilarities2, p, weights1=weights1, weights2=weights2, metric=metric, 
                n_components=n_components, init1=init1, init2=init2, 
                max_iter=max_iter, verbose=verbose, eps=eps, 
                random_state1=seed1, random_state2=seed2)
            for seed1, seed2 in zip(seeds1, seeds2))
        positions1, positions2, stress, n_iters = zip(*results)
        best = np.argmin(stress)
        best_stress = stress[best]
        best_pos1 = positions1[best]
        best_pos2 = positions2[best]
        best_iter = n_iters[best]

    if return_n_iter:
        return best_pos1, best_pos2, best_stress, best_iter
    else:
        return best_pos1, best_pos2, best_stress


class Joint_MDS(BaseEstimator):
    """Multidimensional scaling

    Read more in the :ref:`User Guide <multidimensional_scaling>`.

    Parameters
    ----------
    n_components : int, optional, default: 2
        Number of dimensions in which to immerse the dissimilarities.

    metric : boolean, optional, default: True
        If ``True``, perform metric MDS; otherwise, perform nonmetric MDS.

    n_init : int, optional, default: 4
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress.

    max_iter : int, optional, default: 300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, optional, default: 0
        Level of verbosity.

    eps : float, optional, default: 1e-3
        Relative tolerance with respect to stress at which to declare
        convergence.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging. For ``n_jobs`` below -1,
        (``n_cpus + 1 + n_jobs``) are used. Thus for ``n_jobs = -2``, all CPUs
        but one are used.

    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize the centers. If an integer is given,
        it fixes the seed. Defaults to the global numpy random number
        generator.

    dissimilarity : 'euclidean' | 'precomputed', optional, default: 'euclidean'
        Dissimilarity measure to use:

        - 'euclidean':
            Pairwise Euclidean distances between points in the dataset.

        - 'precomputed':
            Pre-computed dissimilarities are passed directly to ``fit`` and
            ``fit_transform``.

    Attributes
    ----------
    embedding_ : array-like, shape (n_components, n_samples)
        Stores the position of the dataset in the embedding space.

    stress_ : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).


    References
    ----------
    "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
    Groenen P. Springer Series in Statistics (1997)

    "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
    Psychometrika, 29 (1964)

    "Multidimensional scaling by optimizing goodness of fit to a nonmetric
    hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    """
    def __init__(self, n_components=2, weights1=None, weights2=None, p=0, metric=True, n_init=4,
                 max_iter=300, verbose=0, eps=1e-3, n_jobs=1,
                 random_state1=None, random_state2=None, 
                 dissimilarity="euclidean"):
        self.n_components = n_components
        self.weights1 = weights1
        self.weights2 = weights2
        self.p = p
        self.dissimilarity = dissimilarity
        self.metric = metric
        self.n_init = n_init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state1 = random_state1
        self.random_state2 = random_state2

    @property
    def _pairwise(self):
        return self.kernel == "precomputed"

    def fit(self, X1, X2, weights1=None, weights2=None, init=None):
        """
        Computes the position of the points in the embedding space

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        init : ndarray, shape (n_samples,), optional, default: None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.
        """
        self.fit_transform(X1, X2, weights1=weights1, weights2=weights2, init=init)
        return self

    def fit_transform(self, X1, X2, weights1=None, weights2=None, init1=None, init2=None):
        """
        Fit the data from X, and returns the embedded coordinates

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        init : ndarray, shape (n_samples,), optional, default: None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.
        """
        X1 = check_array(X1)
        if X1.shape[0] == X1.shape[1] and self.dissimilarity != "precomputed":
            warnings.warn("The MDS API has changed. ``fit`` now constructs a"
                          " dissimilarity matrix from data. To use a custom "
                          "dissimilarity matrix, set "
                          "``dissimilarity='precomputed'``.")

        if self.dissimilarity == "precomputed":
            self.dissimilarity_matrix1_ = X1
        elif self.dissimilarity == "euclidean":
            self.dissimilarity_matrix1_ = euclidean_distances(X1)
        else:
            raise ValueError("Proximity must be 'precomputed' or 'euclidean'."
                             " Got %s instead" % str(self.dissimilarity))

        X2 = check_array(X2)
        if X2.shape[0] == X2.shape[1] and self.dissimilarity != "precomputed":
            warnings.warn("The MDS API has changed. ``fit`` now constructs a"
                          " dissimilarity matrix from data. To use a custom "
                          "dissimilarity matrix, set "
                          "``dissimilarity='precomputed'``.")

        if self.dissimilarity == "precomputed":
            self.dissimilarity_matrix2_ = X2
        elif self.dissimilarity == "euclidean":
            self.dissimilarity_matrix2_ = euclidean_distances(X2)
        else:
            raise ValueError("Proximity must be 'precomputed' or 'euclidean'."
                             " Got %s instead" % str(self.dissimilarity))

        self.embedding1_, self.embedding2_, self.stress_, self.n_iter_ = smacof(
            self.dissimilarity_matrix1_, self.dissimilarity_matrix2_, p=self.p, weights1=self.weights1, 
            weights2=self.weights2, metric=self.metric, n_components=self.n_components, init1=init1, init2=init2, 
            n_init=self.n_init, n_jobs=self.n_jobs, max_iter=self.max_iter, verbose=self.verbose,
            eps=self.eps, random_state1=self.random_state1, random_state2=self.random_state2,
            return_n_iter=True)

        return self.embedding1_, self.embedding2_
