import matplotlib.pyplot as plt
from sklearn.utils import check_array
from sklearn.utils.validation import FLOAT_DTYPES
from sica._whitening import whitening
from tqdm.notebook import tqdm
import numpy as np
from sica.base import StabilizedICA

def estimate_components(
        df, 
        min_components: int = 10, 
        max_components: int = 60, 
        step: int = 2, 
        min_mean_stability: float = 0.90, 
        n_runs: int = 20,
        algorithm: str = "fastica_par", 
        fun: str = "logcosh", 
        max_iter: int = 2000, 
        n_jobs: int = -1,
        plot_path: str = "mstd.py"
    ) -> int:
    """
    Estimate the number of components for factor analysis. Based on MSTD, see:
    - https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4112-9#Fig1
    - https://github.com/ncaptier/stabilized-ica/blob/master/sica/base.py
    - https://github.com/ncaptier/stabilized-ica/blob/master/examples/transcriptomic_ICA.ipynb

    This code selects the maximal number of components above a mean stability threshold (`min_mean_stability`).
    """

    # reticulate turns these into floats
    min_components, max_components, step, n_runs, max_iter, n_jobs = int(min_components), int(max_components), int(step), int(n_runs), int(max_iter), int(n_jobs)

    # code from MSTD function
    fig, ax = plt.subplots(1, 2, figsize=(20, 7))    
    X_values = check_array(df.values, dtype=FLOAT_DTYPES, accept_sparse=True).T
    mean = []

    # unit variance + PCA
    X_w, _ = whitening(
        X_values,
        n_components=max_components,
        svd_solver="auto",
        chunked=False,
        chunk_size=None,
        zero_center=False,
    )

    # run sICA
    for i in range(min_components, max_components + step, step):

        print(f"Testing {i} components")

        s = StabilizedICA(
            n_components=i,
            n_runs=n_runs,
            algorithm=algorithm,
            fun=fun,
            whiten=False,
            max_iter=max_iter,
            n_jobs=n_jobs
        )
        s.fit(X_w[:, :i].T)
        mean.append(np.mean(s.stability_indexes_))
        ax[0].plot(range(1, len(s.stability_indexes_) + 1), s.stability_indexes_, "k")

    ax[1].plot(range(min_components, max_components + step, step), mean)

    ax[1].set_title("Mean stability")
    ax[1].set_xlabel("Number of components")
    ax[1].axhline(y=min_mean_stability, color='r', linestyle='--')
    ax[1].xaxis.set_major_locator(plt.MultipleLocator(1))
    ax[0].set_title("Index stability distribution")
    ax[0].set_xlabel("Number of components")
    ax[0].axhline(y=min_mean_stability, color='r', linestyle='--')
    ax[0].xaxis.set_major_locator(plt.MultipleLocator(1))

    # get greatest number of components with mean stability > min_mean_stability
    n_components = max([n_components for n_components, m in zip(range(min_components, max_components + step, step), mean) if m > min_mean_stability])
    print(f"{n_components} components estimated above {min_mean_stability} mean stability")
    assert n_components < df.shape[0] * 0.9

    ax[1].axvline(x=n_components, color='b', linestyle='--')
    ax[0].axvline(x=n_components, color='b', linestyle='--')
    plt.savefig(plot_path)

    return n_components
