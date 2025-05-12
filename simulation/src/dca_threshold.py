import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

def find_breakpoint(xs, ys, min_frac=0.9, max_frac=1):
    """
    Brute‐force search over candidate breakpoints in [min_frac, max_frac] fraction
    of the data to find the piecewise‐linear split that minimizes total SSE.
    Returns best index i such that [0:i] is the 'noise' region.
    """
    N = len(xs)
    i_min = int(min_frac * N)
    i_max = int(max_frac * N)
    print('N:' + str(N))
    print('i_min:' + str(i_min))
    print('i_max:' + str(i_max))
    best_sse = np.inf
    best_i = None

    for i in range(i_min, i_max):
        # fit first segment
        xr1 = xs[:i].reshape(-1,1)
        yr1 = ys[:i]
        m1 = LinearRegression().fit(xr1, yr1)
        pred1 = m1.predict(xr1)
        sse1 = np.sum((yr1 - pred1)**2)

        # fit second segment
        xr2 = xs[i:].reshape(-1,1)
        yr2 = ys[i:]
        m2 = LinearRegression().fit(xr2, yr2)
        pred2 = m2.predict(xr2)
        sse2 = np.sum((yr2 - pred2)**2)

        if i_max - i < 10:
            print('i:' + str(i))
            print('sse1:' + str(sse1))
            print('sse2:' + str(sse2))
            print('sse1 + sse2:' + str(sse1 + sse2))
            print('pred2:' + str(pred2))
            print('yr2:' + str(yr2))

        if sse1 + sse2 < best_sse:
            best_sse = sse1 + sse2
            best_i = i

    return best_i

def six_sigma_cutoff(couplings, plot=False):
    # 1. take absolute values and sort
    mags = np.abs(couplings)
    xs = np.sort(mags)
    N = len(xs)
    # 2. compute semi‑log cumulative counts
    counts = np.arange(N, 0, -1)
    log_counts = np.log(counts)

    # 3. find best breakpoint index
    bp_idx = find_breakpoint(xs, log_counts)
    bp_value = xs[bp_idx]

    # 4. estimate noise σ from the 'noise' region [0:bp_idx]
    noise_vals = mags[mags <= bp_value]
    sigma = np.std(noise_vals, ddof=1)

    # 5. threshold = 6σ
    threshold = 6 * sigma
    print('threshold:' + str(threshold))
    print('breakpoint:' + str(bp_value))
    mask = mags >= threshold

    if plot:
        # plot semi‑log CDF and fitted spline
        plt.figure(figsize=(6,4))
        plt.plot(xs, counts, label='empirical CDF')
        plt.yscale('log')
        # overlay the two linear fits
        from sklearn.linear_model import LinearRegression
        # first
        m1 = LinearRegression().fit(xs[:bp_idx].reshape(-1,1), log_counts[:bp_idx])
        plt.plot(xs[:bp_idx], np.exp(m1.predict(xs[:bp_idx].reshape(-1,1))),
                 'r--', label='noise fit')
        # second
        m2 = LinearRegression().fit(xs[bp_idx:].reshape(-1,1), log_counts[bp_idx:])
        plt.plot(xs[bp_idx:], np.exp(m2.predict(xs[bp_idx:].reshape(-1,1))),
                 'g--', label='signal fit')
        plt.axvline(threshold, color='k', linestyle=':', label='6 sigma =' + str(threshold))
        plt.xlabel('Coupling magnitude')
        plt.ylabel('Count ≥ x')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return threshold, mask

def elbow_threshold(couplings):
    # 1. take absolute values and sort
    mags = np.sort(np.abs(couplings))
    # 2. compute the diffs
    diffs = np.diff(mags)
    # 3. find the largest diff
    elbow_idx = np.argmax(diffs)
    # set the threshold as the value at the elbow index
    threshold = mags[elbow_idx]
    return threshold
    
def generate_test_for_threshold(n_noise=1000,
                                noise_scale=1.0,
                                n_signal=50,
                                signal_low=6.0,
                                signal_high=10.0,
                                random_state=42):
    """
    Generate a synthetic dataset of coupling scores with known noise vs signal.
    
    Returns a pandas DataFrame with:
      - 'coupling_score': float value
      - 'is_signal'     : 1 if drawn from the high-magnitude 'signal' distribution
                         0 if drawn from the low-magnitude 'noise' distribution
      - 'label'         : 'signal' or 'noise'
    """
    np.random.seed(random_state)

    # 1) Noise: symmetric exponential around zero
    noise = np.random.exponential(scale=noise_scale, size=n_noise)
    # randomly flip sign
    noise *= np.random.choice([-1, 1], size=n_noise)
    noise_labels = np.zeros(n_noise, dtype=int)

    # 2) Signal: uniform well above the expected 6σ threshold
    signal = np.random.uniform(low=signal_low, high=signal_high, size=n_signal)
    signal *= np.random.choice([-1, 1], size=n_signal)
    signal_labels = np.ones(n_signal, dtype=int)

    # 3) Combine and shuffle
    scores = np.concatenate([noise, signal])
    labels = np.concatenate([noise_labels, signal_labels])
    order = np.random.permutation(len(scores))
    scores = scores[order]
    labels = labels[order]

    # 4) Build DataFrame
    df = pd.DataFrame({
        'coupling_score': scores,
        'is_signal': labels
    })
    df['label'] = df['is_signal'].map({0: 'noise', 1: 'signal'})
    return df


# --- Usage example ---
# suppose `J` is your flat array of all coupling scores
# threshold, keep_mask = six_sigma_cutoff(J)
# significant_J = J[keep_mask]
