# Mathematical methodology and evaluation

## Representation
Let X in R^(n×T), n=4381, T=23. Rows are transcript profiles; columns are ordered log-ratio measurements, not independent biological replicates. No inference about molecular geometry is possible. Raw PCA uses column centering; anomaly PCA uses feature-wise standardization fitted on training genes only.

## Feature engineering
Ten features: mean, population standard deviation, range, minimum, maximum, standard deviation of first differences, maximum absolute difference, mean absolute difference, normalized lag-1 autocovariance, dominant non-DC Fourier power fraction. The last feature is descriptive, not proof of periodicity. Feature scaling uses training mean and variance. No biological labels are inferred.

## Co-expression graph
Pearson r_ij = sum_t[(x_it-mean_i)(x_jt-mean_j)] / sqrt(sum_t(x_it-mean_i)^2 sum_t(x_jt-mean_j)^2). Constant rows yield zero associations. Retain |r|>=0.8 among each node's top ten absolute associations, then union directed selections into an undirected graph. Degrees can exceed ten. Weight w_ij=|r_ij|; original signed r remains in exports. Positive and negative correlations are treated as strong associations for partitioning, not as equivalent biological regulation. Louvain maximizes weighted modularity Q=(1/2m) sum_ij [w_ij - k_i k_j/(2m)] delta(c_i,c_j), with m=sum edge weights. Isolates receive label -1. PageRank uses damping 0.85 (NetworkX default); centrality is not biological essentiality. Density=2E/[n(n-1)]. No p-values are claimed: temporal autocorrelation and multiple testing make naive Pearson significance misleading.

## 3D embeddings
PCA solves covariance eigendecomposition; coordinates Z=(X-column_mean)V_3. Explained variance is the sum of the top three eigenvalues divided by the total. Trustworthiness measures unexpected neighbors in the low-dimensional representation (k=10); it is not biological accuracy. Optional UMAP is exploratory and seeded. Coordinates remain fixed as expression changes over time.

## Unsupervised anomaly models
Isolation Forest isolates unusual profiles with short random-tree path lengths. LOF compares local reachability density with neighbors (k=25, novelty mode). PCA anomaly error = mean_j[(z_ij - reconstructed_z_ij)^2] using 8 components and standardized profiles. Models are fitted on 60% of genes (2,628), calibrated on a disjoint 15% (657), and tested on 25% (1,096). Each raw score is mapped to its empirical percentile against calibration scores; the ensemble is the arithmetic mean of three percentiles. Threshold 95 is a heuristic; ensemble percentiles do not necessarily flag exactly 5%. LOF reference scores are computed on unseen calibration genes. Scores are not probabilities.

## Evaluation
Held-out original profiles serve as proxy negatives; paired corruptions are positives. Three equally assigned perturbations: a single 6-SD spike, linear drift up to 5 SD, shuffled time order, with SD floor .15. These are synthetic stress tests, not biological labels. ROC-AUC measures ranking; average precision summarizes the precision-recall curve; F1=2PR/(P+R) at fixed rank threshold 95. Original profiles may already be atypical. Corruption strength is substantial, and paired samples are dependent. Report per-type and severity metrics and grouped uncertainty in subsequent work. Scores are empirical ranks on a separate calibration split, not calibrated probabilities; external deployment validity is unknown.

## Robustness
Remove contiguous 3-column blocks at four positions, rebuild graph and compare edge Jaccard and adjusted Rand index for shared non-isolates. This probes sensitivity to temporal sampling, not confidence intervals or independent replication. The measured instability is material. Validate against independent experiments, signed community methods, and curated annotations before claiming biological findings. Never equate expression correlation with causal regulation.

## References
Spellman et al. (1998), Comprehensive Identification of Cell Cycle-regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray Hybridization. https://pmc.ncbi.nlm.nih.gov/articles/PMC25624/

Liu et al. (2008), Isolation Forest. Breunig et al. (2000), LOF. Blondel et al. (2008), Fast unfolding of communities in large networks. Jolliffe & Cadima (2016), Principal component analysis: a review and recent developments.

## Explicit detector equations

Isolation Forest: s(x,n)=2^(-E[h(x)]/c(n)), where h is the random-tree path length and c(n) is the expected path length in a random binary search tree. The implementation negates scikit-learn score_samples so larger values are more unusual. 150 trees; default subsampling.

LOF: reach-dist_k(x,o)=max(k-distance(o),d(x,o)); LRD_k(x)=1/mean_o reach-dist_k(x,o); LOF_k(x)=mean_o LRD_k(o)/LRD_k(x). Euclidean distance on ten standardized features, k=25. Novelty mode scores unseen profiles.

PCA: with standardized input z, W contains eight orthonormal principal directions; reconstruction z_hat=mu+WW^T(z-mu). Score is ||z-z_hat||²/23. The separate three-component visualization PCA is fit on the entire dataset for descriptive exploration, not used by the held-out anomaly benchmark.

Calibration: F_j(s)=100/N_cal * sum_i 1[s_ij<=s]. Ensemble = mean_j F_j(s_j). The threshold of 95 is fixed and not optimized. Average ranks need not follow a uniform distribution. Precision=TP/(TP+FP), recall=TP/(TP+FN), F1=2PR/(P+R). AP is the recall-weighted sum of precision values, not trapezoidal PR-AUC. ROC-AUC is the probability of ranking a random synthetic positive above a proxy negative, with half credit for ties.

PageRank: p=(1-d)/n + d P^T p, where P is the weighted transition matrix with dangling nodes redistributed uniformly and d=0.85. Fourier power is |DFT(x-mean(x))|². The graph is computed in O(n²T) time and O(n²) memory, with row sorting O(n² log n). The three detectors are independent of graph features, avoiding use of the full-data graph during held-out detector evaluation.

## Persistence and deployment

The saved bundle contains the exact evaluated detector estimators, both scalers, the calibration reference scores, timepoint schema, and descriptive PCA embedding. It is not refit on the test set. SHA-256 checks detect file changes but do not establish trust; never load unknown Joblib files. Inference expects compatible yeast log-ratio profiles on identical timepoints. Optional UMAP visualization does not alter the saved PCA inference coordinates. No biological causal, diagnostic or clinical claims are supported.
