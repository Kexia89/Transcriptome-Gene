# Transcriptome Atlas

**Temporal Transcriptomics · Mathematical Modeling · Network Science · Cancer-Oriented Computational Research**

<p align="center">
  <img src="artifacts/figures/transcriptome_atlas_banner.png" alt="Transcriptome Atlas — conceptual illustration of gene networks and future cancer research" width="100%">
</p>

*Conceptual artwork; not an experimental visualization or clinical result.*

## Project overview

**Transcriptome Atlas** is an independent computational biology research project focused on the structure, dynamics, and interpretability of temporal gene-expression data. The current study uses the Spellman yeast transcriptomics dataset as a mathematically tractable research substrate for exploring how high-dimensional expression profiles organize into co-expression communities, low-dimensional manifolds, and statistically unusual patterns.

The project combines **graph theory**, **multivariate statistics**, **dimensionality reduction**, and **unsupervised anomaly detection** into a single analytical framework. Rather than presenting a software tutorial, this repository documents a reproducible scientific exploration and its quantitative findings.

Although the present experiments were conducted on **_Saccharomyces cerevisiae_** and **not on cancer cohorts**, the analytical logic offers a potential methodological starting point for oncology. Temporal gene-expression modeling, network-level organization, and anomaly scoring are core ideas that can later be extended to tumor evolution, subtype stratification, treatment response, and systems-level biomarker discovery.

---

## Research focus

The central question of the project is whether temporal expression profiles can reveal:

- coherent **co-expression communities**,
- mathematically meaningful **three-dimensional transcriptome structure**, and
- statistically unusual **expression behaviors** that may warrant biological investigation.

The present version therefore serves as a **research prototype** for transcriptome intelligence: a framework that begins with a controlled temporal dataset and opens a path toward future translational applications in cancer systems biology.

---

## Main results

The current dataset contains **4,381 genes** measured across **23 time points**. Under the present graph construction and anomaly-detection protocol, the pipeline produced the following baseline results:

- **17,977 co-expression edges**
- **78 non-isolated communities**
- **Weighted modularity: 0.6375**
- **3D PCA retained variance: 58.41%**
- **Trustworthiness (k = 10): 0.8893**
- **93 genes flagged** at the ensemble rank threshold of 95

### Anomaly-detection performance

The anomaly benchmark is based on **synthetic perturbations** introduced into held-out profiles. These are not biological anomaly labels and should not be interpreted as disease-detection accuracy.

| Detector | ROC-AUC | Average Precision | F1 at 95 | Precision | Recall |
|---|---:|---:|---:|---:|---:|
| Isolation Forest | 0.7840 | 0.8084 | 0.5254 | 0.9140 | 0.3686 |
| Local Outlier Factor | 0.8363 | 0.8760 | 0.7403 | 0.9236 | 0.6177 |
| **PCA Reconstruction Error** | **0.9469** | **0.9468** | **0.8023** | 0.9367 | **0.7016** |
| Ensemble | 0.8817 | 0.8964 | 0.6851 | **0.9620** | 0.5319 |

The strongest aggregate performer in this experimental setting was **PCA reconstruction error**, suggesting that deviations from the dominant low-dimensional structure are particularly informative for detecting unusual temporal profiles.

### Experimental design

The detectors were fitted on 2,628 gene profiles, calibrated on 657 separate profiles and evaluated on 1,096 held-out genes paired with synthetic perturbations (2,192 evaluation profiles). The random seed was 42. PCA reconstruction achieved the strongest overall ROC-AUC, while the ensemble provided the highest precision at the fixed score threshold. These results apply to this balanced synthetic benchmark only.

The full numerical results are available in [metrics.json](artifacts/metrics.json). The experimental protocol is documented in [Mathematical methodology](METHODOLOGY.md) and the [Model card](MODEL_CARD.md).

---

## Mathematical framework

The repository is built around a mathematical interpretation of transcriptomic structure.

### 1. Co-expression graph
Gene–gene similarity is estimated through Pearson correlation across time. A weighted graph is then constructed by retaining high-magnitude associations under a thresholded top-k strategy. Community structure is estimated with **weighted Louvain modularity optimization**, and node relevance is summarized through graph statistics such as degree and PageRank.

### 2. Transcriptome geometry
A three-dimensional representation of the expression space is generated using **principal component analysis (PCA)**, producing a geometric view of neighborhood structure. This low-dimensional space is not a physical map of molecules or cells, but a mathematical embedding of profile similarity.

### 3. Anomaly detection
Three complementary unsupervised detectors are applied:

- **Isolation Forest** for random-partition anomaly isolation,
- **Local Outlier Factor** for density-based irregularity,
- **PCA reconstruction error** for structural deviation from the dominant linear manifold.

Their outputs are calibrated with empirical percentile scoring and summarized in an ensemble rank.

### 4. Time-series descriptors
The models are supplied with descriptive features derived from each gene profile, including:

- mean and standard deviation,
- amplitude and extrema,
- first-difference variability,
- lag-1 temporal dependence,
- dominant non-DC Fourier power fraction.

This design allows the project to preserve a balance between **interpretability** and **algorithmic expressiveness**.

---

## Scientific interpretation

The present findings indicate that the dataset contains a substantial level of transcriptomic organization. The modularity score describes the community structure of the thresholded graph; without a null-model comparison, it does not establish statistical significance, while the trustworthiness score indicates that the three-dimensional representation preserves a meaningful portion of local neighborhood structure.

At the same time, the block-deletion stability analysis shows that the graph is **materially sensitive to temporal sampling**. This is an important scientific caution: co-expression structure can shift depending on which parts of the temporal trajectory are observed. Consequently, the graph should be interpreted as **exploratory structure**, not as a validated regulatory network.

The anomaly-detection results are encouraging as a modeling exercise, especially for PCA-based reconstruction scoring, but they remain constrained by the use of synthetic perturbations. The current benchmark demonstrates mathematical discrimination, not clinical or biological diagnosis.

---

## Why this matters for cancer research

Even though the current project does not use tumor data, its methodological architecture aligns naturally with several cancer-oriented research directions.

### A. Tumor-state representation
A transcriptome atlas can be used to represent tumor samples in a structured space where proximity reflects molecular similarity. This is relevant for identifying **cancer subtypes**, **cellular states**, and **disease trajectories**.

### B. Network-level oncology
Cancer is increasingly understood as a disease of dysregulated systems rather than isolated genes. Co-expression communities and graph-based representations can help reveal **pathway-level disruption**, **rewiring of regulatory modules**, and **emergent tumor programs**.

### C. Outlier biology
Anomaly scoring can be used to identify transcriptional behaviors that deviate from the dominant background, which may be relevant to:

- rare tumor programs,
- resistant cell populations,
- unusual treatment responses,
- aggressive or transitional phenotypes.

### D. Mathematical oncology
The broader vision of the project is compatible with future integration of:

- longitudinal tumor transcriptomics,
- single-cell RNA-seq,
- spatial transcriptomics,
- dynamical systems models of tumor evolution,
- therapy-response prediction,
- biomarker discovery pipelines.

In this sense, **Transcriptome Atlas** is not merely a data-visualization project; it is a mathematical and computational foundation for future work in **cancer systems biology**.

---

## Future directions

The next phase of this project can move from a controlled yeast setting toward more directly oncological applications. Priority directions include:

1. **External validation on cancer datasets**, such as TCGA, GEO, or curated breast-cancer transcriptome cohorts.
2. **Comparative subtype mapping**, where embedding spaces are used to distinguish tumor classes and intermediate states.
3. **Network perturbation analysis**, to study how communities change between tumor and normal conditions.
4. **Cancer-specific anomaly modeling**, to identify rare transcriptional signatures linked to resistance, metastasis, or relapse.
5. **Integration with survival and clinical variables**, allowing graph and anomaly measures to be related to prognosis.
6. **Extension to single-cell and spatial omics**, where neighborhood structure and outlier states become even more informative.
7. **Hybrid mathematical models**, combining transcriptomic embeddings with differential equations, stochastic processes, or evolutionary modeling.

### Candidate mathematical models for oncology (future work)

These equations describe proposed research extensions, not models fitted or validated in the current yeast experiment.

**Differential expression.** For appropriately normalized positive expression quantities, a descriptive tumor-to-normal contrast is

$$
\Delta_g=\log_2\left(\frac{\mu_{g,T}}{\mu_{g,N}}\right).
$$

Count-based RNA-seq inference would require a suitable statistical model, covariate adjustment and multiple-testing correction.

**Tumor population dynamics.** A possible starting point is a logistic growth model with a time-dependent treatment effect:

$$
\frac{dN}{dt}=rN\left(1-\frac{N}{K}\right)-u(t)N.
$$

Here $N$ is tumor population size, $r$ is growth rate, $K$ is carrying capacity and $u(t)$ is an effective treatment-associated loss rate. Longitudinal data would be needed to estimate and validate these parameters.

**Prognostic modeling.** A Cox proportional hazards model could relate transcriptomic features to time-to-event outcomes:

$$
h(t\mid z)=h_0(t)\exp(\beta^\top z).
$$

Independent patient-level validation, censoring-aware evaluation and assessment of proportional hazards would be essential. None of these proposed models currently provides a clinical prediction.

The long-term objective is to transform transcriptomic structure into a quantitative language for understanding how cancer organizes, adapts, and responds to intervention.

---

## Repository outputs

The project generates:

- trained unsupervised models,
- graph and node-level transcriptomic summaries,
- evaluation metrics,
- publication-ready PNG/SVG figures,
- an interactive 3D dashboard,
- methodological and model documentation.

The emphasis of the repository is therefore not only on code, but on **scientific evidence, mathematical transparency, and interpretable computational results**.

---

## Limitations

Several limitations must be stated clearly.

- The source dataset is **yeast**, not human cancer.
- The 3D embedding is a **mathematical projection**, not a physical biological atlas.
- Co-expression does **not** establish causal regulation.
- The anomaly benchmark uses **synthetic perturbations**, not verified biological outliers.
- Temporal sampling sensitivity remains non-trivial and affects graph stability.
- Redistribution rights for the raw dataset and some derived artifacts should be verified before public release.

These limitations do not diminish the value of the repository; rather, they define its correct scientific scope.

---

## Bibliography

1. Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R., Anders, K., Eisen, M. B., Brown, P. O., Botstein, D., & Futcher, B. (1998). *Comprehensive Identification of Cell Cycle-regulated Genes of the Yeast Saccharomyces cerevisiae by Microarray Hybridization*. Molecular Biology of the Cell, 9(12), 3273–3297. https://pmc.ncbi.nlm.nih.gov/articles/PMC25624/

2. Hanahan, D., & Weinberg, R. A. (2011). *Hallmarks of Cancer: The Next Generation*. Cell, 144(5), 646–674. https://doi.org/10.1016/j.cell.2011.02.013

3. Hoadley, K. A., Yau, C., Wolf, D. M., et al. (2014). *Multiplatform Analysis of 12 Cancer Types Reveals Molecular Classification within and across Tissues of Origin*. Cell, 158(4), 929–944. https://doi.org/10.1016/j.cell.2014.06.049

4. Breunig, M. M., Kriegel, H.-P., Ng, R. T., & Sander, J. (2000). *LOF: Identifying Density-Based Local Outliers*. ACM SIGMOD Record, 29(2), 93–104. https://doi.org/10.1145/335191.335388

5. Liu, F. T., Ting, K. M., & Zhou, Z.-H. (2008). *Isolation Forest*. 2008 Eighth IEEE International Conference on Data Mining, 413–422. https://doi.org/10.1109/ICDM.2008.17

6. Jolliffe, I. T., & Cadima, J. (2016). *Principal Component Analysis: A Review and Recent Developments*. Philosophical Transactions of the Royal Society A, 374(2065), 20150202. https://doi.org/10.1098/rsta.2015.0202

7. Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). *Fast Unfolding of Communities in Large Networks*. Journal of Statistical Mechanics: Theory and Experiment, 2008(10), P10008. https://doi.org/10.1088/1742-5468/2008/10/P10008

8. Cox, D. R. (1972). *Regression Models and Life-Tables*. Journal of the Royal Statistical Society: Series B, 34(2), 187–220. https://doi.org/10.1111/j.2517-6161.1972.tb00899.x

9. Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology, 15, 550. https://doi.org/10.1186/s13059-014-0550-8

---

**Transcriptome Atlas** presents a computational and mathematical perspective on temporal gene-expression structure today, while pointing toward a future in which transcriptome geometry, network organization, and anomaly detection can contribute to the understanding of cancer.
