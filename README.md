# Meningioma Genomic Burden Analysis


**Multi-omics analysis demonstrating TMB and CNA as WHO grade biomarkers in meningiomas**

---

## 🎯 Project Summary

Comprehensive bioinformatics analysis of 123 meningioma samples identifying **tumor mutational burden (TMB)** and **copy number alterations (CNA)** as quantitative biomarkers for WHO grade stratification.

### Key Findings at a Glance

| Metric | Grade I (n=56) | Grade II/III (n=62) | Fold-Change | p-value | Significance |
|--------|----------------|---------------------|-------------|---------|--------------|
| **TMB** (mut/Mb) | 0.222 | 0.450 | **2.03×** | **1.28×10⁻⁴** | *** |
| **CNA burden** (FGA) | 0.098 | 0.166 | **1.70×** | **1.10×10⁻⁴** | *** |
| **ML Classification** | - | - | AUC=**0.879** | Accuracy=**76.5%** | ✓ |

*FGA = Fraction Genome Altered*

---

## 📄 Abstract

**Background:** Meningiomas represent the most common primary intracranial tumor, yet molecular predictors of clinical aggressiveness remain poorly characterized. We conducted a comprehensive multi-omics analysis to investigate whether tumor mutational burden (TMB) and copy number alterations (CNA) correlate with WHO grade classification.

**Methods:** We obtained multi-omics data from 121 meningioma patients (95.0% with somatic mutations, n=115; 100% with copy number data, n=121; 79.3% with protein abundance data, n=96) from the University of Toronto cBioPortal study. Supplementary data included 2 TCGA meningioma samples with whole exome sequencing and 51 meningiomas from Gene Expression Omnibus datasets (GSE43290, n=47; GSE88720, n=14). Clinical annotations included WHO grade (Grade I: n=59; Grade II: n=43; Grade III: n=19; missing: n=2) and recurrence status. Statistical comparisons used the Wilcoxon rank-sum test with significance threshold p<0.05. Machine learning classification was performed using logistic regression, random forest, and XGBoost models on 114 samples with complete genomic data.

**Results:** High-grade meningiomas (WHO Grade II/III, n=62) demonstrated significantly elevated TMB compared to Grade I tumors (n=56): mean 0.450 vs. 0.222 mutations/Mb (fold-change 2.03×, p=1.28×10⁻⁴). Similarly, CNA burden was substantially higher in high-grade tumors: mean fraction genome altered 0.166 vs. 0.098 (fold-change 1.70×, p=1.10×10⁻⁴). Logistic regression achieved the highest classification performance with AUC=0.879 and 76.5% accuracy on the test set (n=34), significantly exceeding the discrimination threshold of AUC>0.70.

**Conclusion:** Genomic variant burden effectively stratifies meningioma aggressiveness and demonstrates clinical utility for grade prediction. These findings support TMB and CNA as quantitative biomarkers for identifying high-risk meningiomas and may guide surveillance protocols and treatment decisions.


## 📊 Results Visualization

### Figure 2: TMB by WHO Grade
[View TMB by Grade Figure](https://github.com/yanny-alt/meningioma-genomic-burden/blob/main/results/figures/figure2_tmb_by_grade.pdf)

**High-grade meningiomas show 2.03× higher TMB (p<0.001)**

### Figure 3: CNA Burden by WHO Grade  
[View CNA Burden Figure](https://github.com/yanny-alt/meningioma-genomic-burden/blob/main/results/figures/figure3_cna_burden_by_grade.pdf)

**High-grade tumors have 1.70× higher CNA burden (p<0.001)**

### Figure 4: Integrated Genomic Landscape
[View Genomic Landscape Figure](https://github.com/yanny-alt/meningioma-genomic-burden/blob/main/results/figures/figure4_integrated_genomic_burden.pdf)

**TMB and CNA correlate with grade severity**

### Figure 5: ML Classification Performance
[View ROC Curves Figure](https://github.com/yanny-alt/meningioma-genomic-burden/blob/main/results/figures/figure5_ml_roc_curves.pdf)

**Logistic regression achieves AUC=0.879 for grade prediction**

---

## 🔬 Methodology

### Data Sources

| Dataset | Samples | Data Types | Source |
|---------|---------|------------|--------|
| **cBioPortal** (UofT) | 121 | Mutations, CNA, Protein | [Study Link](https://www.cbioportal.org/study/summary?id=mng_utoronto_2021) |
| **TCGA** | 2 | WES, CNV segments | [GDC Portal](https://portal.gdc.cancer.gov/) |
| **GSE43290** | 47 | Microarray expression | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43290) |
| **GSE88720** | 14 | Microarray expression | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88720) |

**Total: 123 unique meningioma samples**

### Analysis Pipeline

Data Acquisition → Quality Control → TMB Calculation → CNA Quantification
        ↓                                                      ↓
  Statistical Testing ← Data Integration ← Feature Engineering
        ↓
  ML Classification (Logistic Regression, Random Forest, XGBoost)
        ↓
  Visualization & Interpretation


## 🔬 Key Metrics

**TMB Calculation:**
- Total somatic mutations / coding genome size (38 Mb)
- Reported as mutations per megabase (mut/Mb)

**CNA Burden:**
- Fraction of genome altered (FGA)
- Sum of amplified + deleted segments / total genome

**Statistical Tests:**
- Wilcoxon rank-sum test (non-parametric)
- Significance threshold: p < 0.05

**ML Models:**
- 70-30 train-test split
- 5-fold cross-validation
- Performance metrics: AUC, accuracy, sensitivity, specificity

---

## 📈 Detailed Results

### Primary Hypothesis (H1): TMB Difference by Grade

**Test:** Wilcoxon rank-sum test  
**Result:** p = 1.28×10⁻⁴ (highly significant)

| Grade | n | Mean TMB | Median TMB | SD |
|-------|---|----------|------------|-----|
| Grade I | 56 | 0.222 | 0.184 | 0.156 |
| Grade II/III | 62 | 0.450 | 0.398 | 0.312 |

**Interpretation:** High-grade meningiomas have **2.03× higher** TMB than low-grade tumors.

### Secondary Hypothesis (H2): CNA Burden Difference by Grade

**Test:** Wilcoxon rank-sum test  
**Result:** p = 1.10×10⁻⁴ (highly significant)

| Grade | n | Mean FGA | Median FGA | SD |
|-------|---|----------|------------|-----|
| Grade I | 56 | 0.098 | 0.075 | 0.089 |
| Grade II/III | 62 | 0.166 | 0.142 | 0.124 |

**Interpretation:** High-grade meningiomas have **1.70× higher** CNA burden.

### Tertiary Hypothesis (H3): ML Classification Performance

| Model | AUC | Accuracy | Sensitivity | Specificity |
|-------|-----|----------|-------------|-------------|
| **Logistic Regression** | **0.879** | **76.5%** | 78.9% | 73.3% |
| Random Forest | 0.845 | 73.5% | 73.7% | 73.3% |
| XGBoost | 0.832 | 70.6% | 68.4% | 73.3% |

**Interpretation:** Genomic burden features achieve **excellent discrimination** (AUC>0.85) for grade classification.

---

## 🏢 Internship Context

This project was completed as part of a **Cancer Genomics Research Internship** at **Genomac Institute Inc.** (July-October 2025).

**Responsibilities:**
- Multi-omics data integration (TCGA, cBioPortal, GEO)
- Statistical analysis and hypothesis testing
- Machine learning model development
- Manuscript preparation and figure generation

**Key Achievements:**
- ✅ Identified TMB and CNA as significant biomarkers (p<0.001)
- ✅ Developed ML classification framework (AUC=0.879)
- ✅ Generated publication-ready manuscript
- ✅ Created reproducible analysis pipeline

---

## 📊 Supplementary Materials

### Tables
- **Table S1:** Sample characteristics by dataset and grade
- **Table S2:** Complete TMB and CNA values for all samples
- **Table S3:** Statistical test results (Wilcoxon tests)
- **Table S4:** ML model comparison (3 algorithms)

All supplementary tables available in [`results/tables/`](results/tables/)

### Figures
- **Figure 1:** Study design and cohort overview (planned)
- **Figure 2:** TMB distribution by WHO grade
- **Figure 3:** CNA burden by WHO grade
- **Figure 4:** Integrated genomic burden landscape
- **Figure 5:** ML classification ROC curves

All figures available in [`results/figures/`](results/figures/)

---

## 💡 Clinical Implications

### Diagnostic Utility
- **TMB and CNA** provide quantitative measures of tumor aggressiveness
- **ML classification** offers objective grade prediction (76.5% accuracy)

### Potential Applications
1. **Risk Stratification:** Identify high-risk Grade I tumors needing closer monitoring
2. **Treatment Planning:** Guide adjuvant therapy decisions for borderline cases
3. **Surveillance Protocols:** Personalize follow-up intervals based on genomic burden
4. **Clinical Trials:** Stratify patients for targeted therapy studies

### Future Directions
- Validate in independent cohorts
- Integrate with clinical variables (age, sex, location)
- Explore specific driver mutations
- Investigate RNA-seq signatures
  
