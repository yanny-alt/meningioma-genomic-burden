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
