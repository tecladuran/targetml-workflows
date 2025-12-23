# Untargeted Analysis

This directory contains untargeted analysis workflows applied to the Tarragona samples using GC-IMS data.  
In contrast to targeted analyses, no prior assumptions on biomarkers are made and all GC-IMS features are jointly analyzed.

Two untargeted strategies are implemented depending on whether clinical metadata is incorporated into the analysis.  
All workflows share a common structure in terms of preprocessing, exploratory analysis, feature selection, and performance evaluation.

---

## Analysis Strategies

### [gcims](https://github.com/tecladuran/targetml-workflows/tree/main/docs/untargeted_analysis/gcims)

Untargeted analysis based exclusively on **GC-IMS data**.  
Classification is performed using **PLS-DA models**, without incorporating metadata as predictive variables.

---

### [gcims_metadata](https://github.com/tecladuran/targetml-workflows/tree/main/docs/untargeted_analysis/gcims_metadata)

Untargeted analysis combining **GC-IMS data and clinical metadata**.  
In this case, classification is performed using **Random Forest models**, and metadata variables can be used for stratification during cross-validation.

---

## Feature Selection

All untargeted workflows apply **recursive feature elimination (RFE)** to progressively reduce the dimensionality of the GC-IMS feature space.  
Feature importance is assessed within a nested cross-validation framework to avoid information leakage and ensure robust feature selection.

---

## Model Evaluation

Model performance is evaluated using **nested cross-validation**, with separate outer and inner loops for performance estimation and model tuning.  
This strategy ensures unbiased evaluation of classification performance while integrating feature selection within the training process.

---

## Notes

Normalization strategies (PQN or creatinine-based) are the only differences between workflows within the same analysis strategy, enabling direct comparison of their impact on downstream results.
