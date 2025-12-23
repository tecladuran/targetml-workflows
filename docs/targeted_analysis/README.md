# Targeted Analysis

This directory contains the targeted analysis workflows applied to the Tarragona samples.  
The objective of this analysis is to extract biologically meaningful conclusions from the concentrations of selected analytes that are investigated as potential biomarkers for colorectal cancer (CRC).

The workflow integrates GC-IMS quantitative data of Anisole and 2-Heptanone, quantified analytes from tarragona measurements, and patient metadata.  
Data loading, imputation, preprocessing, exploratory analysis, statistical testing, and classification follow the same pipeline in all workflows.  
The only difference between the workflows is the normalization strategy applied.

Statistical differences between CRC and control groups are evaluated using non-parametric tests, and classification performance is assessed using nested cross-validation to ensure robust and unbiased model evaluation.

---

## Statistical Analysis

Group differences between CRC and control samples are assessed using the **Wilcoxon rank-sum test** applied independently to each compound.  
P-values are adjusted for multiple testing using **false discovery rate (FDR)** correction.

---

## Classification Performance

Classification performance is evaluated using **Random Forest models** with a **nested cross-validation** strategy.  
The dataset consists of 30 samples and is evaluated using an outer cross-validation scheme with **7 folds**, followed by **6 folds** in a second stage.

At each outer iteration, **four samples are held out for testing**, stratified by **condition (CRC vs control)** and **sex** to preserve group balance.  
These stratification variables are configurable and can be easily modified.

---

## Workflows

### [full_workflow_tgn_targeted_pqn.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/targeted_analysis/full_workflow_tgn_targeted_pqn.md)

Complete targeted analysis workflow using **Probabilistic Quotient Normalization (PQN)**.  


---

### [full_workflow_tgn_targeted_creatinine.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/targeted_analysis/full_workflow_tgn_targeted_creatinine.md)

Complete targeted analysis workflow using **creatinine normalization**.  
