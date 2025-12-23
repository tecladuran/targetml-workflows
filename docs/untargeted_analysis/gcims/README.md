# Untargeted Analysis — GC-IMS Only

This directory contains untargeted analysis workflows based exclusively on **GC-IMS data**, without the inclusion of clinical metadata.

The analysis focuses on evaluating the discriminative power of the full GC-IMS peak table using multivariate modeling and machine learning techniques.

Classification is performed using **PLS-DA models**, and feature selection is carried out through **recursive feature elimination (RFE)** within a nested cross-validation framework.

---

## Workflows

### [tgn_untargeted_plsda_pqn.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/untargeted_analysis/gcims/tgn_untargeted_plsda_pqn.md)

Complete untargeted GC-IMS workflow using **Probabilistic Quotient Normalization (PQN)**.

---

### [tgn_untargeted_plsda_creatinine.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/untargeted_analysis/gcims/tgn_untargeted_plsda_creatinine.md)

Complete untargeted GC-IMS workflow using **creatinine normalization**.
