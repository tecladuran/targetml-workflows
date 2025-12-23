# Untargeted Analysis — GC-IMS + Metadata

This directory contains untargeted analysis workflows that combine **GC-IMS data and clinical metadata**.

Classification is performed using **Random Forest models**, while GC-IMS features are selected through **recursive feature elimination (RFE)**.  
Clinical metadata variables can be used to define stratification schemes during cross-validation, allowing flexible control of confounding factors such as sex.

The analytical pipeline is otherwise equivalent to the GC-IMS–only workflows, enabling direct comparison of model performance with and without metadata inclusion.

---

## Workflows

### [tgn_untargeted_rf_pqn.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/untargeted_analysis/gcims_metadata/tgn_untargeted_rf_pqn.md)

Untargeted GC-IMS + metadata workflow using **Probabilistic Quotient Normalization (PQN)**.

---

### [tgn_untargeted_rf_creatinine.md](https://github.com/tecladuran/targetml-workflows/blob/main/docs/untargeted_analysis/gcims_metadata/tgn_untargeted_rf_creatinine.md)

Untargeted GC-IMS + metadata workflow using **creatinine normalization**.
