# targetml-workflows

This repository contains **workflows and resources for TargetML data analysis**.  
Code, data, and documentation are organized in a modular way to facilitate reuse, readability, and reproducibility across projects.

---

## Repository structure

### [R](https://github.com/tecladuran/targetml-workflows/tree/main/R)
R functions grouped by modules.  
This directory contains reusable functions implementing different processing, analysis, and correction steps used throughout the workflows.

---

### [data](https://github.com/tecladuran/targetml-workflows/tree/main/data)
Generated and input tables.  
This directory stores data tables produced by the workflows as well as auxiliary datasets required for analysis and reporting.

---

### [docs](https://github.com/tecladuran/targetml-workflows/tree/main/docs)
Reports with explanations and results.  
This directory contains methodological reports, analyses, and documented results describing the workflows and their outcomes.

---

## Usage

1. Explore the modules in the `R/` directory to identify available functions.
2. Use the tables in `data/` as inputs or outputs of the workflows.
3. Refer to the documents in `docs/` for detailed explanations and interpretation of results.

All scripts source the functions in `R/` and use the data stored in `data/`.
All paths are defined as relative paths, so the workflows can be executed directly after cloning the repository without requiring any path modification.
---

## Author

This repository was created by Tecla Duran Fort.

