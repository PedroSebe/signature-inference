# signature-inference
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.1-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

> This repo is an attempt of recovering mutational signatures information from small panel sequencing data. Still very in-progress!

## 🧬 What are mutational signatures?

The mutations that underlie cancer are caused by many factors, ranging from tobacco smoking to spontaneous chemical reactions in the DNA. These factors leave detectable patterns known as [mutational signatures](https://www.nature.com/articles/s41568-021-00377-7), that can be found using whole-genome or whole-exome sequencing.

## ⚕️ Why is this important?
**Mutational signature analysis can help cancer patients get more effective treatment**, by a strategy known as synthetic lethality. For instance, patients with signature 5 usually respond better to therapy based on platinum compounds or PARP inhibitors, while patients with signature 6 may respond better to immunotherapy. However, expensive whole-genome sequencing is usually necessary for this analysis. This project is an attempt to perform similar analyses for smaller and cheaper panel data, which is more accessible in the clinical setting.

## 💡 Data source
**The core of this project is the whole-genome data available from the consortium [Pan-Cancer Analysis of Whole Genomes (PCAWG)](https://dcc.icgc.org/pcawg).** Using this data, we can understand which combinations of signatures are more likely for each primary site and also simulate small panel data by simple subsetting, in order to relate panel-based mutational catalogues to whole-genome catalogues.

## 💻 Reprodutibility
To reproduce these analyses, you will need:
- `Snakemake>=7.0`;
- `conda>=4.10` ;
> Do not try to reproduce this analysis yet! This is a work in progress, so the main analysis is not yet runnable.

## 📝 License
[MIT](https://choosealicense.com/licenses/mit/)