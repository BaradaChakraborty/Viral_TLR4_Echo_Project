# The Viral TLR4/9 Mitochondrial Echo in Parkinson's Disease

**Computational mapping of RNA virus-induced Rps3-mtDNA extrusion and microglial TLR4/MyD88 hyperactivation in neurodegeneration.**

## Project Overview
Severe RNA viral infections hijack host cellular metabolism, inducing profound mitochondrial fragmentation and oxidative stress. Recent literature demonstrates that heavily oxidized mitochondrial DNA (mtDNA) is recognized by Ribosomal Protein S3 (Rps3) and extruded from neurons via extracellular vesicles. 

This repository contains an independent *in-silico* systems biology pipeline designed to test the **"Viral Echo Hypothesis."** The pipeline computationally investigates whether the transcriptomic footprint of RNA virus-induced mitochondrial stress perfectly aligns with the hyper-inflammatory, MyD88-driven state of Disease-Associated Microglia (DAM) in Parkinson's disease.

## Computational Strategy
By integrating viral infection transcriptomics with human post-mortem Parkinson's disease single-nucleus RNA-sequencing (snRNA-seq) datasets (e.g., GSE140231), this project executes:
* **Custom Module Scoring:** Projecting virally-derived Rps3-mtDNA extrusion signatures onto PD microglial clusters using Seurat.
* **Transcription Factor Inference:** Utilizing `decoupleR` to mathematically decouple the microglial inflammatory regulon, demonstrating that the pathology is driven by the TLR4/9-MyD88 axis (NF-κB/AP-1) rather than the classic cGAS-STING (IRF3) cytosolic DNA-sensing pathway.
