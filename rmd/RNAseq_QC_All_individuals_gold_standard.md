---
title: "RNA-seq QC - All individuals"
author: "Quentin Miagoux"
date: "2026-01-23"
output:
  html_document:
    df_print: paged
    theme: united
    highlight: tango
    keep_md: yes
    number_sections: yes
    toc: yes
    toc_float: yes
    css: "../my_css.css"
params:
  expr_xlsx: ""
  coldata_xlsx: ""
  top_variable_genes: 500
editor_options:
  chunk_output_type: console
---











# QC {.tabset .tabset-fade}



## Gold standard QC metrics {.tabset .tabset-fade}

### Basic sample metrics (counts)

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/qc-sample-metrics-1.png)<!-- -->![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/qc-sample-metrics-2.png)<!-- -->![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/qc-sample-metrics-3.png)<!-- -->

### RLD distributions

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/qc-rld-distribution-1.png)<!-- -->![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/qc-rld-distribution-2.png)<!-- -->

### Mitochondrial fraction (if MT- / mt- genes present)


```
## No mitochondrial genes detected by prefix ^MT- or ^mt- in gene symbols.
```

### DESeq2 outlier diagnostic (Cook's distance)

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/deseq2-diagnostics-1.png)<!-- -->

## All individuals {.tabset .tabset-fade}

### All genes

#### PCA

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/pca-all-genes-1.png)<!-- -->

#### Sample correlation heatmap 

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-sample-correlation-1.png)<!-- -->

#### Sample distance heatmap 

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-sample-distance-1.png)<!-- -->

### Top 500 variable genes

#### PCA

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/pca-top-genes-1.png)<!-- -->

#### Heatmap

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-top-genes-1.png)<!-- -->

## Without WT3 and DMD9 {.tabset .tabset-fade}

### All genes

#### PCA

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/pca-all-genes-no-wt3-dmd9-1.png)<!-- -->

#### Sample correlation heatmap (filtered)

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-sample-correlation-filtered-1.png)<!-- -->

#### Sample distance heatmap (filtered)

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-sample-distance-filtered-1.png)<!-- -->

### Top 500 variable genes

#### PCA

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/pca-top-genes-no-wt3-dmd9-1.png)<!-- -->

#### Heatmap

![](RNAseq_QC_All_individuals_gold_standard_files/figure-html/heatmap-top-genes-no-wt3-dmd9-1.png)<!-- -->
