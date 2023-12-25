# Automatic Single-cell Toolbox (ASCT)

## Targets

- [x] Analyzing single-cell data in Julia, like Seurat v4 in R.
- [x] Automatically choosing parameters for each function to implement auto pipeline.
- [ ] Add more advanced funtions for dimensinal reduction, pseudotime analysis etc.
- [ ] Data exchange to other tools.
- [ ] Large-scaling by Julia's performance/data structure optimization.
- [ ] Web-based app for single-cell beginners and wet-lab users.

## Initial structure

For this version, We use a mutable struct for a single-cell data set. It stores 
UMI data as sparse array and barcodes/features as DataFrame. The struct also 
saves all transformed data, meta information and operation logs.

1. AutomaticSingleCellToolbox.jl: Module entrance.
2. DataObject.jl: Definition the single-cell data struct.
3. DataFetch.jl: Definition the functions for single-cell data reading.
4. QualityCheck.jl: Definition the functions of QC filtering and visualization.
5. DataPreprocess.jl: Definition the functions of data preprocessing.
6. Reductions.jl: Definition the functions of PCA, tSNE, UMAP etc.
7. Clustering.jl/ModularityClustering.jl/SNN.jl/Louvain.jl: Definition the function for clustering.
8. DE.jl: Definition the functions for marker genes identification.
9. DataManipulation.jl: Definition the function for data manipulation.
10. Harmony.jl: Reimplementation of harmony algorithm for data integration. 

## Tutorials
