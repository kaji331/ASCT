module AutomaticSingleCellToolbox

    using Arpack,TSVD
    using MultivariateStats,HDF5,GZip,GLM,Loess,DataFrames,CairoMakie,Dates
    using ColorSchemes,Peaks,Clustering,Random,StatsBase
    using Statistics,RuntimeGeneratedFunctions,Distances,SparseArrays,TSne,UMAP
    using LinearAlgebra,Graphs,SimpleWeightedGraphs,NearestNeighbors
    using HypothesisTests,MultipleTesting,NaturalSort
    import KernelDensity

    export Read10X
    export MergeRawData
    export FeaturePercentage!,ManualFilter!,AutoFilter!
    export NormalizeData!,SelectHVG!,RegressObs!,FeatureScore!
    export PCA!,TSNE!,UMAP!,FastRowScale!
    export Clustering!
    export DE!
    export Harmony!
    export SaveSeuratV4,SaveAnnData
    export FeatureHeat,DimensionPoints,DrawQC,FeatureVariances,ElbowPCA
    export FeatureViolin,FeatureFracDots

    if Threads.nthreads() != 1
        throw(ArgumentError("Only support one thread!"))
    end

    include("DataObject.jl")
    include("DataFetch.jl")
    include("DataManipulation.jl")
    include("QualityCheck.jl")
    include("DataPreprocess.jl")
    include("Reductions.jl")
    include("SNN.jl")
    include("Louvain.jl")
    include("ModularityClustering.jl")
    include("Clustering.jl")
    include("DE.jl")
    include("Harmony.jl")
    include("Visualization.jl")
    include("Utils.jl")

    # https://www.science.org/doi/abs/10.1126/science.aad0501
    # 2019 updated symbols
    const s_genes,g2m_genes = ["MCM5","PCNA","TYMS","FEN1","MCM7","MCM4","RRM1",
                               "UNG","GINS2","MCM6","CDCA7","DTL","PRIM1",
                               "UHRF1","CENPU","HELLS","RFC2","POLR1B","NASP",
                               "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7",
                               "POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45",
                               "CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
                               "USP1","CLSPN","POLA1","CHAF1B","MRPL36","E2F8"],
                              ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2",
                               "TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67",
                               "TMPO","CENPF","TACC3","PIMREG","SMC4","CCNB2",
                               "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E",
                               "TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","JPT1",
                               "CDC20","TTK","CDC25C","KIF2C","RANGAP1",
                               "NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23",
                               "HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5",
                               "CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5",
                               "CENPA"]

end # module
