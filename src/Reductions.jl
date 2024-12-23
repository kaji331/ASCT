"""
    PCA!(obj)

Add the PCA results into WsObj.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `max_pc::Integer = 50`: max principal components would be calculated.
- `use_hvg::Bool = true`: only use highly variable features or not.
- `method::AbstractString = "PCA"`: "Arpack" or "TSVD" might be better for large 
  data.
- `cut::Union{Integer,Symbol} = :LinearRegression`: only 1:cut PCs will be used 
  downstream automatically.
- `seed::Integer = -1`: random seed.
- `ptr::AbstractString = :auto`: calculation on raw data.
"""
function PCA!(obj::WsObj;
        max_pc::Integer = 50,
        use_hvg::Bool = true,
        method::AbstractString = "PCA",
        cut::Union{Integer,Symbol} = :LinearRegression,
        seed::Integer = -1,
        ptr::Union{Symbol,AbstractString} = :auto)

    # Check parameters
    if max_pc < 3 || max_pc > 100
        return "Nothing to do! 'max_pc' is Only support 3 ~ 100 !"
    end

    # Scaling
    if "scale_dat" in keys(obj.dat)
        @info "Using pre-calculated scale data..."
        dat = obj.dat["scale_dat"]
    else
        @info "Scaling data..."
        FastRowScale!(obj;ptr=ptr)
        dat = obj.dat["scale_dat"]
    end

    # High variable genes
    if use_hvg && size(dat,1) > length(obj.meta["hvg_index"])
        dat = @view dat[obj.meta["hvg_index"],:]
    end

    # PCA
    @info "Running PCA ..."
    if !(method in ["PCA","Arpack","TSVD"])
        return "'method' is only support PCA, Arpack or TSVD!"
    end
    if method == "PCA"
        M = MultivariateStats.fit(PCA,dat;maxoutdim=max_pc)
        Y = MultivariateStats.transform(M,dat)' |> Matrix
        var = principalvars(M)
    elseif method == "Arpack"
        # Original matrix m (rows are observation, columns are features) is: 
        #   U * diagm(S) * V'
        # or 
        #   U * Diagonal(S) * V' (for sparse matrix)

        Random.seed!(seed < 0 ? 1984 : seed)
        U,S,V = svds(dat';nsv=max_pc)[1]
        var = S .^ 2
    else
        Random.seed!(seed < 0 ? 1984 : seed)
        U,S,V = tsvd(dat',max_pc)
        var = S .^ 2
    end

    if cut == :LinearRegression
        # Better results for current tests
        @info "Looking for Elbow threshold..."
        cut = FindElbowPointLM(var)
        @info "We recommend top $cut PCs for downstream analysis automatically!"
    elseif cut == :ParallelAnalysis
        # Better performance
        @info "Looking for Elbow threshold..."
        cut = FindElbowPoint(var)
        @info "We recommend top $cut PCs for downstream analysis automatically!"
    elseif typeof(cut) <: Symbol
        return "Nothing to do! 'cut' is only :LinearRegression, " * 
            ":ParallelAnalysis or an positive integer!"
    end

    # Output
    if method == "PCA"
        obj.meta["pca"] = Y
    else
        # https://medium.com/machine-learning-world/linear-algebra-svd-and-pca-5979f739e95a
        obj.meta["pca"] = U * diagm(S)
    end
    # Variance explained, not proportion explained.
    obj.meta["pca_var"] = var ./ sum(var)
    obj.meta["pca_cut"] = cut

    # log
    push!(obj.log,String("PCA!(WsObj;max_pc=$max_pc,use_hvg=$use_hvg," * 
                         "method=$method,seed=$seed,cut=$cut,ptr=$ptr)"))

    return "Finished!"
end

"""
    UMAP!(obj)

Add the UMAP results into WsObj.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `reduction::Union{AbstractString,Symbol} = :auto`: if there's a harmony in the 
  meta of WsObj, it will use harmony data or pca data.
- `use_pc::Union{AbstractString,Integer} = "pca_cut"`: a number of top PCs for 
  calculation or the automatical threshold from PCA! function.
- `map_dims::Integer = 2`: max dimensions of UMAP result.
- `seed::Integer = -1`: random seed.
- `n_neighbors::Integer = 30`: number of k nearest neighbors.
- `metric::SemiMetric = CosineDist()`: the function of distance calculation.
- `min_dist::Real = 0.5`: control the distances among points.
- `n_epochs::Integer = 300`: the number of training epochs.
- `learning_rate::Real = 1`: the initial learning rate.
- `init::Symbol = :spectral`: how to initialize the output embedding. :spectral 
  or :random.
- `spread::Real = 1`: the effective scale.
- `set_operation_ratio::Real = 1`: between 0 and 1. 1 means pure fuzzy union, 
  and 0 means pure fuzzy intersection.
- `local_connectivity::Integer = 1`: the number by assumed local connection of 
  nearest neighbors.
- `repulsion_strength::Real = 1`: the weights of negative samples.
- `neg_sample_rate::Integer = 5`: higer value will increase computational cost 
  but slightly more accuracy.
"""
function UMAP!(obj::WsObj;
        reduction::Union{AbstractString,Symbol} = :auto,
        use_pc::Union{AbstractString,Integer} = "pca_cut",
        map_dims::Integer = 2,
        seed::Integer = -1,
        n_neighbors::Integer = 30,
        metric::SemiMetric = CosineDist(),
        min_dist::Real = 0.5,
        n_epochs::Integer = 300,
        learning_rate::Real = 1,
        init::Symbol = :spectral,
        spread::Real = 2,
        set_operation_ratio::Real = 1,
        local_connectivity::Integer = 1,
        repulsion_strength::Real = 1,
        neg_sample_rate::Integer = 5)

    # Support PCA or Harmony
    if reduction == :auto
        if "harmony" in keys(obj.meta)
            reduction = "harmony"
        else
            reduction = "pca"
        end
    end

    # Check and prepare data
    if typeof(use_pc) <: AbstractString
        use_pc = obj.meta[use_pc]
    end
    if use_pc <= 2
        return "Nothing to do! The 'use_pc' must be larger than 2!"
    end
    if reduction == "pca"
        X = transpose(obj.meta["pca"][:,1:use_pc])
    else
        X = obj.meta[reduction]'
    end

    # Run
    Random.seed!(seed != -1 ? seed : 1984)
    @info "Running UMAP..."
    y = umap(X,map_dims;
             n_neighbors = n_neighbors,
             metric = metric,
             min_dist = min_dist,
             n_epochs = n_epochs,
             learning_rate = learning_rate,
             init = init,
             spread = spread,
             set_operation_ratio = set_operation_ratio,
             local_connectivity = local_connectivity,
             repulsion_strength = repulsion_strength,
             neg_sample_rate = neg_sample_rate)' |> Matrix

    # Update WsObj
    obj.meta["umap"] = y

    # log
    push!(obj.log,String("UMAP!(WsObj;" * 
                         "reduction=$reduction," * 
                         "use_pc=$use_pc," * 
                         "map_dims=$map_dims," * 
                         "seed=$seed," * 
                         "n_neighbors=$n_neighbors," * 
                         "metric=$metric," * 
                         "min_dist=$min_dist," * 
                         "n_epochs=$n_epochs," * 
                         "learning_rate=$learning_rate," * 
                         "init=$init," * 
                         "spread=$spread," * 
                         "set_operation_ratio=$set_operation_ratio," * 
                         "local_connectivity=$local_connectivity," * 
                         "repulsion_strength=$repulsion_strength," * 
                         "neg_sample_rate=$neg_sample_rate)"))

    return "Finished!"
end

function FastTSNE(X::AbstractMatrix;
        use_pc::Union{AbstractString,Integer} = "pca_cut",
        theta::AbstractFloat = 0.5,
        perplexity::Real = 30,
        map_dims::Integer = 2,
        max_iter::Integer = 750,
        stop_early_exag_iter::Integer = 250,
        K::Integer = -1,
        sigma::Real = -1,
        nbody_algo::String = "FFT",
        knn_algo::String = "annoy",
        mom_switch_iter::Integer = 250,
        momentum::AbstractFloat = 0.5,
        final_momentum::AbstractFloat = 0.8,
        learning_rate::Union{Symbol,AbstractFloat} = :auto,
        early_exag_coeff::Real = 12,
        no_momentum_during_exag::Bool = false,
        n_trees::Integer = 50,
        search_k::Union{Symbol,Integer} = :auto,
        start_late_exag_iter::Union{Symbol,Integer} = :auto,
        late_exag_coeff::Real = -1,
        nterms::Integer = 3,
        intervals_per_integer::Real = 1,
        min_num_intervals::Integer = 50,
        seed::Integer = -1,
        initialization::Union{AbstractString,
                              Matrix{T} where T <: AbstractFloat} = "pca",
        load_affinities::Union{Nothing,AbstractString} = nothing,
        perplexity_list::Union{Nothing,
                               Vector{T} where T <: AbstractFloat} = nothing,
        df::Real = 1,
        return_loss::Bool = false,
        nthreads::Integer = -1,
        max_step_norm::Union{Nothing,Integer} = 5)

    # Check fast_tsne
    if Sys.iswindows()
        throw(DomainError("operation system","Does NOT support Windows!"))
    end
    check_fast_tsne = try 
            ignorestatus(`fast_tsne`) |> readchomp
        catch e
            @warn "Please install from https://github.com/KlugerLab/FIt-SNE"
            throw(e)
        end
    version_number = match(r"v(.*) ===",check_fast_tsne)[1]

    if learning_rate == :auto
        learning_rate = maximum([200,size(X,1) / early_exag_coeff])
    end

    if start_late_exag_iter == :auto
        if late_exag_coeff > 0
            start_late_exag_iter = stop_early_exag_iter
        else
            start_late_exag_iter = -1
        end
    end

    if initialization == "pca"
        if !("pca" in keys(obj.meta))
            throw(DomainError(initialization,"Run PCA! function first..."))
        else
            pca_a = obj.meta["pca"][:,1:map_dims]
            pca_l = obj.meta["pca_var"][1:map_dims]
        end
        initialization = sqrt.(pca_l) |> 
            x -> diagm(map_dims,map_dims,x) |> 
            x -> pca_a * x |> 
            x -> x ./ std(x[:,1]) .* 0.0001
    end

    if !isnothing(perplexity_list)
        perplexity = 0
    end

    if sigma > 0 && K > 0
        perplexity = -1
    end

    if search_k == :auto
        if perplexity > 0
            search_k = 3 * perplexity * n_trees
        elseif perplexity == 0
            search_k = 3 * maximum(perplexity_list) * n_trees
        else
            search_k = K * n_trees
        end
    end

    if nbody_algo == "Barnes-Hut"
        nbody_algo = 1
    elseif nbody_algo == "FFT"
        nbody_algo = 2
    else
        throw(DomainError(nbody_algo,
                          "Only 'Barnes-Hut' or 'FFT' are supported for " * 
                          "nbody_algo!"))
    end

    if knn_algo == "annoy"
        knn_algo = 1
    elseif knn_algo == "vp-tree"
        knn_algo = 2
    else
        throw(DomainError(knn_algo,"Only 'annoy' or 'vp-tree' are supported " * 
                          "for knn_algo!"))
    end

    if load_affinities == "load"
        load_affinities = 1
    elseif load_affinities == "save"
        load_affinities= 2
    else
        load_affinities= 0
    end

    if nthreads <= 0
        nthreads = 0
    end

    if isnothing(max_step_norm)
        max_step_norm = -1
    end

    if no_momentum_during_exag
        no_momentum_during_exag = 1
    else
        no_momentum_during_exag = 0
    end

    # Prepare input file and read output file for tsne
    stamp = string(rand(Random.RandomDevice(),1:99999))
    # tn = tempname()
    tn = "."
    input_dat = tn * "/tsne_in_" * Dates.format(now(),"yyyy.m.d.H.M.S") * "-" * 
         stamp * ".dat"
    output_dat = tn * "/tsne_out_" * Dates.format(now(),"yyyy.m.d.H.M.S") * 
        "-" * stamp * ".dat"

    # Write input file
    open(input_dat,"w") do f
        write(f,Int32(size(X,1)))
        write(f,Int32(size(X,2)))
        write(f,Float64(theta))
        write(f,Float64(perplexity))
        if perplexity == 0
            write(f,Int32(length(perplexity_list)))
            for perpl in perplexity_list
                write(f,Float64(perpl))
            end
        end
        write(f,Int32(map_dims))
        write(f,Int32(max_iter))
        write(f,Int32(stop_early_exag_iter))
        write(f,Int32(mom_switch_iter))
        write(f,Float64(momentum))
        write(f,Float64(final_momentum))
        write(f,Float64(learning_rate))
        write(f,Float64(max_step_norm))
        write(f,Int32(K))
        write(f,Float64(sigma))
        write(f,Int32(nbody_algo))
        write(f,Int32(knn_algo))
        write(f,Float64(early_exag_coeff))
        write(f,Int32(no_momentum_during_exag))
        write(f,Int32(n_trees))
        write(f,Int32(search_k))
        write(f,Int32(start_late_exag_iter))
        write(f,Float64(late_exag_coeff))
        write(f,Int32(nterms))
        write(f,Float64(intervals_per_integer))
        write(f,Int32(min_num_intervals))
        write(f,reinterpret(UInt8,vec(X'))) # X is from the obj
        write(f,Int32(seed))
        write(f,Float64(df))
        write(f,Int32(load_affinities))
        if initialization != "random"
            write(f,reinterpret(UInt8,vec(initialization')))
        end
    end

    # Run fast tsne
    run(`fast_tsne $version_number $input_dat $output_dat $nthreads`)

    # Read output file
    row::Int32 = 0
    col::Int32 = 0
    dat = Matrix{Float64}
    loss::Vector{Float64} = []

    open(output_dat,"r") do f
        col = read(f,Int32)
        row = read(f,Int32)
        dat = read(f,row * col * 8) |> x -> reinterpret(Float64,x) |> 
            x -> reshape(Vector(x),Int(row),Int(col))'
        if return_loss
            read(f,Int32)
            loss = read(f,max_iter * 8) |> x -> reinterpret(Float64,x) |> Vector
        end
    end

    rm(input_dat;force=true)
    rm(output_dat;force=true)

    if return_loss
        return (dat,loss)
    else
        return dat
    end
end

"""
    TSNE!(obj)

Add the t-SNE results into WsObj.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `reduction::Union{AbstractString,Symbol} = :auto`: if there's a harmony in the 
  meta of WsObj, it will use harmony data or pca data.
- `use_pc::Union{AbstractString,Integer} = "pca_cut"`: a number of top PCs for 
  calculation or the automatical threshold from PCA! function.
- `fast::Bool = true`: use FIt-SNE or not.
- `theta::AbstractFloat = 0.5`: FIt-SNE precision. 0 for exact.
- `perplexity::Real = 30`: determine the bandwidth of the Gaussian kernel.
- `map_dims::Integer = 2`: max dimensions of t-SNE result.
- `max_iter::Integer = 750`: number of iterations.
- `stop_early_exag_iter::Integer = 250`: when to switch off early exaggeration.
- `K::Integer = -1`: number of k nearest neighbors.
- `sigma::Real = -1`: fixed sigma when perplexity is -1.
- `nbody_algo::String = "FFT"`: "FFT" or "Barnes-Hut".
- `knn_algo::String = "annoy"`: "annoy" or "vp-tree".
- `mom_switch_iter::Integer = 250`: iteration number from momentum to 
  final_momentum.
- `momentum::AbstractFloat = 0.5`: initial momentum.
- `final_momentum::AbstractFloat = 0.8`: momentum used for later optimization.
- `learning_rate::Union{Symbol,AbstractFloat} = :auto`: :auto or a float number.
- `early_exag_coeff::Real = 12`: coefficient for early exaggeration.
- `no_momentum_during_exag::Bool = false`: swith off momentum during the early 
  exaggeration.
- `n_trees::Integer = 50`: number of search trees for annoy.
- `search_k::Union{Symbol,Integer} = :auto`: number of nodes to inspect during 
  search for annoy.
- `start_late_exag_iter::Union{Symbol,Integer} = :auto`: when to start late 
  exaggeration.
- `late_exag_coeff::Real = -1`: coefficient for late exaggeration.
- `nterms::Integer = 3`: number of interpolation points per sub-interval.
- `seed::Integer = -1`: random seed.
- `initialization::Union{AbstractString,Matrix{<: AbstractFloat}} = "pca"`: 
  the matrix of initialization.
- `load_affinities::Union{Nothing,AbstractString} = nothing`: "save" or "load" 
  the input similarities.
- `perplexity_list::Union{Nothing,Vector{<: AbstractFloat}} = nothing`: a list 
  of perplexities to used as a perplexity combination.
- `df::Real = 1`: controls the degree of freedom of t-distribution.
- `nthreads::Integer = -1`: number of threads. -1 for all available threads.
- `max_step_norm::Union{Nothing,Integer} = 5`: max distance that a point is 
  allowed to move on one iteration.
- `verbose::Bool = false`: output informational and diagnostic messages for 
  TSne.jl.
- `progress::Bool = true`: display progress meter during t-SNE optimization.
- `extended_output::Bool = false`: return a tuple of embedded coordinates 
  matrix, point perplexities and final Kullback-Leibler divergence or not.
"""
function TSNE!(obj::WsObj;
        reduction::Union{AbstractString,Symbol} = :auto,
        use_pc::Union{AbstractString,Integer} = "pca_cut",
        fast::Bool = true,
        theta::AbstractFloat = 0.5,
        perplexity::Real = 30,
        map_dims::Integer = 2,
        max_iter::Integer = 750,
        stop_early_exag_iter::Integer = 250,
        K::Integer = -1,
        sigma::Real = -1,
        nbody_algo::String = "FFT",
        knn_algo::String = "annoy",
        mom_switch_iter::Integer = 250,
        momentum::AbstractFloat = 0.5,
        final_momentum::AbstractFloat = 0.8,
        learning_rate::Union{Symbol,AbstractFloat} = :auto,
        early_exag_coeff::Real = 12,
        no_momentum_during_exag::Bool = false,
        n_trees::Integer = 50,
        search_k::Union{Symbol,Integer} = :auto,
        start_late_exag_iter::Union{Symbol,Integer} = :auto,
        late_exag_coeff::Real = -1,
        nterms::Integer = 3,
        intervals_per_integer::Real = 1,
        min_num_intervals::Integer = 50,
        seed::Integer = -1,
        initialization::Union{AbstractString,
                              Matrix{T} where T <: AbstractFloat} = "pca",
        load_affinities::Union{Nothing,AbstractString} = nothing,
        perplexity_list::Union{Nothing,
                               Vector{T} where T <: AbstractFloat} = nothing,
        df::Real = 1,
        nthreads::Integer = -1,
        max_step_norm::Union{Nothing,Integer} = 5,
        verbose::Bool = false,
        progress::Bool = true,
        min_gain::AbstractFloat = 0.01,
        extended_output::Bool = false)

    # Support PCA or Harmony
    if reduction == :auto
        if "harmony" in keys(obj.meta)
            reduction = "harmony"
        else
            reduction = "pca"
        end
    end

    # Choose PCs
    if typeof(use_pc) <: AbstractString
        use_pc = obj.meta[use_pc]
    end
    if use_pc <= 2
        return "The 'use_pc' must be larger than 2!"
    end
    if reduction == "pca"
        X = obj.meta["pca"][:,1:use_pc]
        initialization = "random"
    else
        X = obj.meta[reduction]
        initialization = "random"
    end

    # Check algorithm
    if fast
        @info "Running FIt-SNE ..."
        y = FastTSNE(X;
                     use_pc = use_pc,
                     theta = theta,
                     perplexity = perplexity,
                     map_dims = map_dims,
                     max_iter = max_iter,
                     stop_early_exag_iter = stop_early_exag_iter,
                     K = K,
                     sigma = sigma,
                     nbody_algo = nbody_algo,
                     knn_algo = knn_algo,
                     mom_switch_iter = mom_switch_iter,
                     momentum = momentum,
                     final_momentum = final_momentum,
                     learning_rate = learning_rate,
                     early_exag_coeff = early_exag_coeff,
                     no_momentum_during_exag = no_momentum_during_exag,
                     n_trees = n_trees,
                     search_k = search_k,
                     start_late_exag_iter = start_late_exag_iter,
                     late_exag_coeff = late_exag_coeff,
                     nterms = nterms,
                     intervals_per_integer = intervals_per_integer,
                     min_num_intervals = min_num_intervals,
                     seed = seed,
                     initialization = initialization,
                     load_affinities = load_affinities,
                     perplexity_list = perplexity_list,
                     df = df,
                     return_loss = extended_output,
                     nthreads = nthreads,
                     max_step_norm = max_step_norm)

    else
        Random.seed!(seed != -1 ? seed : 1984)
        @info "Running t-SNE ..."
        # -1 means 'Do NOT run PCA' in the 'tsne' function internally
        y = tsne(X,map_dims,-1,max_iter,perplexity;
                 pca_init = initialization == "pca" ? true : false,
                 verbose = verbose,
                 progress = progress,
                 min_gain = min_gain,
                 eta = maximum([200,size(X,1) / early_exag_coeff]),
                 initial_momentum = momentum,
                 final_momentum = final_momentum,
                 momentum_switch_iter = mom_switch_iter,
                 stop_cheat_iter = stop_early_exag_iter,
                 cheat_scale = early_exag_coeff,
                 extended_output = extended_output)
    end

    # Update WsObj
    if extended_output
        obj.meta["tsne"] = y[1] |> Matrix
        if fast
            obj.meta["tsne_loss"] = y[2]
        else
            obj.meta["tsne_point_perplexity"] = y[2]
        end
    else
        obj.meta["tsne"] = y |> Matrix
    end
    
    # log
    push!(obj.log,String("TSNE!(WsObj;" * 
                         "reduction=$reduction," * 
                         "use_pc=$use_pc," * 
                         "fast=$fast," * 
                         "theta=$theta," * 
                         "perplexity=$perplexity," * 
                         "map_dims=$map_dims," * 
                         "max_iter=$max_iter," * 
                         "stop_early_exag_iter=$stop_early_exag_iter," * 
                         "K=$K," * 
                         "sigma=$sigma," * 
                         "nbody_algo=$nbody_algo," * 
                         "knn_algo=$knn_algo," * 
                         "mom_switch_iter=$mom_switch_iter," * 
                         "momentum=$momentum," * 
                         "final_momentum=$final_momentum," * 
                         "learning_rate=$learning_rate," * 
                         "early_exag_coeff=$early_exag_coeff," * 
                         "no_momentum_during_exag=$no_momentum_during_exag," * 
                         "n_trees=$n_trees," * 
                         "search_k=$search_k," * 
                         "start_late_exag_iter=$start_late_exag_iter," * 
                         "late_exag_coeff=$late_exag_coeff," * 
                         "nterms=$nterms," * 
                         "intervals_per_integer=$intervals_per_integer," * 
                         "min_num_intervals=$min_num_intervals," * 
                         "seed=$seed," * 
                         "initialization=$initialization," * 
                         "load_affinities=$load_affinities," * 
                         "perplexity_list=$perplexity_list," * 
                         "df=$df," * 
                         "nthreads=$nthreads," * 
                         "max_step_norm=$max_step_norm," * 
                         "verbose=$verbose," * 
                         "progress=$progress," * 
                         "min_gain=$min_gain," * 
                         "extended_output=$extended_output)"))

    return "Finished!"
end

"""
    FastRowScale!(obj)

Add the scale data into WsObj's dat.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `center::Bool = true`: centralize the data or not.
- `scale::Bool = true`: scale the data or not
- `scale_max::Real = 10`: 10 is the default value, like Seurat/Scanpy.
- `ptr::Union{Symbol,AbstractString} = :auto`: use "norm\\_dat", "regress\\_dat" 
  or :auto.
"""
function FastRowScale!(obj::WsObj;
        center::Bool = true,
        scale::Bool = true,
        scale_max::Real = 10,
        ptr::Union{Symbol,AbstractString} = :auto)

    if ptr == :auto && "regress_dat" in keys(obj.dat)
        ptr = "regress"
    elseif ptr == :auto && !("regress_dat" in keys(obj.dat))
        ptr = "norm"
    elseif typeof(ptr) <: Symbol && ptr != :auto
        return "Nothing to do! 'ptr' needs :auto or a specific name of dat!"
    end

    mat = copy(obj.dat[ptr * "_dat"])' |> sparse
    if center
        row_mean = mean(mat;dims=1)'
    end
    if scale
        if center
            row_sd = std(mat;dims=1)'
        else
            row_sd = sqrt.((mat .^ 2)' * ones(Int32,mat.m) ./ (mat.m - 1))
        end
    end
    mat = Matrix{Float64}(mat')
    if center
        mat .-= row_mean
    end
    if scale
        mat ./= row_sd
    end
    if !isinf(scale_max)
        mat[mat .> scale_max] .= scale_max
    end
    mat[isnan.(mat)] .= 0

    # Output
    obj.dat["scale_dat"] = mat
    # log
    push!(obj.log,String("FastRowScale!(WsObj;center=$center,scale=$scale" * 
                         "scale_max=$scale_max,ptr=$ptr)"))

    return "Finished!"
end

function FastRowScale(m::AbstractSparseMatrix;
        center::Bool = true,
        scale::Bool = true,
        scale_max::Real = 10)

    mat = copy(m)' |> sparse
    if center
        row_mean = mean(mat;dims=1)'
    end
    if scale
        if center
            row_sd = std(mat;dims=1)'
        else
            row_sd = sqrt.((mat .^ 2) * ones(Int32,mat.m) ./ (mat.m - 1))
        end
    end
    mat = Matrix{Float64}(mat')
    if center
        mat .-= row_mean
    end
    if scale
        mat ./= row_sd
    end
    if !isinf(scale_max)
        mat[mat .> scale_max] .= scale_max
    end
    mat[isnan.(mat)] .= 0

    # Scaled matrix
    return mat
end

# From kevinblighe/PCAtools for R language in github.com
# Use variances, not standard deviations
function FindElbowPoint(variances::Vector{<: AbstractFloat})
    lv = length(variances)
    dx = lv - 1
    dy = variances[end] - variances[1]
    l2 = sqrt(dx ^ 2 + dy ^ 2)
    dx = dx / l2
    dy = dy / l2

    dx0 = eachindex(variances) .- 1
    dy0 = variances .- variances[1]

    parallel_l2 = sqrt.((dx0 .* dx) .^ 2 .+ (dy0 .* dy) .^ 2)
    normal_x = dx0 .- dx .* parallel_l2
    normal_y = dy0 .- dy .* parallel_l2
    normal_l2 = sqrt.(normal_x .^ 2 .+ normal_y .^ 2)

    below_line = (normal_x .< 0) .& (normal_y .< 0)
    if !any(below_line)
        return lv
    else
        cut = findall(below_line)[findmax(normal_l2[below_line])[2]] + 3
        if cut <= lv
            return cut
        else
            return lv
        end
    end
end

function FindElbowPointLM(variances::Vector{<: AbstractFloat})
    value = log10.(variances .^ 2 / sum(variances .^ 2))
    lv = length(value)
    dif = DataFrame(axis_x=[1,lv],axis_y=[value[1],value[end]]) |> 
        x -> lm(@formula(axis_y ~ axis_x),x) |> 
        x -> GLM.predict(x,DataFrame(axis_x=1:lv)) |> 
        x -> x .- value
    cut = findmax(dif)[2] + 1
    if cut <= lv
        return cut
    else
        return lv
    end
end
