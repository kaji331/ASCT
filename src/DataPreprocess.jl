function NormalizeData!(obj::WsObj;
        ptr::AbstractString = "raw")

    dat = copy(obj.dat[ptr * "_dat"]) |> x -> Float32.(x)

    # Fast way for SparseArrays
    for col in 1:dat.n
        cs = sum(dat.nzval[dat.colptr[col]:dat.colptr[col + 1] - 1])
        dat.nzval[dat.colptr[col]:dat.colptr[col + 1] - 1] .= 
            log1p.(dat.nzval[dat.colptr[col]:dat.colptr[col + 1] - 1] ./ 
                   cs * 1e4)
    end

    # Output
    obj.dat["norm_dat"] = dat

    # log
    push!(obj.log,String("NormalizeData!(WsObj;ptr=$ptr)"))

    return "Finished!"
end

function SelectHVG!(obj::WsObj;
        hvg_number::Union{Symbol,Integer} = 2000,
        ptr::AbstractString = "raw",
        method::Symbol = :simple)

    dat = obj.dat[ptr * "_dat"]' |> sparse
    clip_max = sqrt(dat.m)
    # 'mean' is faster for columns, slower for rows to matrix multiplication
    avg = mean(dat;dims=1)' 
    variances = var(dat;dims=1)' |> vec |> Vector
    variances_exp = zeros(Float32,dat.n)
    variances_std = zeros(Float32,dat.n)
    idx_more_zero = variances .> 0
    x = log10.(avg[idx_more_zero])
    y = log10.(variances[idx_more_zero])
    # Locally regression, span is 0.3 in Seurat
    m = loess(x,y;span=0.3)
    variances_exp[idx_more_zero] = 10 .^ Loess.predict(m,x)
    std_exp = sqrt.(variances_exp)

    idx = [ i for i in 1:dat.n ][std_exp .!= 0]
    for col in idx
        s = sum([ minimum([clip_max,
                           (dat.nzval[i] - avg[col]) / std_exp[col]]) ^ 2 
                 for i in dat.colptr[col]:dat.colptr[col + 1] - 1 ])
        s += ((0 - avg[col]) / std_exp[col]) ^ 2 * 
            (dat.m - dat.colptr[col + 1] + dat.colptr[col])
        variances_std[col] = s / (dat.m - 1)
    end

    if hvg_number == :auto
        if method == :loess
            cutoff = FindCutoff1(variances_std)
        elseif method == :simple
            cutoff = FindCutoff2(avg,variances_std;method=method)
        elseif method == :dbscan
            cutoff = FindCutoff2(avg,variances_std;method=method)
        else
            @warn "Only ':simple', ':dbscan' and ':loess' were supported for " * 
                " 'method'! ':simple' is selected for wrong method!"
            cutoff = FindCutoff2(avg,variances_std;method=method)
        end
        idx = findall(x -> x > cutoff,variances_std)
    elseif typeof(hvg_number) <: Integer
        idx = findall(x -> x <= hvg_number,
                      sortperm(variances_std;rev=true) |> invperm)
    else
        @warn "Only support an integer or :auto for 'hvg_number'!" * 
            " Now, all genes were selected...long time running!"
        idx = 1:length(variances_std)
    end
    @info "$(length(idx)) HVGs were selected automatically!"

    # Output
    if isnothing(obj.meta)
        obj.meta = Dict("hvg_name" => obj.var.name[idx],
                        "hvg_index" => idx,
                        "hvg_mean" => avg,
                        "hvg_var_std" => variances_std)
    else
        obj.meta["hvg_name"] = obj.var.name[idx]
        obj.meta["hvg_index"] = idx
        obj.meta["hvg_mean"] = avg
        obj.meta["hvg_var_std"] = variances_std
    end

    # log
    push!(obj.log,String("SelectHVG!(WsObj;ptr=$ptr)"))

    return "Finished!"
end

function RegressObs!(obj::WsObj;
        var::AbstractString = "percentage_mt",
        use_hvg::Bool = true,
        ptr::AbstractString = "norm")

    if var in keys(obj.meta)
        if typeof(obj.meta[var]) == DataFrame
            latent_data = copy(obj.meta[var])
        else
            latent_data = DataFrame(obj.meta[var],:auto)
        end
    elseif var in names(obj.obs)
        latent_data = DataFrame()
        latent_data[!,var] = obj.obs[:,var]
    end

    if !@isdefined(latent_data)
        throw(DomainError(var,"No var found in WsObj!"))
    else
        model = term("y") ~ sum(term.(names(latent_data)))
    end

    if use_hvg
        dat = obj.dat[ptr * "_dat"][obj.meta["hvg_index"],:]' |> sparse
    else
        dat = obj.dat[ptr * "_dat"]' |> sparse
    end

    resid = Vector{Vector{Float32}}()
    for col in 1:dat.n
        latent_data.y = Vector(dat[:,col])
        push!(resid,
              lm(model,latent_data) |> residuals |> x -> x .-= minimum(x))
    end

    # dat
    obj.dat["regress_dat"] = log1p.(hcat(resid...)') |> sparse

    # log
    push!(obj.log,String("RegressMeta!(WsObj;var=$var,use_hvg=$use_hvg," * 
                         "ptr=$ptr)"))

    return "Finished!"
end

function FeatureScore!(obj::WsObj;
        features::Union{Nothing,
                        Dict{AbstractString,
                             Vector{AbstractString}}} = nothing,
        meta_name::AbstractString = "cell_cycle",
        seed::Integer = 1,
        ptr::AbstractString = "norm")

    # No features providing would calculate cell cycling scores
    if isnothing(features)
        features = Dict("s" => s_genes,"g2m" => g2m_genes)
        meta_name = "cell_cycle"
    end
    
    dat = obj.dat[ptr * "_dat"]
    gene_names = obj.var.name
    # Check genes
    tmp = Dict{eltype(features).types[1],Vector{Int64}}()
    for k in keys(features)
        idx = findall(x -> x in features[k],gene_names)
        if !isempty(idx)
            push!(tmp,k => idx)
        end
    end
    if !isempty(tmp)
        features = tmp
        tmp = nothing
    else
        return "Nothing to do! No enough features in the WsObj!"
    end

    # Generate control data
    avg_genes = dat * ones(Int32,dat.n) ./ dat.n
    GC.gc()
    order_avg = sortperm(avg_genes)
    gene_names_ordered = gene_names[order_avg]
    # Cut genes to regions, like Seurat
    width = floor(Int,dat.m / 24)
    seq = repeat(1:24;inner=width)
    if dat.m % 24 != 0
        index = 24
        @simd for i in 1:(dat.m % 24)
            if i % 2 == 1
                push!(seq,index)
                index -= 1
            else
                push!(seq,24 - index)
            end
        end
        seq = sort(seq)
    end

    scores = DataFrame()
    @inbounds for k in keys(features)
        d = dat[features[k],:]
        if typeof(d) <: AbstractSparseVector
            avg_feature_cell = d
        else
            avg_feature_cell = d' * ones(Int32,d.m) ./ d.m
        end
        ctrl_idx = Vector{Vector{Int}}()
        Random.seed!(seed < 0 ? 1984 : seed)
        @simd for j in features[k]
            idx = findfirst(x -> x == gene_names[j],
                            gene_names_ordered)
            idx = findall(x -> x == seq[idx],seq)
            idx = sample(idx,100)
            push!(ctrl_idx,
                  findall(x -> x in gene_names_ordered[idx],gene_names))
        end
        d = dat[unique(vcat(ctrl_idx...)),:]
        avg_ctrl_cell = d' * ones(Int32,d.m) ./ d.m
        scores[!,k] = avg_feature_cell - avg_ctrl_cell
    end

    if "s" in keys(features) && "g2m" in keys(features)
        assignment = [ Judge(scores[row,:]) for row in 1:size(scores,1) ]
        assignment[assignment .== "unknown"] .= "g1"
    else
        assignment = [ Judge(scores[row,:]) for row in 1:size(scores,1) ]
    end

    # Output
    if isnothing(obj.meta)
        obj.meta = Dict(meta_name => scores)
    else
        obj.meta[meta_name] = scores
    end
    if "s" in keys(features) && "g2m" in keys(features)
        obj.meta[meta_name * "_phase"] = assignment
    else
        obj.meta[meta_name * "_assignment"] = assignment
    end

    # log
    push!(obj.log,String("FeatureScore!(WsObj;features=$features," * 
                         "meta_name=$meta_name,ptr=$ptr)"))

    return "Finished!"
end

# Judge score status
function Judge(x::DataFrameRow)
    loc = findmax(x)
    if loc[1] > 0
        return string(loc[2])
    else
        return "unknown"
    end
end

# Find cutoff for HVGs selection
function FindCutoff1(variances_std::Vector{Float32})
    # 1000 windows and window peak
    bin = range(minimum(variances_std),maximum(variances_std);length=1001)
    dense = Int[]
    for i in 2:1001
        push!(dense,count(bin[i - 1] .< variances_std .< bin[i]))
    end
    pk = findlast(x -> x == maximum(dense),dense)
    # LOESS regression
    loess_model = loess(pk:1000,dense[pk:1000];span=0.05)
    smooth_dense = Loess.predict(loess_model,Float64.(pk:1000))
    return bin[findmin(abs.(smooth_dense .- smooth_dense[1] / 2.0))[2] + pk - 1]
end

# Find cutoff for HVGs selection by a simple, might be unreliable method
function FindCutoff2(avg::Vector{Float64},
        variances_std::Vector{Float32};
        method::Symbol = :simple)
    if method == :dbscan
        # After DBSCAN, auto-selection by top 1/30 avg
        # Good result, but too many HVGs to run fast
        clusters = dbscan(hcat(avg,variances_std)',5e-2)
        big_one = findmax(clusters.counts)[2]
        avg_big = avg[clusters.assignments .== big_one]
        N = round(Int,length(avg) / 30)
        top_avg = sortperm(avg_big;rev=true) |> invperm |> x -> x .<= N
        variances_std_big = variances_std[clusters.assignments .== big_one]
        return mean(variances_std_big[top_avg])
    elseif method == :simple
        # Auto-selection by top 1/30 avg, faster without DBSCAN algorithm
        N = round(Int,length(avg) / 30)
        top_avg = findall(x -> x <= N,sortperm(avg;rev=true) |> invperm)
        return mean(variances_std[top_avg])
    end
    @warn "The method mush be :simple or :quality!"
    return Inf
end
