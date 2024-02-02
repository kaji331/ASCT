"""
    DE!(obj)

Find markers of cell groups.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `group_name::AbstractString = "clusters_latest"`: a obs column to distinguish 
  cell groups.
- `latent_name::Union{Nothing,AbstractString} = nothing`: a obs column 
  represents latent variables.
- `method::AbstractString = "cosg"`: hypothesis testing algorithm.
- `trt::Union{Nothing,AbstractString,Vector{<: AbstractString}} = nothing`: 
  label of treatment groups.
- `ctl::Union{Nothing,AbstractString,Vector{<: AbstractString}} = nothing`: 
  label of control groups.
- `seed::Integer = -1`: random seed.
- `min_frac::AbstractFloat = 0.1`: only test genes are detected in more than 
  this percentage of cells.
- `only_positive::Bool = true`: only return positive markers.
- `fc_thres::AbstractFloat = 0.25`: return genes with higher log2-fold change.
- `features::Union{Nothing,Vector{<: AbstractString}} = nothing`: feature name 
  vector for test.
- `adjust_method::AbstractString = "Bonferroni"`: algorithm for adjusted 
  p-value.
- `μ::Integer = 1`: cosg penalty factor to gene expression in cells not 
  belonging to the group of interest.
- `cosg_top::Integer = 50`: keep top results of cosg algorithm.
- `sub_group::Tuple{Vector{<: AbstractString},Vector{<: AbstractString}} = 
  (String[],[])`: set sub-groups for test.
- `ptr::Union{Symbol,AbstractString} = :auto`: use "norm\\_dat" or 
  "regress\\_dat".

# Example
```julia-repl
julia> DE!(obj;method="lr",latent_name="WsObj_batch")
julia> DE!(obj;group_name="clusters_MC_0.4",sub_group=(["WsObj_batch"],["1"]))
```
"""
function DE!(obj::WsObj;
        group_name::AbstractString = "clusters_latest",
        latent_name::Union{Nothing,AbstractString} = nothing,
        method::AbstractString = "cosg",
        trt::Union{Nothing,AbstractString,Vector{<: AbstractString}} = nothing,
        ctl::Union{Nothing,AbstractString,Vector{<: AbstractString}} = nothing,
        seed::Integer = -1,
        min_frac::AbstractFloat = 0.1,
        only_positive::Bool = true,
        fc_thres::AbstractFloat = 0.25,
        features::Union{Nothing,Vector{<: AbstractString}} = nothing,
        adjust_method::AbstractString = "Bonferroni",
        μ::Integer = 1,
        cosg_top::Integer = 50,
        sub_group::Tuple{Vector{<: AbstractString},Vector{<: AbstractString}} = 
            (String[],String[]),
        ptr::Union{Symbol,AbstractString} = :auto)

    Random.seed!(seed < 0 ? 1984 : seed)
    # Check parameters
    if ptr == :auto && "regress_dat" in keys(obj.dat)
        ptr = "regress"
    elseif ptr == :auto && !("regress_dat" in keys(obj.dat))
        ptr = "norm"
    end
    if !(0 <= min_frac <= 1)
        return "Nothing to do! 'min_frac' is only support between 0 and 1!"
    end
    if fc_thres < 0
        return "Nothing to do! 'fc_thres' must be >= 0!"
    end
    if μ < 1
        return "Nothing to do! 'μ' must be >= 1!"
    end
    if cosg_top < 1
        return "Nothing to do! 'cosg_top' must be >= 1!"
    end

    dat = obj.dat[ptr * "_dat"]

    # Latent data
    if !isnothing(latent_name)
        if method != "lr"
            @warn "Only 'lr' method supports 'latent_name'!" * 
                "Other methods will ignore this variable!"
        end
        if latent_name in names(obj.obs)
            latent = obj.obs[!,latent_name]
        else
            return "Nothing to do! No latent in the 'obs'!"
        end
    else
        latent = nothing
    end

    # Features
    if !isnothing(features)
        idx = findall(x -> x in obj.var.name,features)
        dat = dat[idx,:]
        genes = features[idx]
    else
        idx = Colon()
        genes = obj.var.name
    end

    # Method
    if !(method in ["wilcox","ttest","lr","cosg"])
        return "Nothing to do! The method only support " * "
            'wilcox','ttest', 'lr' and 'cosg'!"
    end

    # Calculation
    test_groups = string.(obj.obs[!,group_name])

    # If you wanna compare clusters only in library 1 and Female, you can make 
    # 'sub_group' be (["library","gender"],["1","Female"])...
    if !isempty(sub_group[1])
        pre_sub_idx = nothing
        curr_sub_idx = nothing
        for i in eachindex(sub_group[1])
            if i == 1
                curr_sub_idx = 
                    string.(obj.obs[!,sub_group[1][i]]) .== 
                    string(sub_group[2][i])
                pre_sub_idx = curr_sub_idx
            else
                curr_sub_idx = 
                    string.(obj.obj[!,sub_group[1][i]]) .== 
                    string(sub_group[2][i]) |> x -> x .& pre_sub_idx
                if i != length(sub_group[1])
                    pre_sub_idx = curr_sub_idx
                end
            end
        end
        test_groups = test_groups[curr_sub_idx]
        dat = dat[:,curr_sub_idx]
    end

    if method != "cosg"
        res = DoTest(method,
                     trt,
                     ctl,
                     test_groups,
                     latent,
                     idx,
                     dat,
                     genes,
                     min_frac,
                     fc_thres,
                     adjust_method,
                     only_positive)
    else
        if !isnothing(trt) || !isnothing(ctl)
            @warn "The 'cosg' method does not support defined 'trt'/'ctl'!"
            return "Nothing to do!"
        end
        res = COSG(test_groups,
                   idx,
                   dat,
                   genes,
                   min_frac,
                   fc_thres,
                   μ,
                   cosg_top,
                   only_positive)
    end

    # Output
    obj.meta[group_name * "_DE"] = res

    # log
    push!(obj.log,String("DE!(WsObj;" * 
                         "group_name = $group_name," * 
                         "latent_name = $latent_name," * 
                         "method = $method," *
                         "trt = $trt," * 
                         "ctl = $ctl," * 
                         "seed = $seed," * 
                         "min_frac = $min_frac," * 
                         "only_positive = $only_positive," * 
                         "fc_thres = $fc_thres," * 
                         "features = $features," * 
                         "adjust_method = $adjust_method," * 
                         "μ = $μ," * 
                         "cosg_top = $cosg_top," * 
                         "ptr = $ptr)"))

    return "Finished!"
end

@inline function DoTest(method::AbstractString,
        trt_grp::Union{Nothing,AbstractString,Vector{<: AbstractString}},
        ctl_grp::Union{Nothing,AbstractString,Vector{<: AbstractString}},
        test_groups::Vector{String},
        latent::Union{Nothing,AbstractVector},
        idx_features::Union{Colon,Vector{<: Integer}},
        dat::AbstractSparseMatrix,
        genes::Vector{<: AbstractString},
        min_frac::AbstractFloat,
        fc_thres::AbstractFloat,
        adj_method::AbstractString,
        only_positive::Bool)

    res = DataFrame[]

    if isnothing(trt_grp) && isnothing(ctl_grp)
        # Each group compares REST groups
        @inbounds for trt_grp in sort(unique(test_groups);lt=natural)
            idx_trt = findall(x -> x == trt_grp,test_groups)
            idx_ctl = findall(x -> x != trt_grp,test_groups)
            cell_number_trt = length(idx_trt)
            cell_number_ctl = length(idx_ctl)
            dat_trt = Matrix(dat[idx_features,idx_trt])
            dat_ctl = Matrix(dat[idx_features,idx_ctl])
            latent = isnothing(latent) ? nothing : latent[vcat(idx_trt,idx_ctl)]
            rn = size(dat_trt,1)
            # Filtering genes by cell percentages like Seurat
            pval = Float64[]
            idx_feature_pct = Int[]
            fc_pct = Float64[]
            pct1_pct = Float64[]
            pct2_pct = Float64[]
            @simd for idx_feature in 1:rn
                trt = @view dat_trt[idx_feature,:]
                ctl = @view dat_ctl[idx_feature,:]
                pct1 = round(count(x -> x != 0,trt) / cell_number_trt;
                             digits=3)
                pct2 = round(count(x -> x != 0,ctl) / cell_number_ctl;
                             digits=3)
                if (pct1 > min_frac || pct2 > min_frac) && pct1 != pct2
                    # Calculate and filter fold change
                    # Default method from Seurat v4
                    fc = log2(mean(expm1.(trt)) + 1) - 
                        log2(mean(expm1.(ctl)) + 1)
                    if only_positive && fc >= fc_thres
                        push!(idx_feature_pct,idx_feature)
                        push!(fc_pct,fc)
                        push!(pct1_pct,pct1)
                        push!(pct2_pct,pct2)
                        push!(pval,PvalFun(method,trt,ctl;latent))
                    elseif !only_positive && abs(fc) >= fc_thres
                        push!(idx_feature_pct,idx_feature)
                        push!(fc_pct,fc)
                        push!(pct1_pct,pct1)
                        push!(pct2_pct,pct2)
                        push!(pval,PvalFun(method,trt,ctl;latent))
                    end
                end
            end
            padj = MultipleTesting.adjust(pval,getfield(MultipleTesting,
                                                        Symbol(adj_method))())
            push!(res,sort(DataFrame(gene=genes[idx_feature_pct],
                                     log2fc=fc_pct,
                                     pval=pval,
                                     padj=padj,
                                     pct1=pct1_pct,
                                     pct2=pct2_pct,
                                     group=trt_grp * " .vs. REST"),
                           "log2fc";rev=true))
        end
    elseif !isnothing(trt_grp) && isnothing(ctl_grp)
        # One group compares a REST group
        idx_trt = findall(x -> any(x .== trt_grp),test_groups)
        idx_ctl = findall(x -> all(x .!= trt_grp),test_groups)
        cell_number_trt = length(idx_trt)
        cell_number_ctl = length(idx_ctl)
        dat_trt = Matrix(dat[idx_features,idx_trt])
        dat_ctl = Matrix(dat[idx_features,idx_ctl])
        latent = isnothing(latent) ? nothing : latent[vcat(idx_trt,idx_ctl)]
        rn = size(dat_trt,1)
        # Filtering genes by cell percentages like Seurat
        pval = Float64[]
        idx_feature_pct = Int[]
        fc_pct = Float64[]
        pct1_pct = Float64[]
        pct2_pct = Float64[]
        @simd for idx_feature in 1:rn
            trt = @view dat_trt[idx_feature,:]
            ctl = @view dat_ctl[idx_feature,:]
            pct1 = round(count(x -> x != 0,trt) / cell_number_trt;
                         digits=3)
            pct2 = round(count(x -> x != 0,ctl) / cell_number_ctl;
                         digits=3)
            if (pct1 > min_frac || pct2 > min_frac) && pct1 != pct2
                # Calculate and filter fold change
                fc = log2(mean(expm1.(trt)) + 1) - 
                    log2(mean(expm1.(ctl)) + 1)
                if only_positive && fc >= fc_thres
                    push!(idx_feature_pct,idx_feature)
                    push!(fc_pct,fc)
                    push!(pct1_pct,pct1)
                    push!(pct2_pct,pct2)
                    push!(pval,PvalFun(method,trt,ctl;latent))
                elseif !only_positive && abs(fc) >= fc_thres
                    push!(idx_feature_pct,idx_feature)
                    push!(fc_pct,fc)
                    push!(pct1_pct,pct1)
                    push!(pct2_pct,pct2)
                    push!(pval,PvalFun(method,trt,ctl;latent))
                end
            end
        end
        padj = MultipleTesting.adjust(pval,getfield(MultipleTesting,
                                                    Symbol(adj_method))())
        push!(res,sort(DataFrame(gene=genes[idx_feature_pct],
                                 log2fc=fc_pct,
                                 pval=pval,
                                 padj=padj,
                                 pct1=pct1_pct,
                                 pct2=pct2_pct,
                                 group=join(trt_grp,"|") * " .vs. REST"),
                       "log2fc";rev=true))
    elseif isnothing(trt_grp) && !isnothing(ctl_grp)
        # Each group compares one control group
        grp = unique(test_groups) |> x -> sort(x;lt=natural)
        grp = grp[map(x -> all(x .!= ctl_grp),grp)]
        @inbounds for trt_grp in grp
            idx_trt = findall(x -> x == trt_grp,test_groups)
            idx_ctl = findall(x -> any(x .== ctl_grp),test_groups)
            cell_number_trt = length(idx_trt)
            cell_number_ctl = length(idx_ctl)
            dat_trt = Matrix(dat[idx_features,idx_trt])
            dat_ctl = Matrix(dat[idx_features,idx_ctl])
            latent = isnothing(latent) ? nothing : latent[vcat(idx_trt,idx_ctl)]
            rn = size(dat_trt,1)
            # Filtering genes by cell percentages like Seurat
            pval = Float64[]
            idx_feature_pct = Int[]
            fc_pct = Float64[]
            pct1_pct = Float64[]
            pct2_pct = Float64[]
            @simd for idx_feature in 1:rn
                trt = @view dat_trt[idx_feature,:]
                ctl = @view dat_ctl[idx_feature,:]
                pct1 = round(count(x -> x != 0,trt) / cell_number_trt;
                             digits=3)
                pct2 = round(count(x -> x != 0,ctl) / cell_number_ctl;
                             digits=3)
                if (pct1 > min_frac || pct2 > min_frac) && pct1 != pct2
                    # Calculate and filter fold change
                    fc = log2(mean(expm1.(trt)) + 1) - 
                        log2(mean(expm1.(ctl)) + 1)
                    if only_positive && fc >= fc_thres
                        push!(idx_feature_pct,idx_feature)
                        push!(fc_pct,fc)
                        push!(pct1_pct,pct1)
                        push!(pct2_pct,pct2)
                        push!(pval,PvalFun(method,trt,ctl;latent))
                    elseif !only_positive && abs(fc) >= fc_thres
                        push!(idx_feature_pct,idx_feature)
                        push!(fc_pct,fc)
                        push!(pct1_pct,pct1)
                        push!(pct2_pct,pct2)
                        push!(pval,PvalFun(method,trt,ctl;latent))
                    end
                end
            end
            padj = MultipleTesting.adjust(pval,getfield(MultipleTesting,
                                                        Symbol(adj_method))())
            push!(res,sort(DataFrame(gene=genes[idx_feature_pct],
                                     log2fc=fc_pct,
                                     pval=pval,
                                     padj=padj,
                                     pct1=pct1_pct,
                                     pct2=pct2_pct,
                                     group=trt_grp * " .vs. " * 
                                        join(ctl_grp,"|")),
                           "log2fc";rev=true))
        end
    else
        # One group compares one control group
        idx_trt = findall(x -> any(x .== trt_grp),test_groups)
        idx_ctl = findall(x -> any(x .== ctl_grp),test_groups)
        cell_number_trt = length(idx_trt)
        cell_number_ctl = length(idx_ctl)
        dat_trt = Matrix(dat[idx_features,idx_trt])
        dat_ctl = Matrix(dat[idx_features,idx_ctl])
        latent = isnothing(latent) ? nothing : latent[vcat(idx_trt,idx_ctl)]
        rn = size(dat_trt,1)
        # Filtering genes by cell percentages like Seurat
        pval = Float64[]
        idx_feature_pct = Int[]
        fc_pct = Float64[]
        pct1_pct = Float64[]
        pct2_pct = Float64[]
        for idx_feature in 1:rn
            trt = @view dat_trt[idx_feature,:]
            ctl = @view dat_ctl[idx_feature,:]
            pct1 = round(count(x -> x != 0,trt) / cell_number_trt;
                         digits=3)
            pct2 = round(count(x -> x != 0,ctl) / cell_number_ctl;
                         digits=3)
            if (pct1 > min_frac || pct2 > min_frac) && pct1 != pct2
                # Calculate and filter fold change
                fc = log2(mean(expm1.(trt)) + 1) - 
                    log2(mean(expm1.(ctl)) + 1)
                if only_positive && fc >= fc_thres
                    push!(idx_feature_pct,idx_feature)
                    push!(fc_pct,fc)
                    push!(pct1_pct,pct1)
                    push!(pct2_pct,pct2)
                    push!(pval,PvalFun(method,trt,ctl;latent))
                elseif !only_positive && abs(fc) >= fc_thres
                    push!(idx_feature_pct,idx_feature)
                    push!(fc_pct,fc)
                    push!(pct1_pct,pct1)
                    push!(pct2_pct,pct2)
                    push!(pval,PvalFun(method,trt,ctl;latent))
                end
            end
        end
        padj = MultipleTesting.adjust(pval,getfield(MultipleTesting,
                                                    Symbol(adj_method))())
        push!(res,sort(DataFrame(gene=genes[idx_feature_pct],
                                 log2fc=fc_pct,
                                 pval=pval,
                                 padj=padj,
                                 pct1=pct1_pct,
                                 pct2=pct2_pct,
                                 group=join(trt_grp,"|") * " .vs. " * 
                                    join(ctl_grp,"|")),
                       "log2fc";rev=true))

    end

    return vcat(res...)
end

function PvalFun(method::AbstractString,
        trt::AbstractVector,
        ctl::AbstractVector;
        latent::Union{Nothing,AbstractVector} = nothing)::Float64

    if method == "wilcox"
        return MannWhitneyUTest(trt,ctl) |> pvalue
    elseif method == "ttest"
        return UnequalVarianceTTest(trt,ctl) |> pvalue
    elseif method == "lr"
        if isnothing(latent)
            d = DataFrame(gene=vcat(trt,ctl),
                          group=vcat(repeat([0],length(trt)),
                                     repeat([1],length(ctl))))
            nullmodel = glm(@formula(group ~ 1),d,Binomial(),LogitLink())
            model = glm(@formula(group ~ gene),d,Binomial(),LogitLink())
            return lrtest(model,nullmodel).pval[end]
        else
            d = DataFrame(gene=vcat(trt,ctl),
                          group=vcat(repeat([0],length(trt)),
                                     repeat([1],length(ctl))),
                          latent=latent)
            nullmodel = glm(@formula(group ~ latent),d,Binomial(),
                            LogitLink())
            model = glm(@formula(group ~ gene + latent),d,Binomial(),
                        LogitLink())
            return lrtest(model,nullmodel).pval[end]
        end
    end
end

@inline function COSG(test_groups::Vector{String},
        idx_features::Union{Colon,Vector{Int}},
        dat::AbstractSparseMatrix,
        genes::Vector{String},
        min_frac::AbstractFloat,
        fc_thres::AbstractFloat,
        μ::Integer,
        cosg_top::Integer,
        only_positive::Bool)

    res = DataFrame[]

    clusters_label = sort(unique(test_groups);lt=natural)
    clusters_number = length(clusters_label)
    clusters_matrix = zeros(clusters_number,dat.n)
    for cluster_idx in eachindex(clusters_label)
        clusters_matrix[cluster_idx,
                        test_groups .== clusters_label[cluster_idx]] .= 1
    end
    similarity = 1.0 .- pairwise(CosineDist(),Matrix(dat)',clusters_matrix')
    λ = similarity .^ 2
    ϵ = λ * ones(Int,clusters_number)
    λ = λ ./ ((1 - μ) .* λ .+ μ .* repeat([ϵ;],outer=(1,clusters_number)))
    λ .*= similarity
    @inbounds for cluster_idx in eachindex(clusters_label)
        dat_cluster = dat[:,test_groups .== clusters_label[cluster_idx]]' |> 
            sparse
        # percentage
        pct = [ dat_cluster.colptr[col + 1] - dat_cluster.colptr[col] 
               for col in 1:dat_cluster.n ] ./ dat_cluster.m
        λ[pct .< min_frac,cluster_idx] .= -1
        df = DataFrame(gene=genes,
                       score=λ[:,cluster_idx],
                       pct1=pct,
                       group=clusters_label[cluster_idx])
        push!(res,df)
    end

    # Sort
    res = sort(vcat(res...),["group","score"];rev=[false,true],lt=natural)
    # Filter low cell percentage
    filter!(x -> x.score != -1,res)

    # Calculate and filter log2fc, bottleneck (now @btime about 2s, @btime 
    # about 1.5s without log2fc calculation)
    all_fc = Float64[]
    @inbounds for cluster in clusters_label
        row_idx = indexin(res[res.group .== cluster,:gene],genes) |> 
            x -> x[BitVector(1 .- isnothing.(x))] |> x -> Int.(x)
        dat_cluster = dat[row_idx,test_groups .== cluster]' |> sparse
        dat_not_cluster = dat[row_idx,test_groups .!= cluster]' |> sparse
        fc = [ log2(mean(expm1.(Vector(dat_cluster[:,col]))) + 1) - 
              log2(mean(expm1.(Vector(dat_not_cluster[:,col]))) + 1) 
              for col in 1:dat_cluster.n ]
        append!(all_fc,fc)
    end
    res[!,:log2fc] = all_fc
    if only_positive
        filter!(x -> x.log2fc >= fc_thres,res)
    else
        filter!(x -> abs(x.log2fc) >= fc_thres,res)
    end
    # Resort
    res = sort(res,["group","log2fc"];rev=[false,true],lt=natural)

    # Top
    res = groupby(res,:group) |> x -> DataFrames.combine(x) do sdf
        first(sdf,cosg_top)
    end

    return res
end
