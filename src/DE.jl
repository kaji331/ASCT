function DE!(obj::WsObj;
        meta::AbstractString = "clusters_latest",
        latent_name::Union{Nothing,AbstractString} = nothing,
        method::AbstractString = "cosg",
        trt::Union{Nothing,Integer,AbstractString,AbstractVector} = nothing,
        ctl::Union{Nothing,Integer,AbstractString,AbstractVector} = nothing,
        seed::Integer = -1,
        min_frac::AbstractFloat = 0.1,
        only_positive::Bool = true,
        fc_thres::AbstractFloat = 0.25,
        features::Union{Nothing,Vector{<: AbstractString}} = nothing,
        adjust_method::AbstractString = "Bonferroni",
        μ::Integer = 1,
        cosg_top::Integer = 50,
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
        overlap = intersect(keys(obj.meta),names(obj.obs))
        if length(overlap) != 0
            throw(DomainError(latent_name,"Overlapping between 'meta' and " * 
                              "'obs'. Bad WsObj!"))
        end
        if latent_name in names(obj.obs)
            latent = obj.obs[!,latent_name]
        elseif latent_name in keys(obj.meta)
            latent = obj.meta[latent_name]
        else
            throw(DomainError(latent_name,"No latent in 'meta' or 'obs'!"))
        end
    else
        latent = nothing
    end

    # Group
    test_groups = obj.meta[meta]
    if !(typeof(test_groups) <: AbstractVector)
        return "Nothing to do! " * 
            "The type of 'meta' content must belong to AbstractVector"
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
    # Convert trt/ctl and the elements of meta to String!
    trt = isnothing(trt) ? trt : string.(trt)
    ctl = isnothing(ctl) ? ctl : string.(ctl)
    test_groups = string.(obj.meta[meta])
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
    obj.meta[meta * "_DE"] = res

    # log
    push!(obj.log,String("DE!(WsObj;" * 
                         "meta = $meta," * 
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
                         "ptr = $ptr"))

    return "Finished!"
end

@inline function DoTest(method::AbstractString,
        trt_grp::Union{Nothing,Integer,AbstractString,AbstractVector},
        ctl_grp::Union{Nothing,Integer,AbstractString,AbstractVector},
        test_groups::Vector{String},
        latent::Union{Nothing,AbstractVector},
        idx_features::Union{Colon,Vector{Int}},
        dat::AbstractSparseMatrix,
        genes::Vector{String},
        min_frac::AbstractFloat,
        fc_thres::AbstractFloat,
        adj_method::AbstractString,
        only_positive::Bool)

    res = DataFrame[]

    if isnothing(trt_grp) && isnothing(ctl_grp)
        # Each group compares REST groups
        @inbounds for trt_grp in sort(unique(test_groups))
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
                           "padj"))
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
                       "padj"))
    elseif isnothing(trt_grp) && !isnothing(ctl_grp)
        # Each group compares one control group
        grp = unique(test_groups) |> sort
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
                           "padj"))
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
                       "padj"))

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

    clusters_label = sort(unique(test_groups))
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
    res = sort(vcat(res...),["group","score"];rev=[false,true])
    # Filter low cell percentage
    filter!(x -> x.score != -1,res)

    # Calculate and filter log2fc, bottleneck
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

    # Top
    res = groupby(res,:group) |> x -> DataFrames.combine(x) do sdf
        first(sdf,cosg_top)
    end

    return res
end
