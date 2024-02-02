mutable struct HarmonyObj
    ρ::AbstractMatrix
    z_orig::AbstractMatrix
    z_corr::AbstractMatrix
    z_cos::AbstractMatrix
    γ::AbstractMatrix
    ϕ::AbstractMatrix
    ϕ_moe::AbstractMatrix
    batches_ratio::AbstractVector
    θ::AbstractVector
    batches_number::AbstractVector
    σ::AbstractVector
    λ::Diagonal
    obj_harmony::AbstractVector
    obj_kmeans::AbstractVector
    obj_kmeans_dist::AbstractVector
    obj_kmeans_entropy::AbstractVector
    obj_kmeans_cross::AbstractVector
    kmeans_rounds::AbstractVector
    block_size::AbstractFloat
    ϵ_kmeans::AbstractFloat
    ϵ_harmony::AbstractFloat
    nclust::Integer
    β::Integer
    cells_number::Integer
    max_iter_kmeans::Integer
    win_size::Integer
    scale_dist::AbstractMatrix
    dist_mat::AbstractMatrix
    Ω::AbstractMatrix
    Ε::AbstractMatrix
    ϕ_rk::AbstractMatrix
    ω::AbstractMatrix
    ran_setup::Bool
    ran_init::Bool
end

"""
    Harmony!(objs)

Integrate single-cell datasets using harmony algorithm.

# Arguments
- `objs::Union{WsObj,Vector{WsObj}}`: a single-cell WsObj struct or WsObj vector.

# Keyword Arguments
- `project_name::AbstractString = "WsObj_batch"`: a column name in obs for 
  integrated data.
- `batches::Union{AbstractString,Vector{<: AbstractString},Symbol} = :auto`: 
  :auto will make each WsObj as a batch, one string would find the column name 
  in obs for just one WsObj and the vector will represent the batch of each 
  WsObj or cell based the length of vector.
- `ref_batches::Union{Vector{<: AbstractString},Nothing} = nothing`: the batch 
  names for optional reference batches.
- `hvg_number::Union{Symbol,Integer} = 2000`: the top number of features for 
  highly variable selection or the symbol :auto for auto-determine of numbers.
- `hvg_method::Symbol = :simple`: the symbol of :simple, :loess and :dbscan for 
  different method used by highly variable number selection.
- `max_pca::Integer = 50`: max principal components would be calculated.
- `pca_method::AbstractString = "PCA"`: "Arpack" or "TSVD" might be better for 
  large data.
- `use_pc::Union{Integer,Symbol} = :auto`: only 1:cut PCs will be used 
  downstream automatically.
- `min_features::Union{Nothing,Integer} = nothing`: drop cells containing 
  features less than this number.
- `min_cells::Union{Nothing,Integer} = nothing`: drop features containing 
  cells less than this number.
"""
function Harmony!(objs::Union{WsObj,Vector{WsObj}};
        project_name::AbstractString = "WsObj_batch",
        batches::Union{AbstractString,Vector{<: AbstractString},Symbol} = :auto,
        ref_batches::Union{Vector{<: AbstractString},Nothing} = nothing,
        hvg_number::Union{Symbol,Integer} = 2000,
        hvg_method::Symbol = :simple,
        max_pca::Integer = 50,
        pca_method::AbstractString = "PCA",
        use_pc::Union{Symbol,Integer} = :auto,
        min_features::Union{Nothing,Integer} = nothing,
        min_cells::Union{Nothing,Integer} = nothing)

    # Check batches variable
    if batches != :auto && 
        length(batches) != 1 && 
        (typeof(objs) <: AbstractVector ? (length(batches) != length(objs) && 
                                           length(batches) != 
                                            sum([ obj.dat["raw_dat"].n 
                                                 for obj in objs ])) : 
         length(batches) != objs.dat["raw_dat"].n)
        return "Nothing to do! One 'batches' for a meta name in WsObj(s)," * 
            "or the length of batches should be equal to the length of WsObj " * 
            "vector or the cell numbers of WsObj(s)!"
    end

    # Groups
    mark = 0
    if batches == :auto
        if typeof(objs) <: AbstractVector
            grp = string.(vcat([ repeat([i];inner=size(objs[i].obs,1)) 
                                for i in eachindex(objs) ]...))
        else
            @warn "One WsObj should explicitly confirm the 'batches'! " * 
                "Now try to use the 'project_name'..."
            grp = objs.obs[!,project_name]
        end
    elseif length(batches) == 1
        # Example: ["batch_labels"]
        if typeof(objs) <: AbstractVector
            grp = vcat([ obj.obs[!,batches] for obj in objs ]...)
        else
            grp = objs.obs[!,batches]
        end
    else
        # Example: ["trt","trt","ctl","ctl"]
        if typeof(objs) <: AbstractVector && length(batches) == length(objs)
            grp = vcat([ batches[i] in names(objs[i].obs) ? 
                        objs[i].obs[!,batches[i]] : 
                        repeat([batches[i]],
                               objs[i].dat["raw_dat"].n) 
                        for i in eachindex(objs) ]...)
        elseif typeof(objs) <: AbstractVector && 
            length(batches) == sum([ obj.dat["raw_dat"].n for obj in objs ])
            grp = batches
        elseif typeof(objs) <: WsObj && length(batches) == objs.dat["raw_dat"].n
            grp = batches
            mark = 1
        else
            return "Nothing to do! Length of 'batches' is not corrected!"
        end
    end

    if length(unique(grp)) <= 1
        return "Nothing to do! Only one batch???"
    end

    # Check ref_batches variable
    if !isnothing(ref_batches)
        idx = indexin(ref_batches,unique(grp))
        if !(findfirst(x -> isnothing(x),idx) |> x -> isnothing(x))
            return "Nothing to do! Reference batches were not existed!"
        end
    end

    # Check for merging data
    if typeof(objs) <: AbstractVector
        obj = MergeRawData(objs,
                           Dict{AbstractString,
                                AbstractVector}(project_name => grp);
                           min_features=min_features,min_cells=min_cells)
        grp = DataFrame(project_name => obj.obs[!,project_name])
    elseif !(isnothing(min_features) && isnothing(min_cells))
        obj = objs
        if mark == 1
            obj.obs[!,project_name] = grp
            SimpleFilter!(obj,min_features,min_cells)
            grp = DataFrame(project_name => obj.obs[!,project_name])
        elseif batches == :auto
            SimpleFilter!(obj,min_features,min_cells)
            grp = DataFrame(project_name => obj.obs[!,project_name])
        else
            SimpleFilter!(obj,min_features,min_cells)
            grp = DataFrame(project_name => obj.obs[!,batches])
        end
        mark = 1
    else
        obj = objs
        if mark == 1
            obj.obs[!,project_name] = grp
        end
        grp = DataFrame(project_name => grp)
        mark = 1
    end
    objs = nothing
    GC.gc()

    # Normalization
    @info "Normalizing..."
    NormalizeData!(obj)
    # HVGs
    @info "Looking for HVGs..."
    SelectHVG!(obj;hvg_number=hvg_number,method=hvg_method)
    # PCA, choose PCs automatically
    @info "Doing PCA..."
    PCA!(obj;max_pc=max_pca,method=pca_method,cut=use_pc)

    # Harmony matrix
    @info "Doing Harmony..."
    hm = HarmonyMatrix(obj.meta["pca"][:,1:obj.meta["pca_cut"]]',grp,names(grp);
                       ref_batches=ref_batches)'

    # Modify WsObj
    if !isnothing(hm)
        obj.meta["harmony"] = hm
    else
        @warn "No harmony result: clustering error!"
    end

    if mark == 0
        return obj
    end
end

function HarmonyMatrix(em::AbstractMatrix,
        grp::DataFrame,
        grp_names::Vector{String};
        σ::AbstractFloat = 0.1,
        block_size::AbstractFloat = 0.05,
        max_iter_harmony::Integer = 10,
        max_iter_cluster::Integer = 20,
        ϵ_cluster::AbstractFloat = 1e-5,
        ϵ_harmony::AbstractFloat = 1e-4,
        ref_batches::Union{Nothing,String,Vector{String}} = nothing)

    # Now the row is PCs, and the column is Cells...
    pc_num,cells_number = size(em)
    groups = grp_names
    nclust = minimum([round(Int,cells_number / 30),100])
    θ = repeat([2];inner=length(groups))
    λ = repeat([1];inner=length(groups))
    if nclust > 1
        σ = repeat([σ];inner=nclust)
    end

    # One-hot
    # https://discourse.julialang.org/t/all-the-ways-to-do-one-hot-encoding/64807
    ϕ = vcat([ unique(grp[:,g]) .== permutedims(grp[:,g]) for g in groups ]...)
    # model = ConstantTerm(1) ~ sum(term.(groups))
    # ϕ = modelmatrix(model,grp) |> x -> hcat(Int.(x) .⊻ 1,Int.(x))'

    # Max of Int32 is 2147483647.
    batches_number = ϕ * ones(Int32,size(ϕ,2))
    batches_ratio = batches_number ./ cells_number
    batches_class = [ length(unique(col)) for col in eachcol(grp) ]
    θ = vcat([ repeat([θ[i]];inner=batches_class[i]) 
              for i in eachindex(batches_class) ]...)
    τ = 0
    θ = θ .* (1 .- exp.(-(batches_number ./ (nclust * τ)) .^ 2)) # Still θ ?
    λ = vcat([ repeat([λ[i]];inner=batches_class[i]) 
              for i in eachindex(batches_class) ]...)
    λ_mat = Diagonal(vcat(0,λ))

    # Reference
    if isnothing(ref_batches)
        ϕ_moe = vcat(repeat([1],cells_number)',ϕ)
    else
        idx = indexin(ref_batches,grp.groups)
        ref_i = ((ϕ[idx,:] .== 1)' * ones(Int32,length(idx))) .>= 1
        ϕ_moe = ϕ'[:,Not(idx)]' |> Matrix
        # ϕ_moe[:,ref_i] .= 0
        ϕ_moe = vcat(repeat([1],cells_number)',ϕ_moe)
        λ_idx = vcat(1,eachindex(grp.groups)[Not(idx)] .+ 1)
        λ_mat = λ_mat[λ_idx,λ_idx]
    end

    # Run
    # In Armadillo, % means .*, as_scalar just makes 1 element vector to scalar, 
    # mean/sum have dimension option, accu equals to sum, repmat means 
    # repeat with outer option
    z_cos = copy(em)
    CosNormal!(z_cos,0,true)
    β = size(ϕ,1)
    Random.seed!(1984)
    γ = kmeans(z_cos,nclust;maxiter=25).centers
    CosNormal!(γ,0,false)
    dist_mat = 2 .* (1 .- γ' * z_cos)
    ρ = -dist_mat ./ σ
    cm = [ maximum(col) for col in eachcol(ρ) ]
    @inbounds @simd for i in 1:size(ρ,1)
        ρ'[:,i] .-= cm
    end
    ρ = exp.(ρ)
    cs = [ sum(col) for col in eachcol(ρ) ]
    @inbounds @simd for i in 1:size(ρ,1)
        ρ'[:,i] ./= cs
    end
    rs = [ sum(row) for row in eachrow(ρ) ]
    Ε = rs * batches_ratio'
    Ω = ρ * ϕ'

    kmeans_error = sum(ρ .* dist_mat)
    safe_entropy = ρ .* log.(ρ)
    safe_entropy[isnan.(safe_entropy)] .= 0
    safe_entropy[isinf.(safe_entropy)] .= 0
    entropy = sum(safe_entropy .* σ)
    cross_entropy = sum((ρ .* σ) .* 
                        ((repeat(θ';outer=(nclust,1)) .* 
                          log.(Float32.((Ω .+ 1) ./ (Ε .+ 1)))) * 
                         ϕ))
    obj_kmeans = [kmeans_error + entropy + cross_entropy]
    obj_kmeans_dist = [kmeans_error]
    obj_kmeans_entropy = [entropy]
    obj_kmeans_cross = [cross_entropy]
    obj_harmony = [obj_kmeans[end]]

    ho = HarmonyObj(ρ,em,copy(em),z_cos,γ,ϕ,ϕ_moe,batches_ratio,θ,
                    batches_number,σ,λ_mat,obj_harmony,obj_kmeans,
                    obj_kmeans_dist,obj_kmeans_entropy,obj_kmeans_cross,Int[],
                    block_size,ϵ_cluster,ϵ_harmony,nclust,β,cells_number,
                    max_iter_cluster,3,zeros(nclust,cells_number),dist_mat,Ω,Ε,
                    zeros(β + 1,cells_number),zeros(β + 1,cells_number),true,
                    true)

    if Harmonize!(ho,max_iter_harmony) == 2
        return nothing
    end

    return ho.z_corr
end

function CosNormal!(m::AbstractMatrix,
        margin::Integer,
        do_safe::Bool)
    if margin == 1
        @inbounds for row in eachrow(m)
            if do_safe
                # Need a temporary variable for corrected result
                M = maximum(row)
                replace!(x -> x / M,row)
            end
            # Need a temporary variable for corrected result
            N = norm(row,2)
            replace!(x -> x / N,row)
        end
        return "Finished!"
    elseif margin == 0
        @inbounds for col in eachcol(m)
            if do_safe
                # Need a temporary variable for corrected result
                M = maximum(col)
                replace!(x -> x / M,col)
            end
            # Need a temporary variable for corrected result
            N = norm(col,2)
            replace!(x -> x / N,col)
        end
        return "Finished!"
    end

    @warn "margin mush be 0 or 1! Nothing is changed!"
end

function ComputeΓ(z_cos::AbstractMatrix,
        ρ::AbstractMatrix)
    m = z_cos * ρ'
    @inbounds for col in eachcol(m)
        N = norm(col)
        replace!(x -> x / N,col)
    end
    return m
end

function HarmonyPow(m::AbstractMatrix,
        t::AbstractVector)
    n = copy(m)
    @inbounds for i in 1:size(n,2)
        n[:,i] .^= t[i]
    end
    return n
end

function UpdateΡ!(ho::HarmonyObj)
    Random.seed!(1984)
    update_order = shuffle(1:ho.cells_number)
    scale_dist = -ho.dist_mat
    scale_dist ./= ho.σ
    col_max = [ maximum(col) for col in eachcol(scale_dist) ]
    scale_dist = scale_dist'
    @inbounds for i in 1:size(scale_dist,2)
        scale_dist[:,i] .-= col_max
    end
    scale_dist = exp.(scale_dist')

    # Online updates
    n_blocks = Int(ceil(1 / ho.block_size))
    cells_per_block = trunc(Int,ho.cells_number / n_blocks)
    # Gather cell updates indices
    idx_min = [ i * cells_per_block for i in 1:n_blocks ]
    idx_max = [ minimum([tmp + cells_per_block,ho.cells_number]) for tmp in idx_min ]
    if any(idx_min .> idx_max)
        pos = findfirst(idx_min .> idx_max) - 1
    else
        pos = n_blocks
    end
    @inbounds for i in 1:pos
        idx_list = Int.(trunc.(range(idx_min[i],idx_max[i];
                                     length=idx_max[i] - idx_min[i])))
        cells_update = update_order[idx_list]

        # Step1. remove cells
        ho.Ε .-= ho.ρ[:,cells_update] * ones(Int32,length(cells_update)) * 
            ho.batches_ratio'
        ho.Ω .-= ho.ρ[:,cells_update] * ho.ϕ[:,cells_update]'

        # Step2. re-compute ρ for moved cells
        ho.ρ[:,cells_update] = scale_dist[:,cells_update]
        ho.ρ[:,cells_update] .*= HarmonyPow(Float32.((ho.Ε .+ 1) ./ 
                                                     (ho.Ω .+ 1)),ho.θ) * 
            ho.ϕ[:,cells_update]
        @inbounds for col in eachcol(@view ho.ρ[:,cells_update])
            l1 = norm(col,1)
            replace!(x -> x / l1,col)
        end

        # Step3. put cells back
        ho.Ε .+= ho.ρ[:,cells_update] * ones(Int32,length(cells_update)) * 
            ho.batches_ratio'
        ho.Ω .+= ho.ρ[:,cells_update] * ho.ϕ[:,cells_update]'
    end

    return false
end

function CheckConvergence(ho::HarmonyObj,
        t::Integer)

    if !(t in [0,1])
        return true
    end

    if t == 0
        # Clustering
        # Compute new window mean
        oo = 0
        on = 0
        @inbounds for i in 1:ho.win_size
            oo += ho.obj_kmeans[length(ho.obj_kmeans) - i]
            on += ho.obj_kmeans[length(ho.obj_kmeans) - i + 1]
        end
        if ((oo - on) / abs(oo)) < ho.ϵ_kmeans
            return true
        else
            return false
        end
    else
        # Harmony
        oo = ho.obj_harmony[length(ho.obj_harmony) - 1]
        on = ho.obj_harmony[length(ho.obj_harmony)]
        if ((oo - on) / abs(oo)) < ho.ϵ_harmony
            return true
        else
            return false
        end
    end
end

function Cluster!(ho::HarmonyObj)

    if ho.ran_init
        err_status = false
        ho.dist_mat = 2 .* (1 .- ho.γ' * ho.z_cos)

        iter = 1
        @inbounds for i in 1:ho.max_iter_kmeans
            # Update γ
            ho.γ = ComputeΓ(ho.z_cos,ho.ρ)
            ho.dist_mat = 2 .* (1 .- ho.γ' * ho.z_cos)

            # Update ρ
            err_status = UpdateΡ!(ho)
            if err_status
                return err_status
            end

            kmeans_error = sum(ho.ρ .* ho.dist_mat)
            safe_entropy = ho.ρ .* log.(ho.ρ)
            safe_entropy[isnan.(safe_entropy)] .= 0
            safe_entropy[isinf.(safe_entropy)] .= 0
            entropy = sum(safe_entropy .* ho.ρ)
            cross_entropy = sum((ho.ρ .* ho.σ) .* 
                                ((repeat(ho.θ';outer=(ho.nclust,1)) .* 
                                  log.(Float32.((ho.Ω .+ 1) ./ (ho.Ε .+ 1)))) * 
                                 ho.ϕ))
            push!(ho.obj_kmeans,kmeans_error + entropy,cross_entropy)
            push!(ho.obj_kmeans_dist,kmeans_error)
            push!(ho.obj_kmeans_entropy,entropy)
            push!(ho.obj_kmeans_cross,cross_entropy)

            iter = i
            if i > ho.win_size + 1 && CheckConvergence(ho,0)
                break
            end
        end

        push!(ho.kmeans_rounds,iter)
        push!(ho.obj_harmony,ho.obj_kmeans[end])

        return false
    else
        return true
    end
end

function MoeCorrectRidge!(ho::HarmonyObj)
    ho.z_corr = copy(ho.z_orig)
    @inbounds for i in 1:ho.nclust
        ho.ϕ_rk = ho.ϕ_moe * Diagonal(ho.ρ[i,:])
        W = inv(ho.ϕ_rk * ho.ϕ_moe' .+ ho.λ) * ho.ϕ_rk * ho.z_orig'
        W[1,:] .= 0
        ho.z_corr .-= W' * ho.ϕ_rk
    end
    ho.z_cos = copy(ho.z_corr)
    @inbounds for col in eachcol(ho.z_cos)
        N = norm(col)
        replace!(x -> x / N,col)
    end
end

function Harmonize!(ho::HarmonyObj,
        iter_num::Integer)
    # return 0: complete all iterations
    # return 1: break iterations by convergence
    # return 2: clustering error
    @inbounds for i in 1:iter_num
        @info "Harmonizing $i times..."
        # Step1. clustering
        err_status = Cluster!(ho)
        if err_status
            return 2
        end

        # Step2. regressing out covariates
        MoeCorrectRidge!(ho)

        # Step3. check for convergence
        if CheckConvergence(ho,1)
            return 1
        end
    end
    return 0
end
