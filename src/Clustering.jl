# We provide graph-based and kmedoids methods for clustering.

function Clustering!(obj::WsObj;
        method::AbstractString = "mc",
        reduction::Union{AbstractString,Symbol} = :auto,
        use_pc::Union{AbstractString,Integer} = "pca_cut",
        tree_K::Integer = 20,
        resolution::Union{Symbol,Real,AbstractRange} = :auto,
        cluster_K::Union{Nothing,Integer} = nothing,
        cluster_K_max::Union{Nothing,Integer} = 30,
        dist::AbstractString = "Euclidean",
        network::AbstractString = "SNN",
        random_starts_number::Integer = 10,
        iter_number::Integer = 10,
        prune::AbstractFloat = 1 /  15,
        seed::Integer = -1)

    # Support PCA or Harmony
    if reduction == :auto
        if "harmony" in keys(obj.meta)
            reduction = "harmony"
        else
            reduction = "pca"
        end
    end

    # Check prune
    if prune > 1 || prune < 0
        throw(DomainError(prune,"Must be 0 ~ 1."))
    end
    # Check method
    # Leiden package has some problems. It should be re-written in future!
    if !(method in ["mc","kmedoids"])
        return "Nothing to do! 'method' is only support \"mc\"|\"kmedoids\"!"
    end

    # Check PCs
    if typeof(use_pc) <: AbstractString
        use_pc = obj.meta[use_pc]
    end

    # Check resolution
    if resolution == :auto
        resolution = 0.2:0.1:2.0
    elseif typeof(resolution) <: Real
        resolution = resolution
    end

    # Do clustering
    if method == "mc"
        ModularityClustering!(obj;
                              reduction = reduction,
                              num_pc = use_pc,
                              network = network,
                              K = tree_K,
                              resolution = resolution,
                              prune = prune,
                              random_starts_number = random_starts_number,
                              iter_number = iter_number,
                              seed = seed)
    else
        KClustering!(obj;
                     reduction = reduction,
                     num_pc = use_pc,
                     dist_method = dist,
                     K = cluster_K,
                     max_K = cluster_K_max,
                     seed = seed)
    end

    # log
    push!(obj.log,String("Clustering!(WsObj;" * 
                         "method = $method," * 
                         "reduction = $reduction," * 
                         "use_pc = $use_pc," * 
                         "tree_K = $tree_K," * 
                         "resolution = $resolution," * 
                         "cluster_K = $cluster_K," * 
                         "cluster_K_max = $cluster_K_max," *
                         "dist = $dist," * 
                         "network = $network," * 
                         "random_starts_number = $random_starts_number," * 
                         "iter_number = $iter_number," * 
                         "prune = $prune," * 
                         "seed = $seed)"))

    return "Finished!"
end

function KClustering!(obj::WsObj;
        reduction::AbstractString = "pca",
        num_pc::Integer = 5,
        dist_method::AbstractString = "CorrDist",
        K::Union{Nothing,Integer} = nothing,
        max_K::Union{Nothing,Integer} = 50,
        seed::Integer = -1)

    if dist_method in ["Euclidean",
                       "SqEuclidean",
                       "Cityblock",
                       "Chebyshev",
                       "Jaccard",
                       "CosineDist",
                       "CorrDist",
                       "RMSDeviation"]
        if reduction == "pca"
            dist = pairwise(getfield(Main,Symbol(dist_method))(),
                            obj.meta["pca"][:,1:num_pc]')
        else
            dist = pairwise(getfield(Main,Symbol(dist_method))(),
                            obj.meta[reduction]')
        end
    end

    if isnothing(K) && isnothing(max_K)
        return "There must be a \"cluster_K\" or \"cluster_K_max\"!"
    end

    if !isnothing(K)
        Random.seed!(seed < 0 ? 1984 : seed)
        cell_assignments = kmedoids(dist,K).assignments
    else
        # Automatically choose resolution, based the silhouettes coefficients 
        # and elbow point by total costs
        clustering_vector = Vector{KmedoidsResult{Float64}}()
        @inbounds @simd for k in 2:max_K
            Random.seed!(seed < 0 ? 1984 : seed)
            push!(clustering_vector,kmedoids(dist,k))
        end

        # Check silhouettes coefficients
        sil_check = Vector{Tuple}()
        @inbounds for c in clustering_vector
            sh = silhouettes(c,dist)
            sh_avg = mean(sh)
            sh_std = std(sh)
            check_max = Vector{Bool}()
            check_min = Vector{Bool}()
            @simd for i in unique(c.assignments)
                idx = findall(x -> x == i,c.assignments)
                # Check all clusters have the max coefficient more than whole 
                # average
                push!(check_max,maximum(sh[idx]) > sh_avg)
                # Check all clusters have no coefficient less than zero
                push!(check_min,minimum(sh[idx]) < 0)
            end
            check_max = all(check_max)
            check_min = all(check_min)
            push!(sil_check,(sh_avg,sh_std,check_max,check_min))
        end
        # Order and find the best K
        idx = findall(x -> all(x[3:4]),sil_check)
        res = [idx .+ 1 [ sil_check[i][1] 
                         for i in idx ] [ sil_check[i][2] 
                                         for i in idx ]]
        sh_coef_order = res[sortperm(res[:,2]),1]
        sh_std_order = res[sortperm(res[:,3];rev=true),1]
        # Find the top order of average of silhouettes coefficients and bottom 
        # order of standard deviation of silhouettes coefficients
        silhouettes_K = [ (i + findall(x -> x == sh_coef_order[i],
                                       sh_std_order)[1],
                           sh_coef_order[i]) 
                         for i in eachindex(sh_coef_order) ] |> 
            x -> x[findmax([ i[1] for i in x])[2]][2] |> 
            Int

        # Find elbow of total costs
        tc_fit = glm(@formula(y ~ x + x ^ 2 + x ^ 3),
                     DataFrame(x=2.0:max_K,
                               y=[c.totalcost for c in clustering_vector ]),
                     Gamma(),
                     InverseLink())
        tc_fit2 = lm(@formula(y ~ x),
                     DataFrame(x=[2,max_K],
                               y=[clustering_vector[1].totalcost,
                                  clustering_vector[end].totalcost]))
        totalcost_K = findmax(GLM.predict(tc_fit2,DataFrame(x=2:max_K)) |> 
                                x -> x - [c.totalcost 
                                          for c in clustering_vector ])[2]

        K = round(Int,(silhouettes_K + totalcost_K) / 2)
        @info "Automatically choose K: $K"
        cell_assignments = clustering_vector[K - 1].assignments
    end

    SortAssignments!(cell_assignments)
    # Output
    obj.meta["clusters_" * dist_method * "_" * string(K)] = cell_assignments
    obj.meta["clusters_latest"] = cell_assignments
end

function ModularityClustering!(obj::WsObj;
        modularity_fun::Integer = 1,
        algorithm::Integer = 1,
        reduction::AbstractString = "pca",
        num_pc::Integer = 5,
        network::AbstractString = "SNN",
        K::Integer = 20,
        resolution::Union{AbstractRange,Real} = 0.2:0.1:1.5,
        prune::AbstractFloat = 1 / 15,
        random_starts_number::Integer = 10,
        iter_number::Integer = 10,
        seed::Integer = -1)

    if !(modularity_fun in [1,2])
        return "Nothing to do! 'modularity_fun' must be 1 or 2!"
    end
    if !(algorithm in [1,2,3,4])
        return "Nothing to do! 'algorithm' must be 1, 2, 3 or 4!"
    end
    if random_starts_number < 1
        return "Nothing to do! " * 
            "'random_starts_number must be more than or equal to 1!"
    end
    if iter_number < 1
        return "Nothing to do!, 'iter_number' must be more than or equal to 1!"
    end
    if any(resolution .> 1) && modularity_fun == 2
        return "Nothing to do! 'resolution' must be less than or equal to 1, " * 
            "when 'modularity_fun' is 2"
    end

    seed = seed < 0 ? 1984 : seed

    if reduction == "pca"
        d = @view obj.meta["pca"][:,1:num_pc]
    else
        d = obj.meta[reduction]
    end

    tree = KDTree(d')
    res = knn(tree,d',K,true)[1]

    if network == "KNN"
        # KNN clustering is faster than SNN clustering!

        # # From Seurat, faster but not adjacent/symmetric matrix!
        # # Clustering result may be not good!
        # i = vcat(res...)
        # knn_matrix = sparse((eachindex(i) . - 1) .÷ K .+ 1,i,1)

        # Adjacent matrix, but slower!
        # Clustering result may be better!
        g = SimpleWeightedGraph(size(d,1))
        @inbounds for v in vertices(g)
            @simd for k in 1:K
                add_edge!(g,v,res[v][k],1)
            end
        end
        m = g.weights
    elseif network == "SNN"
        m = SNN(hcat(res...)',prune)
    else
        return "Nothing to do! 'network' must be 'KNN' or 'SNN'!"
    end

    # Check empty network
    empty_network = true
    @inbounds for i in eachindex(m)
        if m[i] != 0
            empty_network = false
            break
        end
    end
    if empty_network
        return "Nothing to do! This matrix has no network inside!"
    end

    if typeof(resolution) <: Real
        cell_assignments = (resolution,
                            ModClustering(m;
                                          resolution=resolution,
                                          rsn=random_starts_number,
                                          itn=iter_number,
                                          seed=seed,
                                          modularity_fun=modularity_fun,
                                          algorithm=algorithm).clusters)
    else
        clustering_vector = Tuple{Float32,Vector{Int}}[]
        slice_time = 0
        @inbounds for i in eachindex(resolution)
            push!(clustering_vector,
                  (resolution[i],
                   ModClustering(m;
                                 resolution=resolution[i],
                                 rsn=random_starts_number,
                                 itn=iter_number,
                                 seed=seed,
                                 modularity_fun=modularity_fun,
                                 algorithm=algorithm).clusters))
            if i == 1
                continue
            end
            slice_number = 0
            @simd for j in unique(clustering_vector[i - 1][2])
                idx = findall(x -> x == j,clustering_vector[i - 1][2])
                counter = StatsBase.counts(clustering_vector[i][2][idx])
                percentage = [ i / sum(counter) for i in counter ]
                # # Frankly, I think the commited condition is reasonable!
                # if isnothing(findfirst(x -> x > 0.99,percentage))
                if !isnothing(findfirst(x -> x > 0.1,percentage))
                    slice_number += 1
                end
            end
            if slice_number >= 2 && 
                length(unique(clustering_vector[i][2])) == 
                    length(unique(clustering_vector[i - 1][2]))
                slice_time += 1
            end
            if slice_time == 2
                @info "Recommended resolution is $(resolution[i - 1])"
                break
            end
        end
        pos = length(clustering_vector) - 1
        cell_assignments = (clustering_vector[pos][1],
                            clustering_vector[pos][2])
    end

    SortAssignments!(cell_assignments[2])
    # Output
    obj.meta["clusters_MC_" * string(cell_assignments[1])] = cell_assignments[2]
    obj.meta["clusters_latest"] = cell_assignments[2]
end

function SortAssignments!(ca::AbstractVector)
    x = Int[]
    m = countmap(ca)
    for i in 2:length(m)
        mm = findmax(m)[2]
        push!(x,mm)
        pop!(m,mm)
    end
    idx = [ ca .== v for v in x ]
    for (i,id) in enumerate(idx)
        ca[id] .= i
    end
end
