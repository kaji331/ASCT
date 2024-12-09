# --- All structures ---
mutable struct NetWork
    nodes_number::Integer
    edges_number::Integer
    node_weights::Vector{Float32}
    first_neighbor_ind::Vector{Int}
    neighbor::Vector{Int}
    edge_weights::Vector{Float32}
    total_weights_self::AbstractFloat
end

mutable struct MClustering
    nodes_number::Integer
    clusters_number::Integer
    clusters::Vector{Int}
end

mutable struct JavaRandom
    seed::UInt
end

mutable struct VOSct
    network::NetWork
    clustering::MClustering
    resolution::AbstractFloat
end

# --- NetWork ---
# 1
function NetWork(nodes_number::Integer,
        node_weights::Vector{Float32},
        first_neighbor_ind::Vector{Int},
        neighbor::Vector{Int},
        edge_weights::Vector{Float32})
    if isempty(node_weights)
        node_weights_tmp = zeros(Float32,nodes_number)
    end
    if isempty(edge_weights)
        edge_weights = ones(Float32,length(neighbor))
    end
    network = NetWork(nodes_number,length(neighbor),node_weights_tmp,
                      first_neighbor_ind,neighbor,edge_weights,0)

    if isempty(node_weights)
        network.node_weights = TotalEdgeWeightsPerNode(network)
    end
    return network
end

function NetWork(nodes_number::Integer,
        first_neighbor_ind::Vector{Int},
        neighbor::Vector{Int},
        edge_weights::Vector{Float32})
    NetWork(nodes_number,Float32[],first_neighbor_ind,neighbor,edge_weights)
end

function TotalEdgeWeights(network::NetWork)
    sum(network.edge_weights) / 2
end

function TotalEdgeWeights(network::NetWork,
        node::Integer)
    idx1 = network.first_neighbor_ind[node]
    idx2 = network.first_neighbor_ind[node + 1]
    s = sum(network.edge_weights[idx1:(idx2 - 1)])
    return s
end

function TotalEdgeWeightsPerNode(network::NetWork)
    total_edge_weights_per_node = zeros(Float32,network.nodes_number)
    map(i -> total_edge_weights_per_node[i] = TotalEdgeWeights(network,i),
        1:network.nodes_number)
    return total_edge_weights_per_node
end

@inline function reduce_weights!(m::Int,
        i::Int,
        j::Base.RefValue,
        network::NetWork,
        clustering::MClustering,
        reduced_network_weights2::Base.RefValue,
        reduced_network_neighbor2::Vector,
        total_weights_self::Base.RefValue)
    n = clustering.clusters[network.neighbor[m]]
    if n != i
        if reduced_network_weights2[][n] == 0
            reduced_network_neighbor2[j[]] = n
            j[] += 1
        end
        reduced_network_weights2[][n] += network.edge_weights[m]
    else
        total_weights_self[] += network.edge_weights[m]
    end
end

@inline function reduce_weights2!(k::Int,
        edges_number::Int,
        reduced_network_neighbor1::Vector,
        reduced_network_neighbor2::Vector,
        reduced_network_weights1::Vector,
        reduced_network_weights2::Base.RefValue)
    reduced_network_neighbor1[edges_number + k] = 
        reduced_network_neighbor2[k]
    reduced_network_weights1[edges_number + k] = 
        reduced_network_weights2[][reduced_network_neighbor2[k]]
    reduced_network_weights2[][reduced_network_neighbor2[k]] = 0
end

@inline function ReducedNet(network::NetWork,
        clustering::MClustering)
    edges_number = 0
    node_weights = zeros(Float32,clustering.clusters_number)
    total_weights_self = Ref{Float32}(network.total_weights_self)
    first_neighbor_ind = ones(Int,clustering.clusters_number + 1)
    reduced_network_neighbor1 = ones(Int,network.edges_number)
    reduced_network_weights1 = zeros(Float32,network.edges_number)
    reduced_network_neighbor2 = ones(Int,clustering.clusters_number - 1)
    reduced_network_weights2 = Ref(zeros(Float32,clustering.clusters_number))
    nodes_per_cluster = NodesPerCluster(clustering)

    @inbounds @fastmath @simd for i in 1:clustering.clusters_number
        j = Ref{Int}(1)
        @inbounds @fastmath @simd for k in 1:length(nodes_per_cluster[i])
            l = nodes_per_cluster[i][k]
            node_weights[i] += network.node_weights[l]
            map(m -> reduce_weights!(m,i,j,network,clustering,
                                     reduced_network_weights2,
                                     reduced_network_neighbor2,
                                     total_weights_self),
                network.first_neighbor_ind[l]:
                (network.first_neighbor_ind[l + 1] - 1))
        end
        map(k -> reduce_weights2!(k,edges_number,reduced_network_neighbor1,
                                  reduced_network_neighbor2,
                                  reduced_network_weights1,
                                  reduced_network_weights2),1:(j[] - 1))
        edges_number += j[] - 1
        first_neighbor_ind[i + 1] = edges_number + 1
    end
    neighbor = reduced_network_neighbor1[begin:(begin + edges_number)]
    edge_weights = reduced_network_weights1[begin:(begin + edges_number)]

    NetWork(clustering.clusters_number,edges_number,node_weights,
            first_neighbor_ind,neighbor,edge_weights,total_weights_self[])
end

# --- MClustering ---
function MClustering(nodes_number::Integer)
    MClustering(nodes_number,1,ones(Int,nodes_number))
end

function NodesNumberPerCluster(clustering::MClustering)
    nodes_number_per_cluster = zeros(Int,clustering.clusters_number)
    @inbounds @fastmath @simd for clust in clustering.clusters
        nodes_number_per_cluster[clust] += 1
    end
    return nodes_number_per_cluster
end

function cluster_position!(i::Integer,
        pos_per_cluster::Vector{Int},
        nodes_per_cluster::Vector{Vector{Int}},
        clustering::MClustering)
    position = pos_per_cluster[clustering.clusters[i]]
    nodes_per_cluster[clustering.clusters[i]][position] = i
    pos_per_cluster[clustering.clusters[i]] += 1
end

function NodesPerCluster(clustering::MClustering)
    nodes_per_cluster = [ Int[] for _ in 1:clustering.clusters_number ]
    nodes_number_per_cluster = NodesNumberPerCluster(clustering)
    pos_per_cluster = ones(Int,clustering.clusters_number)
    @inbounds @fastmath @simd for i in 1:clustering.clusters_number
        resize!(nodes_per_cluster[i],nodes_number_per_cluster[i])
    end
    @inbounds @fastmath @simd for i in 1:clustering.nodes_number
        # position is just a copied number
        position = pos_per_cluster[clustering.clusters[i]]
        nodes_per_cluster[clustering.clusters[i]][position] = i
        pos_per_cluster[clustering.clusters[i]] += 1
    end

    return nodes_per_cluster
end

function InitSingletonClusters!(clustering::MClustering)
    @inbounds @fastmath @simd for i in 1:clustering.nodes_number
        clustering.clusters[i] = i
    end
    clustering.clusters_number = clustering.nodes_number
end

function MergeClusters!(clustering::MClustering,
        clustering_other::MClustering)
    @inbounds @fastmath @simd for i in 1:clustering.nodes_number
        clustering.clusters[i] = 
            clustering_other.clusters[clustering.clusters[i]]
    end
    clustering.clusters_number = clustering_other.clusters_number
end

# --- JavaRandom ---
function Next!(java_random::JavaRandom,
        bits::Integer)
    java_random.seed = ((java_random.seed * UInt(0x5DEECE66D) + UInt(0xB)) & 
                        ((UInt(1) << 48) - 1))
    return Int(java_random.seed >> (48 - bits))
end

function NextInt!(java_random::JavaRandom,
        n::Integer)
    # No need to check if n <= 0
    # Throw error might be slow
    if (n & -n) == n
        return Int((UInt(n) * UInt(Next!(java_random,31))) >> 31)
    end
    bits = nothing
    val = nothing
    bits = Next!(java_random,31)
    val = bits % n
    while bits - val + (n - 1) < 0
        bits = Next!(java_random,31)
        val = bits % n
    end
    
    return val
end

function Seed!(java_random::JavaRandom,
        seed::Integer)
    # \xor
    java_random.seed = (UInt(seed) ⊻ UInt(0x5DEECE66D)) & ((UInt(1) << 48) - 1)
end

# --- VOSClusteringTechnique ---
function VOSct(network::NetWork,
        resolution::AbstractFloat)
    clustering = MClustering(network.nodes_number)
    InitSingletonClusters!(clustering)
    return VOSct(network,clustering,resolution)
end

@inline function quality_fun_sum!(k::Integer,
        j::Integer,
        vosct::VOSct,
        quality_fun::Base.RefValue)
    if vosct.clustering.clusters[vosct.network.neighbor[k]] == j
        quality_fun[] += vosct.network.edge_weights[k]
    end
end

@inline function clust_weights_sum!(i::Integer,
        cluster_weights::Base.RefValue,
        vosct::VOSct)
    cluster_weights[][vosct.clustering.clusters[i]] += 
        vosct.network.node_weights[i]
end

@inline function Quality(vosct::VOSct)
    quality_fun = Ref{Float32}(0.0)
    @inbounds @fastmath @simd for i in 1:vosct.network.nodes_number
        j = vosct.clustering.clusters[i]
        map(k -> quality_fun_sum!(k,j,vosct,quality_fun),
            (vosct.network.first_neighbor_ind[i]:
             vosct.network.first_neighbor_ind[i + 1] - 1))
    end
    quality_fun[] += vosct.network.total_weights_self

    cluster_weights = Ref(zeros(Float32,vosct.clustering.clusters_number))
    map(i -> clust_weights_sum!(i,cluster_weights,vosct),
        1:vosct.network.nodes_number)
    map(i -> quality_fun[] -= cluster_weights[][i] ^ 2 * vosct.resolution,
        1:vosct.clustering.clusters_number)
    quality_fun[] /= 2 * TotalEdgeWeights(vosct.network) + 
        vosct.network.total_weights_self

    return quality_fun[]
end

@inline function neighbor_moving!(k::Int,
        vcc::Base.RefValue,
        vnn::Base.RefValue,
        vnew::Base.RefValue,
        edge_weights_per_cluster::Base.RefValue,
        neighboring_clusters::Base.RefValue,
        neighboring_clusters_number::Base.RefValue)
    l = vcc[][vnn[][k]]
    if edge_weights_per_cluster[][l] == 0
        neighboring_clusters_number[] += 1
        neighboring_clusters[][neighboring_clusters_number[]] = l
    end
    edge_weights_per_cluster[][l] += vnew[][k]
end

@inline function clust_node!(i::Integer,
        vcc::Base.RefValue,
        vnnw::Base.RefValue,
        cluster_weights::Vector,
        nodes_number_per_cluster::Vector)
    cluster_weights[vcc[][i]] += vnnw[][i]
    nodes_number_per_cluster[vcc[][i]] += 1
end

@inline function unused!(i::Integer,
        nodes_number_per_cluster::Vector,
        unused_clusters::Vector,
        unused_clusters_number::Vector)
    if nodes_number_per_cluster[i] == 0
        unused_clusters[unused_clusters_number[1] + 1] = i
        unused_clusters_number[1] += 1
    end
end

@inline function max_qual!(k::Integer,
        j,
        neighboring_clusters::Base.RefValue,
        edge_weights_per_cluster::Base.RefValue,
        vnnw::Base.RefValue,
        cluster_weights::Vector,
        vosct::VOSct,
        max_quality_fun::Base.RefValue,
        best_cluster::Base.RefValue)
    l = neighboring_clusters[][k]
    local quality_fun = edge_weights_per_cluster[][l] - vnnw[][j] * 
        cluster_weights[l] * vosct.resolution
    if quality_fun > max_quality_fun[] || 
        (quality_fun == max_quality_fun[] && l < best_cluster[])
        best_cluster[] = l
        max_quality_fun[] = quality_fun
    end
    edge_weights_per_cluster[][l] = 0
end

@inline function vosct_clust!(i::Integer,
        nodes_number_per_cluster::Vector,
        new_clusters::Vector,
        vosct::VOSct)
    if nodes_number_per_cluster[i] > 0
        new_clusters[i] = vosct.clustering.clusters_number + 1
        vosct.clustering.clusters_number += 1
    end
end

@inline function LocalMoving!(vosct::VOSct,
        random::JavaRandom)
    update = false
    cluster_weights = zeros(Float32,vosct.network.nodes_number)
    nodes_number_per_cluster = zeros(Int,vosct.network.nodes_number)
    vcc = Ref(vosct.clustering.clusters)
    vnnw = Ref(vosct.network.node_weights)

    if vosct.network.nodes_number == 1
        return false
    end
    map(i -> clust_node!(i,vcc,vnnw,cluster_weights,nodes_number_per_cluster),
        1:vosct.network.nodes_number)
    
    unused_clusters_number = [0]
    unused_clusters = zeros(Int,vosct.network.nodes_number)
    map(i -> unused!(i,nodes_number_per_cluster,unused_clusters,
                    unused_clusters_number),
        1:vosct.network.nodes_number)

    nodes_permutation = GenRandomPermutation!(vosct.network.nodes_number,random)
    edge_weights_per_cluster = Ref(zeros(Float32,vosct.network.nodes_number))
    neighboring_clusters = Ref(ones(Int,vosct.network.nodes_number - 1))
    stable_nodes_number = 0
    i = 1
    vnew = Ref(vosct.network.edge_weights)
    vnn = Ref(vosct.network.neighbor)
    while true
        j = nodes_permutation[i]
        neighboring_clusters_number = Ref{Int}(0)

        map(k -> neighbor_moving!(k,vcc,vnn,vnew,edge_weights_per_cluster,
                                  neighboring_clusters,
                                  neighboring_clusters_number),
            vosct.network.first_neighbor_ind[j]:
            (vosct.network.first_neighbor_ind[j + 1] - 1))

        cluster_weights[vcc[][j]] -= vnnw[][j]
        nodes_number_per_cluster[vcc[][j]] -= 1
        if nodes_number_per_cluster[vcc[][j]] == 0
            unused_clusters[unused_clusters_number[1] + 1] = vcc[][j]
            unused_clusters_number[1] += 1
        end
        best_cluster = Ref{Integer}(-1)
        max_quality_fun = Ref(0.)
        map(k -> max_qual!(k,j,neighboring_clusters,edge_weights_per_cluster,
                           vnnw,cluster_weights,vosct,max_quality_fun,
                           best_cluster),1:neighboring_clusters_number[])
        if max_quality_fun[] == 0
            best_cluster[] = unused_clusters[unused_clusters_number[1]]
            unused_clusters_number[1] -= 1
        end
        cluster_weights[best_cluster[]] += vnnw[][j]
        nodes_number_per_cluster[best_cluster[]] += 1
        if best_cluster[] == vcc[][j]
            stable_nodes_number += 1
        else
            vcc[][j] = best_cluster[]
            stable_nodes_number = 1
            update = true
        end
        i = i < vosct.network.nodes_number ? i + 1 : 1
        stable_nodes_number < vosct.network.nodes_number || break
    end

    new_clusters = ones(Int,vosct.network.nodes_number)
    vosct.clustering.clusters_number = 0
    map(i -> vosct_clust!(i,nodes_number_per_cluster,new_clusters,vosct),
        1:vosct.network.nodes_number)
    @inbounds @fastmath @simd for i in 1:vosct.network.nodes_number
        vcc[][i] = new_clusters[vcc[][i]]
    end

    return update
end

# --- Louvain ---
# Recursive Louvain algorithm
function Louvain!(vosct::VOSct,
        random::JavaRandom)
    if vosct.network.nodes_number == 1
        return false
    end
    update = LocalMoving!(vosct,random)
    if vosct.clustering.clusters_number < vosct.network.nodes_number
        vosct_tmp = VOSct(ReducedNet(vosct.network,vosct.clustering),
                          vosct.resolution)
        update2 = Louvain!(vosct_tmp,random)
        if update2
            update = true
            MergeClusters!(vosct.clustering,vosct_tmp.clustering)
        end
    end

    return update
end

# --- Arrays2 ---
@inline function permute!(i::Integer,
        num_ele::Integer,
        java_random::JavaRandom,
        permutation::Vector)
    j = NextInt!(java_random,num_ele)
    k = permutation[i]
    permutation[i] = permutation[j + 1]
    permutation[j + 1] = k
end

function GenRandomPermutation!(num_ele::Integer,
        java_random::JavaRandom)
    permutation = collect(0:(num_ele - 1))
    map(i -> permute!(i,num_ele,java_random,permutation),1:num_ele)
    return permutation .+ 1
end

function MatrixToNetwork(node1::Vector{Int},
        node2::Vector{Int},
        edge_weights::Vector{Float32},
        fun::Integer,
        nodes_number::Integer)

    neighbors_number = zeros(Int,nodes_number)
    @inbounds @fastmath @simd for i in eachindex(node1)
        if node1[i] < node2[i]
            neighbors_number[node1[i]] += 1
            neighbors_number[node2[i]] += 1
        end
    end

    first_neighbor_ind = zeros(Int,nodes_number + 1)
    edges_number = 0
    @inbounds @fastmath @simd for i in 1:nodes_number
        first_neighbor_ind[i] = edges_number + 1
        edges_number += neighbors_number[i]
    end
    first_neighbor_ind[nodes_number + 1] = edges_number + 1

    neighbor = zeros(Int,edges_number)
    edge_weights2 = zeros(Float32,edges_number)
    neighbors_number = zeros(Int,nodes_number)
    @inbounds @fastmath @simd for i in eachindex(node1)
        if node1[i] < node2[i]
            j = first_neighbor_ind[node1[i]] + neighbors_number[node1[i]]
            neighbor[j] = node2[i]
            edge_weights2[j] = edge_weights[i]
            neighbors_number[node1[i]] += 1
            j = first_neighbor_ind[node2[i]] + neighbors_number[node2[i]]
            neighbor[j] = node1[i]
            edge_weights2[j] = edge_weights[i]
            neighbors_number[node2[i]] += 1
        end
    end

    if fun == 1
        return NetWork(nodes_number,first_neighbor_ind,neighbor,edge_weights2)
    else
        @warn "Now only support 1 for modularity function!"
        return NetWork(nodes_number,first_neighbor_ind,neighbor,edge_weights2)
    end
end
