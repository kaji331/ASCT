function MergeRawData(objs::Vector{WsObj},
        grp::Dict{AbstractString,AbstractVector};
        min_features::Union{Nothing,Integer} = nothing,
        min_cells::Union{Nothing,Integer} = nothing)
    # Merge dat using DataFrame. Need to optimize in future!
    @info "Data merging..."
    all_dat = [ hcat(DataFrame(gene=obj.var.name),
                     DataFrame(id=obj.var.id),
                     DataFrame(Matrix(obj.dat["raw_dat"]),
                               obj.obs.barcodes .* "_$i")) 
               for (i,obj) in enumerate(objs) ] |> 
        x -> outerjoin(x...,on=:gene,validate=(true,true),makeunique=true)
    id = match.(r"^id",names(all_dat)) |> x -> all_dat[:,x .!= nothing]
    id = [ unique(id[i,:]) |> x -> join(x[typeof.(x) .!= Missing],"_") 
          for i in 1:size(id,1) ]
    all_dat = match.(r"^id",names(all_dat)) |> x -> all_dat[:,x .== nothing]
    all_dat.id = id
    @inbounds @simd for col in eachcol(all_dat)
        replace!(col,missing => 0)
    end
    dat = all_dat[!,Not(:gene,:id)] |> Matrix |> SparseMatrixCSC{Int32,Int32}
    features = Dict{AbstractString,
                    AbstractVector}("id" => String.(all_dat.id),
                                    "name" => String.(all_dat.gene))
    barcodes = names(all_dat)[2:(dat.n + 1)]
    all_dat = nothing
    GC.gc()

    # Create WsObj
    # Basic calculation
    @info "Creating WsObj..."
    cell_counts = dat' * ones(Int,size(dat,1))
    cell_features = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
    dat = dat' |> sparse
    cells = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
    features["feature_cells"] = cells
    not_zero = cells .!= 0
    for k in keys(features)
        features[k] = features[k][not_zero]
    end
    push!(grp,
          "barcodes" => barcodes,
          "cell_counts" => cell_counts,
          "cell_features" => cell_features)
    obj = WsObj(Dict("raw_dat" => dat[:,not_zero]' |> sparse),
                DataFrame(grp),
                DataFrame(features),
                String[],
                nothing)

   if !(isnothing(min_features) && isnothing(min_cells))
       SimpleFilter!(obj,min_features,min_cells)
   end

   return obj
end
