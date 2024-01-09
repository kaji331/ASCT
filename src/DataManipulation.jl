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

function SaveSeuratV4(obj::WsObj,
        file::String;
        project::String = "ASCT")
    h5 = h5open(file,"w")

    # Attributions
    attributes(h5)["active.assay"] = ["RNA"]
    attributes(h5)["project"] = [project]
    attributes(h5)["version"] = ["4.0.0"]

    # Datasets

    # Create active.ident
    h5["active.ident/levels"] = [project]
    h5["active.ident/values"] = repeat([1],obj.dat["raw_dat"].n)

    # Create assays
    h5["assays/RNA/counts/data"] = obj.dat["raw_dat"].nzval
    h5["assays/RNA/counts/indices"] = obj.dat["raw_dat"].rowval .- 1
    h5["assays/RNA/counts/indptr"] = obj.dat["raw_dat"].colptr .- 1
    attributes(h5["assays/RNA"])["key"] = ["rna_"]
    attributes(h5["assays/RNA/counts"])["dims"] = 
        [obj.dat["raw_dat"].m,obj.dat["raw_dat"].n]

    if "norm_dat" in keys(obj.dat)
        h5["assays/RNA/data/data"] = obj.dat["norm_dat"].nzval
        h5["assays/RNA/data/indices"] = obj.dat["norm_dat"].rowval .- 1
        h5["assays/RNA/data/indptr"] = obj.dat["norm_dat"].colptr .- 1
        attributes(h5["assays/RNA/data"])["dims"] = 
            [obj.dat["norm_dat"].m,obj.dat["norm_dat"].n]
    else
        h5["assays/RNA/data/data"] = obj.dat["raw_dat"].nzval
        h5["assays/RNA/data/indices"] = obj.dat["raw_dat"].rowval .- 1
        h5["assays/RNA/data/indptr"] = obj.dat["raw_dat"].colptr .- 1
        attributes(h5["assays/RNA/data"])["dims"] = 
            [obj.dat["raw_dat"].m,obj.dat["raw_dat"].n]
    end

    gn = replace.(obj.var.name,r"_" => s"-")
    h5["assays/RNA/features"] = gn

    if "scale_dat" in keys(obj.dat)
        h5["assays/RNA/scale.data"] = obj.dat["scale_dat"]
        h5["assays/RNA/scaled.features"] = gn
    end

    create_group(h5,"assays/RNA/misc")

    # Create cell.names
    h5["cell.names"] = obj.obs.barcodes

    # Create commands
    create_group(h5,"commands")
    current_time = now() |> string |> splitext |> 
        x -> replace(x[1],"T" => " ")
    for command in obj.log
        # assay_used = match(r"ptr=(.*?)\)",command)
        # assay_used = isnothing(assay_used) ? "NA" : assay_used.captures[1]
        assay_used = "RNA"
        call_string = command
        name = match(r"(.*?)\((.*)",command).captures[1]
        para = findall(r"([ a-zA-Z_\-0-9]+=[ a-zA-Z\"_\-0-9:\.]+)",command) |> 
            x -> [ command[m] for m in x ] |> x -> split.(x,"=")
        create_group(h5,"commands/" * name)
        attributes(h5["commands/" * name])["assay.used"] = [assay_used]
        attributes(h5["commands/" * name])["call.string"] = [call_string]
        attributes(h5["commands/" * name])["name"] = [name]
        names = isempty(para) ? [""] : [ n[1] for n in para ]
        attributes(h5["commands/" * name])["names"] = names
        attributes(h5["commands/" * name])["time.stamp"] = [current_time]
        for p in para
            h5["commands/" * name * "/" * p[1]] = [p[2]]
            if p[2] in ["true","false"]
                attributes(h5["commands/" * name * "/" * p[1]])["s3class"] = 
                    "logical"
            end
        end
    end

    # Create graphs
    create_group(h5,"graphs")

    # Create images
    create_group(h5,"images")

    # Create misc
    create_group(h5,"misc")

    # Create neighbors
    create_group(h5,"neighbors")

    # Create reductions
    create_group(h5,"reductions")

    # Create tools
    create_group(h5,"tools")

    # Create meta.data
    create_group(h5,"meta.data")
    attributes(h5["meta.data"])["_index"] = ["_index"]
    meta_colnames1 = ["nCount_RNA","nFeature_RNA"]
    meta_colnames2 = names(obj.obs)[findall(x -> x ∉ ["barcodes",
                                                      "cell_counts",
                                                      "cell_features"],
                                            names(obj.obs))]
    h5["meta.data/orig.ident/levels"] = [project]
    h5["meta.data/orig.ident/values"] = repeat([1],obj.dat["raw_dat"].n)
    h5["meta.data/nCount_RNA"] = obj.obs.cell_counts
    h5["meta.data/nFeature_RNA"] = obj.obs.cell_features
    for m in meta_colnames2
        h5["meta.data/" * m] = obj.obs[!,m]
    end
    if isnothing(obj.meta)
        attributes(h5["meta.data"])["colnames"] = 
            vcat("orig.ident",meta_colnames1,meta_colnames2)
        close(h5)
        return "Finished!"
    else
        x = keys(obj.meta) |> x -> string.(x)
        i = @. !isnothing(match(r"^clusters_",x))
        if !all(i .== 0)
            meta_colnames3 = x[i]
        end
        attributes(h5["meta.data"])["colnames"] = 
            vcat("orig.ident",meta_colnames1,meta_colnames2,meta_colnames3)
    end
    for m in meta_colnames3
        if typeof(obj.meta[m]) == DataFrame
            # DEGs
            create_group(h5,"misc/" * m)
            for col in names(obj.meta[m])
                h5["misc/" * m * "/" * col] = obj.meta[m][!,col]
            end
            attributes(h5["misc/$m"])["colnames"] = names(obj.meta[m])
        else
            # Clusters
            create_group(h5,"meta.data/" * m)
            h5["meta.data/" * m * "/levels"] = 
                obj.meta[m] |> unique |> sort |> x -> string.(x)
            h5["meta.data/" * m * "/values"] = obj.meta[m]
        end
    end

    # Other info along with meta
    # ======

    # Fill HVGs
    if "hvg_index" in keys(obj.meta)
        h5["assays/RNA/meta.features/_index"] = gn
        h5["assays/RNA/meta.features/vst.mean"] = obj.meta["hvg_mean"]
        h5["assays/RNA/meta.features/vst.variance.standardized"] = 
            obj.meta["hvg_var_std"]
    end

    if "hvg_index" in keys(obj.meta)
        h5["assays/RNA/variable.features"] = 
            replace.(obj.meta["hvg_name"],r"_" => s"-")
    end

    # Fill reductions
    if "pca" in keys(obj.meta)
        create_group(h5,"reductions/pca")
        attributes(h5["reductions/pca"])["active.assay"] = ["RNA"]
        attributes(h5["reductions/pca"])["global"] = [0]
        attributes(h5["reductions/pca"])["key"] = ["PC_"]
        h5["reductions/pca/cell.embeddings"] = obj.meta["pca"]
        h5["reductions/pca/features"] = obj.meta["hvg_name"]
        h5["reductions/pca/misc/pca_cut"] = [obj.meta["pca_cut"]]
        attributes(h5["reductions/pca/misc"])["names"] = ["pca_cut"]
        h5["reductions/pca/pca_var"] = obj.meta["pca_var"]
    end

    if "tsne" in keys(obj.meta)
        create_group(h5,"reductions/tsne")
        attributes(h5["reductions/tsne"])["active.assay"] = ["RNA"]
        attributes(h5["reductions/tsne"])["global"] = [1]
        attributes(h5["reductions/tsne"])["key"] = ["TSNE_"]
        h5["reductions/tsne/cell.embeddings"] = obj.meta["tsne"]
        create_group(h5,"reductions/tsne/misc")
    end

    if "umap" in keys(obj.meta)
        create_group(h5,"reductions/umap")
        attributes(h5["reductions/umap"])["active.assay"] = ["RNA"]
        attributes(h5["reductions/umap"])["global"] = [1]
        attributes(h5["reductions/umap"])["key"] = ["UMAP_"]
        h5["reductions/umap/cell.embeddings"] = obj.meta["umap"]
        create_group(h5,"reductions/umap/misc")
    end

    if "harmony" in keys(obj.meta)
        create_group(h5,"reductions/harmony")
        attributes(h5["reductions/harmony"])["active.assay"] = ["RNA"]
        attributes(h5["reductions/harmony"])["global"] = [1]
        attributes(h5["reductions/harmony"])["key"] = ["Harmony_"]
        h5["reductions/harmony/cell.embeddings"] = obj.meta["harmony"]
        create_group(h5,"reductions/harmony/misc")
    end

    close(h5)

    return "Finished!"
end
