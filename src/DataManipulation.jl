"""
    MergeRawData(objs,grp)

Merge raw counts from different WsObj structs into a single WsObj.

# Arguments
- `objs::Vector{WsObj}`: a vector containing different WsObj structs.
- `grp::Dict{<: AbstractString,<: AbstractVector}`: a dictionary. The key means 
  a column name of observation, and the vector means the labels of different 
  datasets.

# Keyword Arguments
- `min_features::Union{Nothing,Integer} = nothing`: drop cells containing 
  features less than this number.
- `min_cells::Union{Nothing,Integer} = nothing`: drop features containing 
  cells less than this number.

# Example
```julia-repl
julia> new_data = MergeRawData([data_200cells,data_300cells],
                               Dict("batch" => vcat(repeat(["1"],200),
                                                    repeat(["2"],300))))
```
"""
function MergeRawData(objs::Vector{WsObj},
        grp::Dict{<: AbstractString,<: AbstractVector};
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
    id = match.(r"^id",names(all_dat)) |> x -> all_dat[!,x .!= nothing]
    id = [ unique(id[i,:]) |> x -> join(x[typeof.(x) .!= Missing],"_") 
          for i in 1:size(id,1) ]
    all_dat = match.(r"^id",names(all_dat)) |> x -> all_dat[!,x .== nothing]
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

"""
    SaveSeuratV4(obj,file)

Save a WsObj to the h5seurat file.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.
- `file::AbstractString`: file name and path.
"""
function SaveSeuratV4(obj::WsObj,
        file::AbstractString;
        project::AbstractString = "ASCT")

    # Check suffix
    if splitext(file)[2] != ".h5seurat"
        @warn "You must save data to a h5seurat file!"
        return nothing
    end

    h5 = h5open(file,"w")

    # Attributions
    HDF5.attributes(h5)["active.assay"] = ["RNA"]
    HDF5.attributes(h5)["project"] = [project]
    HDF5.attributes(h5)["version"] = ["4.0.0"]

    # Datasets

    # Create active.ident
    h5["active.ident/levels"] = [project]
    h5["active.ident/values"] = repeat([1],obj.dat["raw_dat"].n)

    # Create assays
    h5["assays/RNA/counts/data"] = obj.dat["raw_dat"].nzval
    h5["assays/RNA/counts/indices"] = obj.dat["raw_dat"].rowval .- 1
    h5["assays/RNA/counts/indptr"] = obj.dat["raw_dat"].colptr .- 1
    HDF5.attributes(h5["assays/RNA"])["key"] = ["rna_"]
    HDF5.attributes(h5["assays/RNA/counts"])["dims"] = 
        [obj.dat["raw_dat"].m,obj.dat["raw_dat"].n]

    if "norm_dat" in keys(obj.dat)
        h5["assays/RNA/data/data"] = obj.dat["norm_dat"].nzval
        h5["assays/RNA/data/indices"] = obj.dat["norm_dat"].rowval .- 1
        h5["assays/RNA/data/indptr"] = obj.dat["norm_dat"].colptr .- 1
        HDF5.attributes(h5["assays/RNA/data"])["dims"] = 
            [obj.dat["norm_dat"].m,obj.dat["norm_dat"].n]
    else
        h5["assays/RNA/data/data"] = obj.dat["raw_dat"].nzval
        h5["assays/RNA/data/indices"] = obj.dat["raw_dat"].rowval .- 1
        h5["assays/RNA/data/indptr"] = obj.dat["raw_dat"].colptr .- 1
        HDF5.attributes(h5["assays/RNA/data"])["dims"] = 
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
        para = findall(r"([ a-zA-Zα-ωΑ-Ω_\-0-9]+=[ a-zA-Zα-ωΑ-Ω\"_\-0-9:\.]+)",
                       command) |> 
            x -> [ command[m] for m in x ] |> x -> split.(x,"=") |> 
            x -> [ strip.(e) for e in x ]
        create_group(h5,"commands/" * name)
        HDF5.attributes(h5["commands/" * name])["assay.used"] = [assay_used]
        HDF5.attributes(h5["commands/" * name])["call.string"] = [call_string]
        HDF5.attributes(h5["commands/" * name])["name"] = [name]
        names = isempty(para) ? [""] : [ n[1] for n in para ]
        HDF5.attributes(h5["commands/" * name])["names"] = names
        HDF5.attributes(h5["commands/" * name])["time.stamp"] = [current_time]
        for p in para
            h5["commands/" * name * "/" * p[1]] = [p[2]]
            if p[2] in ["true","false"]
                HDF5.attributes(h5["commands/" * name * "/" * 
                                   p[1]])["s3class"] = "logical"
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
    HDF5.attributes(h5["meta.data"])["_index"] = ["_index"]
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
        if eltype(obj.obs[!,m]) <: AbstractFloat
            h5["meta.data/" * m] = obj.obs[!,m]
        else
            create_group(h5,"meta.data/" * m)
            h5["meta.data/" * m * "/levels"] = 
                obj.obs[!,m] |> unique |> sort |> x -> string.(x)
            h5["meta.data/" * m * "/values"] = obj.obs[!,m]
        end
    end
    if isnothing(obj.meta)
        HDF5.attributes(h5["meta.data"])["colnames"] = 
            vcat("orig.ident",meta_colnames1,meta_colnames2)
        close(h5)
        return "Finished!"
    end

    # Other info along with meta
    # ======

    # Fill misc
    k = string.(keys(obj.meta))
    # DEGs
    for i in k[findall(x -> occursin(r"_DE$",x),k)]
        create_group(h5,"misc/" * i)
        for col in names(obj.meta[i])
            h5["misc/" * i * "/" * col] = obj.meta[i][!,col]
        end
        HDF5.attributes(h5["misc/" * i])["colnames"] = names(obj.meta[i])
    end

    # Fill HVGs
    if "hvg_index" in keys(obj.meta)
        h5["assays/RNA/meta.features/_index"] = gn
        h5["assays/RNA/meta.features/vst.mean"] = obj.var.hvg_mean
        h5["assays/RNA/meta.features/vst.variance.standardized"] = 
            obj.var.hvg_var_std
    end

    if "hvg_index" in keys(obj.meta)
        h5["assays/RNA/variable.features"] = 
            replace.(obj.meta["hvg_name"],r"_" => s"-")
    end

    # Fill reductions
    if "pca" in keys(obj.meta)
        create_group(h5,"reductions/pca")
        HDF5.attributes(h5["reductions/pca"])["active.assay"] = ["RNA"]
        HDF5.attributes(h5["reductions/pca"])["global"] = [0]
        HDF5.attributes(h5["reductions/pca"])["key"] = ["PC_"]
        h5["reductions/pca/cell.embeddings"] = obj.meta["pca"]
        h5["reductions/pca/features"] = obj.meta["hvg_name"]
        h5["reductions/pca/misc/pca_cut"] = [obj.meta["pca_cut"]]
        HDF5.attributes(h5["reductions/pca/misc"])["names"] = ["pca_cut"]
        h5["reductions/pca/pca_var"] = obj.meta["pca_var"]
    end

    if "tsne" in keys(obj.meta)
        create_group(h5,"reductions/tsne")
        HDF5.attributes(h5["reductions/tsne"])["active.assay"] = ["RNA"]
        HDF5.attributes(h5["reductions/tsne"])["global"] = [1]
        HDF5.attributes(h5["reductions/tsne"])["key"] = ["TSNE_"]
        h5["reductions/tsne/cell.embeddings"] = obj.meta["tsne"]
        create_group(h5,"reductions/tsne/misc")
    end

    if "umap" in keys(obj.meta)
        create_group(h5,"reductions/umap")
        HDF5.attributes(h5["reductions/umap"])["active.assay"] = ["RNA"]
        HDF5.attributes(h5["reductions/umap"])["global"] = [1]
        HDF5.attributes(h5["reductions/umap"])["key"] = ["UMAP_"]
        h5["reductions/umap/cell.embeddings"] = obj.meta["umap"]
        create_group(h5,"reductions/umap/misc")
    end

    if "harmony" in keys(obj.meta)
        create_group(h5,"reductions/harmony")
        HDF5.attributes(h5["reductions/harmony"])["active.assay"] = ["RNA"]
        HDF5.attributes(h5["reductions/harmony"])["global"] = [1]
        HDF5.attributes(h5["reductions/harmony"])["key"] = ["Harmony_"]
        h5["reductions/harmony/cell.embeddings"] = obj.meta["harmony"]
        create_group(h5,"reductions/harmony/misc")
    end

    close(h5)

    return "Finished!"
end

"""
    SaveAnnData(obj,file)

Save a WsObj to the h5ad file.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.
- `file::AbstractString`: file name and path.
"""
function SaveAnnData(obj::WsObj,
        file::AbstractString)

    # Check suffix
    if splitext(file)[2] != ".h5ad"
        @warn "You must save data to a h5ad file!"
        return nothing
    end

    h5 = h5open(file,"w")
    HDF5.attributes(h5)["encoding-type"] = "anndata"
    HDF5.attributes(h5)["encoding-version"] = "0.1.0"

    # X
    if "scale_dat" in keys(obj.dat)
        @info "Following scanpy's thoughts, it will drop UMI counts..."
        h5["X"] = obj.dat["scale_dat"]
        EnTypeVersion!(h5,"X","array","0.2.0")
        d = obj.dat["norm_dat"]' |> sparse
        pos = "raw"
    else
        if "norm_dat" in keys(obj.dat)
            @info "Following scanpy's thoughts, it will drop UMI counts..."
            d = obj.dat["norm_dat"]' |> sparse
        else
            d = obj.dat["raw_dat"]' |> sparse
        end
        pos = ""
    end

    h5[pos * "/X/data"] = d.nzval
    h5[pos * "/X/indices"] = d.rowval .- 1
    h5[pos * "/X/indptr"] = d.colptr .- 1
    EnTypeVersion!(h5,pos,"raw","0.1.0")
    EnTypeVersion!(h5,pos * "/X","csc_matrix","0.1.0")
    HDF5.attributes(h5[pos * "/X"])["shape"] = [d.m,d.n]

    # varm(s)
    create_group(h5,pos * "/varm")
    EnTypeVersion!(h5,pos * "/varm","dict","0.1.0")
    if pos == "raw"
        # Drop pca of HVGs to varm/PCs
        create_group(h5,"varm")
        EnTypeVersion!(h5,"varm","dict","0.1.0")
    end

    # var(s)
    for i in 1:2
        h5[pos * "/var/_index"] = obj.var.name
        EnTypeVersion!(h5,pos * "/var/_index","string-array","0.2.0")
        h5[pos * "/var/gene_ids"] = obj.var.id
        EnTypeVersion!(h5,pos * "/var/gene_ids","string-array","0.2.0")

        column_order = replace(names(obj.var),
                               "id" => "gene_ids","name" => "_index")

        # Remove 'hvg_mean' and 'hvg_var_std' from raw/var
        idx = findall(x -> x in ["hvg_mean","hvg_var_std"],column_order)
        deleteat!(column_order,idx)

        HDF5.attributes(h5[pos * "/var"])["_index"] = "_index"
        HDF5.attributes(h5[pos * "/var"])["column-order"] = column_order
        EnTypeVersion!(h5,pos * "/var","dataframe","0.2.0")

        for i in column_order[findall(x -> x ∉ ["gene_ids","_index"],
                                      column_order)]
            h5[pos * "/var/" * i] = obj.var[!,i]
            EnTypeVersion!(h5,pos * "/var/" * i,"array","0.2.0")
        end


        if pos != ""
            pos = ""
        else
            if "hvg_mean" in names(obj.var)
                h5["var/mean"] = obj.var.hvg_mean
                EnTypeVersion!(h5,"var/mean","array","0.2.0")
                h5["var/var_std"] = obj.var.hvg_var_std
                EnTypeVersion!(h5,"var/var_std","array","0.2.0")
            end
            break
        end
    end

    # layers
    create_group(h5,"layers")
    EnTypeVersion!(h5,"layers","dict","0.1.0")

    # obs
    h5["obs/_index"] = obj.obs.barcodes
    EnTypeVersion!(h5,"obs/_index","string-array","0.2.0")

    column_order = replace(names(obj.obs),"barcodes" => "_index")
    if "clusters_latest" in names(obj.obs)
        categories = (obj.obs.clusters_latest .- 1) |> 
            unique |> 
            sort |> 
            x -> string.(x)
        codes = obj.obs.clusters_latest .- 1
    end
    HDF5.attributes(h5["obs"])["_index"] = "_index"
    HDF5.attributes(h5["obs"])["column-order"] = column_order
    EnTypeVersion!(h5,"obs","dataframe","0.2.0")

    for i in column_order[findall(x -> x ∉ ["_index"],column_order)]
        if i == "clusters_latest"
            h5["obs/clusters_latest/categories"] = categories
            h5["obs/clusters_latest/codes"] = codes
            EnTypeVersion!(h5,"obs/clusters_latest","categorical","0.2.0")
            HDF5.attributes(h5["obs/clusters_latest"])["ordered"] = 0
            EnTypeVersion!(h5,"obs/clusters_latest/categories","string-array",
                           "0.2.0")
            EnTypeVersion!(h5,"obs/clusters_latest/codes","array","0.2.0")
        else
            h5["obs/" * i] = obj.obs[!,i]
            EnTypeVersion!(h5,"obs/" * i,"array","0.2.0")
        end
    end

    # obsm
    create_group(h5,"obsm")
    EnTypeVersion!(h5,"obsm","dict","0.1.0")
    if "pca" in keys(obj.meta)
        h5["obsm/X_pca"] = obj.meta["pca"]' |> Matrix
        EnTypeVersion!(h5,"obsm/X_pca","array","0.2.0")
    end
    if "tsne" in keys(obj.meta)
        h5["obsm/X_tsne"] = obj.meta["tsne"]' |> Matrix
        EnTypeVersion!(h5,"obsm/X_tsne","array","0.2.0")
    end
    if "umap" in keys(obj.meta)
        h5["obsm/X_umap"] = obj.meta["umap"]' |> Matrix
        EnTypeVersion!(h5,"obsm/X_umap","array","0.2.0")
    end
    if "harmony" in keys(obj.meta)
        h5["obsm/X_harmony"] = obj.meta["harmony"]' |> Matrix
        EnTypeVersion!(h5,"obsm/X_harmony","array","0.2.0")
    end

    # obsp
    create_group(h5,"obsp")
    EnTypeVersion!(h5,"obsp","dict","0.1.0")

    # uns
    create_group(h5,"uns")
    EnTypeVersion!(h5,"uns","dict","0.1.0")
    for command in obj.log
        name = match(r"(.*?)\((.*)",command).captures[1]
        para = findall(r"([ a-zA-Zα-ωΑ-Ω_\-0-9]+=[ a-zA-Zα-ωΑ-Ω\"_\-0-9:\.]+)",
                       command) |> 
            x -> [ command[m] for m in x ] |> x -> split.(x,"=") |> 
            x -> [ strip.(e) for e in x ]
        # Might be some different 'DE'
        if occursin(r"^DE",name)
            group_name = para[findfirst(x -> x[1] == "group_name",para)][2]
            name = name * "." * group_name
        end
        create_group(h5,"uns/" * name)
        EnTypeVersion!(h5,"uns/" * name,"dict","0.1.0")
        create_group(h5,"uns/" * name * "/params")
        EnTypeVersion!(h5,"uns/" * name * "/params","dict","0.1.0")
        for p in para
            h5["uns/" * name * "/params/" * p[1]] = p[2]
            EnTypeVersion!(h5,"uns/" * name * "/params/" * p[1],"string",
                           "0.2.0")
        end
        if occursin(r"^PCA",name)
            h5["uns/" * name * "/pca_var"] = obj.meta["pca_var"]
            EnTypeVersion!(h5,"uns/" * name * "/pca_var","array","0.2.0")
            h5["uns/" * name * "/pca_cut"] = string(obj.meta["pca_cut"])
            EnTypeVersion!(h5,"uns/" * name * "/pca_cut","string",
                           "0.2.0")
        end
        if occursin(r"^SelectHVG",name)
            h5["uns/" * name * "/hvg_name"] = obj.meta["hvg_name"]
            EnTypeVersion!(h5,"uns/" * name * "/hvg_name","array","0.2.0")
            h5["uns/" * name * "/hvg_index"] = obj.meta["hvg_index"]
            EnTypeVersion!(h5,"uns/" * name * "/hvg_index","array","0.2.0")
        end
        if occursin(r"^DE",name)
            for n in names(obj.meta[group_name * "_DE"])
                h5["uns/" * name * "/" * n] = obj.meta[group_name * "_DE"][!,n]
                EnTypeVersion!(h5,"uns/" * name * "/" * n,"array","0.2.0")
            end
        end
    end

    # varp
    create_group(h5,"varp")
    EnTypeVersion!(h5,"varp","dict","0.1.0")

    close(h5)

    return "Finished!"
end

function EnTypeVersion!(h5::HDF5.File,
        parent::AbstractString,
        type::AbstractString,
        version::AbstractString)
    HDF5.attributes(h5[parent])["encoding-type"] = type
    HDF5.attributes(h5[parent])["encoding-version"] = version
end
