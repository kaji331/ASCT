function Read10X(Path::AbstractString;
        min_features::Union{Nothing,Integer} = nothing,
        min_cells::Union{Nothing,Integer} = nothing)
    if HDF5.ishdf5(Path)
        @info "Reading H5 file ..."
        try
            h5open(Path,"r") do h5
                key = keys(h5)[1]
                obj = GenerateObjH5(h5,key)

                if !(isnothing(min_features) && isnothing(min_cells))
                    SimpleFilter!(obj,min_features,min_cells)
                end

                # Make unique gene name
                dup = [ key for (key,val) in countmap(obj.var.name) 
                       if val != 1 ]
                if length(dup) != 0
                    @info "There are some duplicated gene names!" * 
                        " Add '.numbers' automatically!"
                    for d in dup
                        for (i,p) in enumerate(findall(x -> x == d,
                                                       obj.var.name))
                            obj.var.name[p] = obj.var.name[p] * 
                                (i > 1 ? ".$(i - 1)" : "")
                        end
                    end
                end

                return obj
            end
        catch
            throw(DomainError(Path,"Could not open the HDF5 file!"))
        end
    else
        if "matrix.mtx.gz" in readdir(Path) || "matrix.mtx" in readdir(Path)
            m = Ref(GenerateMatrix(Path))
            dat = m[][1]
            shape = m[][2]
        else 
            # Other quantification softwares, like kb-python, might generate a 
            # transposed matrix with different file names. We will implement the 
            # reading methods in other functions.
            throw(DomainError(Path,"No Matrix Market files!"))
        end
        if "barcodes.tsv.gz" in readdir(Path) || "barcodes.tsv" in readdir(Path)
            barcodes = GenerateBarcodes(Path)
        else 
            throw(DomainError(Path,"No Barcodes files!"))
        end
        if "features.tsv.gz" in readdir(Path) || 
            "genes.tsv" in readdir(Path)
            features = GenerateFeatures(Path)
        else 
            throw(DomainError(Path,"No Features/Genes files!"))
        end

        # Basic calculation
        @info "Gathering basic information"
        cell_counts = dat' * ones(Int32,shape[1])
        cell_features = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
        dat = dat' |> sparse
        cells = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
        features["feature_cells"] = cells
        not_zero = cells .!= 0
        for k in keys(features)
            features[k] = features[k][not_zero]
        end
        obj = WsObj(Dict("raw_dat" => dat[:,not_zero]' |> sparse),
                    DataFrame("barcodes" => barcodes,
                              "cell_counts" => cell_counts,
                              "cell_features" => cell_features),
                    DataFrame(features),
                    String[],
                    nothing)

        if !(isnothing(min_features) && isnothing(min_cells))
            SimpleFilter!(obj,min_features,min_cells)
        end

        # Make unique gene name
        dup = [ key for (key,val) in countmap(obj.var.name) if val != 1 ]
        if length(dup) != 0
            @info "There are some duplicated gene names!" * 
                " Add '.numbers' automatically!"
            for d in dup
                for (i,p) in enumerate(findall(x -> x == d,obj.var.name))
                    obj.var.name[p] = obj.var.name[p] * 
                        (i > 1 ? ".$(i - 1)" : "")
                end
            end
        end

        return obj
    end
end

function SimpleFilter!(obj::WsObj,
        min_features::Union{Nothing,Integer},
        min_cells::Union{Nothing,Integer})

    if isnothing(min_features)
        idx_col = Colon()
    else
        idx_col = obj.obs.cell_features .>= min_features
    end
    if isnothing(min_cells)
        idx_row = Colon()
    else
        idx_row = obj.var.feature_cells .>= min_cells
    end
    obj.dat["raw_dat"] = obj.dat["raw_dat"][idx_row,idx_col]
    obj.obs = @view obj.obs[idx_col,:]
    obj.var = @view obj.var[idx_row,:]
end

function GenerateObjH5(h5,key)
    dat = read(h5["$key/data"])
    indices = read(h5["$key/indices"])
    indptr = read(h5["$key/indptr"])
    shape = read(h5["$key/shape"])
    indices .+= 1
    indptr .+= 1

    # Column index, row index, and corresponding value
    #
    # append! will modify original object, might be better 
    # performance than vcat
    @info "Generating matrix ..."
    r = [ repeat([row],length(indptr[row]:indptr[row + 1] - 1)) 
         for row in 1:shape[2] ] |> x -> vcat(x...)
    dat = sparse(indices,r,dat,shape[1],shape[2])
    indices,indptr,r = nothing,nothing,nothing
    GC.gc()

    barcodes = read(h5["$key/barcodes"])
    if "features" in keys(h5["$key"])
        features = read(h5["$key/features"])
    elseif "genes" in keys(h5["$key"]) && "gene_names" in keys(h5["$key"])
        features = Dict{AbstractString,
                        AbstractVector}("id" => read(h5["$key/genes"]),
                                        "name" => read(h5["$key/gene_names"]))
    end
    chemistry = read_attribute(h5,"chemistry_description")

    # Basic calculation
    @info "Gathering basic information"
    cell_counts = dat' * ones(Int32,shape[1])
    cell_features = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
    dat = dat' |> sparse
    cells = [ dat.colptr[col + 1] - dat.colptr[col] for col in 1:dat.n ]
    features["feature_cells"] = cells
    not_zero = cells .!= 0
    for k in keys(features)
        if length(features[k]) == length(cells)
            features[k] = features[k][not_zero]
        else
            features[k] = repeat(features[k],count(not_zero))
        end
    end

    return WsObj(Dict("raw_dat" => dat[:,not_zero]' |> sparse),
                 DataFrame("barcodes" => barcodes,
                           "cell_counts" => cell_counts,
                           "cell_features" => cell_features),
                 DataFrame(features),
                 String[],
                 Dict("chemistry_description" => chemistry))
end

function GenerateMatrix(Path)
    @info "Currently, only support Cellranger's result!!!"
    @info "Reading matrix directory ..."
    shape = Int32[]
    c = Int32[]
    r = Int32[]
    v = Int32[]
    @info "Generating matrix ..."
    try
        GZip.open(Path * "/matrix.mtx.gz","r") do f
            counter = 1
            for line in eachline(f)
                if startswith(line,"%")
                    continue
                else
                    if counter == 1
                        counter = -1
                        append!(shape,parse.(Int32,split(line)))
                    else
                        ele = split(line)
                        push!(c,parse(Int32,ele[1]))
                        push!(r,parse(Int32,ele[2]))
                        push!(v,parse(Int32,ele[3]))
                    end
                end
            end
        end
    catch
        open(Path * "/matrix.mtx","r") do f
            counter = 1
            for line in eachline(f)
                if startswith(line,"%")
                    continue
                else
                    if counter == 1
                        counter = -1
                        append!(shape,parse.(Int32,split(line)))
                    else
                        ele = split(line)
                        push!(c,parse(Int32,ele[1]))
                        push!(r,parse(Int32,ele[2]))
                        push!(v,parse(Int32,ele[3]))
                    end
                end
            end
        end
    end

    dat = sparse(c,r,v,shape[1],shape[2])
    return (dat,shape)
end

function GenerateBarcodes(Path)
    barcodes = String[]
    @info "Reading barcodes ..."
    try
        GZip.open(Path * "/barcodes.tsv.gz","r") do f
            append!(barcodes,readlines(f))
        end
    catch
        open(Path * "/barcodes.tsv","r") do f
            append!(barcodes,readlines(f))
        end
    end
    
    return barcodes
end

function GenerateFeatures(Path)
    id = String[]
    name = String[]
    ft = String[]
    @info "Reading features/genes ..."
    try
        GZip.open(Path * "/features.tsv.gz","r") do f
            for line in eachline(f)
                l = split(line,"\t")
                push!(id,l[1])
                push!(name,l[2])
                push!(ft,l[3])
            end
        end
    catch
        open(Path * "/genes.tsv","r") do f
            for line in eachline(f)
                l = split(line,"\t")
                push!(id,l[1])
                push!(name,l[2])
            end
        end
    end
    if isempty(ft)
        features = Dict{AbstractString,AbstractVector}("id" => id,
                                                       "name" => name)
    else
        features = Dict{AbstractString,
                        AbstractVector}("id" => id,
                                        "name" => name,
                                        "feature_type" => ft)
    end

    return features
end
