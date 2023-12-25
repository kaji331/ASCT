RuntimeGeneratedFunctions.init(@__MODULE__)

# Calculation the percentage of gene
function FeaturePercentage!(obj::WsObj;
        regex::Union{Regex,Nothing} = r"^[m|M][t|T]-",
        gene::Union{Vector{T} where T <: AbstractString,Nothing} = nothing,
        ptr::AbstractString = "raw",
        obs_name::AbstractString = "percentage_mt")

    # Check obs
    if !isnothing(regex) && !isnothing(gene)
        @warn "The 'regex' and 'gene' should not be defined simutaneously!" * 
            "Now just use the 'gene' parameter!"
    end
    if isnothing(regex) && isnothing(gene)
        @warn "Either 'regex' or 'gene' should be defined!"
        return "Nothing to do!"
    end

    mt = isnothing(gene) ? match.(regex,obj.var.name) |> 
        y -> findall(x -> !isnothing(x),y) : 
        findall(x -> x in gene,obj.var.name)

    if isempty(mt)
        return "Nothing to do! No features were selected!"
    end

    # Calculation
    # Fastest way
    d = obj.dat[ptr * "_dat"]' |> sparse
    percentage = d[:,mt] * ones(Int32,length(mt)) ./ (d * ones(Int32,d.n)) .* 
        100

    obj.obs[:,obs_name] = percentage

    # log
    push!(obj.log,String("FeaturePercentage!(WsObj;regex=$regex," * 
                         "gene=$gene,ptr=$ptr,obs_name=$obs_name"))

    return "Finished!"
end

# Filtering cells manually
function ManualFilter!(obj::WsObj;
        obs_name::Union{Nothing,AbstractString} = nothing,
        var_name::Union{Nothing,AbstractString} = nothing,
        expr::Union{AbstractString,Nothing} = nothing,
        ptr::AbstractString = "raw")

    if !isnothing(obs_name) && !isnothing(var_name)
        return "Nothing to do! Only need either 'obs_name' or 'var_name'!"
    end
    if isnothing(obs_name) && isnothing(var_name)
        return "Nothing to do! Need one of 'obs_name' and 'var_name'!"
    end

    if obs_name in names(obj.obs)
        if !isnothing(expr)
            f = Meta.parse("x -> x " * expr) |> 
                x -> @RuntimeGeneratedFunction(x)
            idx = findall(f,obj.obs[:,obs_name])
            # dat
            obj.dat[ptr * "_dat"] = obj.dat[ptr * "_dat"][:,idx]
            # obs
            obj.obs = obj.obs[idx,:]
        else
            return "Nothing to do! Must set filter expression!"
        end
    else
        if !isnothing(obs_name)
            return "'obs_name' is Not corrected!"
        end
    end

    if var_name in names(obj.var)
        if !isnothing(expr)
            f = Meta.parse("x -> x " * expr) |> 
                x -> @RuntimeGeneratedFunction(x)
            idx = findall(f,obj.var[:,var_name])
            # dat
            obj.dat[ptr * "_dat"] = obj.dat[ptr * "_dat"][idx,:]
            # var
            obj.var = obj.var[idx,:]
        else
            return "Nothing to do! Must set filter expression!"
        end
    else
        if !isnothing(var_name)
            return "'var_name' is Not corrected!"
        end
    end

    # log
    push!(obj.log,String("ManualFilter!(WsObj;obs_name=$obs_name," * 
                         "var_name=$var_name,expr=$expr,ptr=$ptr"))

    return "Finished!"
end

function AutoFilter!(obj::WsObj;
        obs_name::Union{Nothing,AbstractString} = "percentage_mt",
        operator::AbstractString = "<",
        filter_UmiFeature::Bool = true,
        ptr::AbstractString = "raw")

    # Filter cells by the elbow point of percentage_mt
    # ------
    if !(operator in ["<",">"])
        return "Nothing to do! 'operator' Needs \"<\" or \">\"!"
    end
    if obs_name in names(obj.obs)
        dat = obj.obs[:,obs_name]
        @info "Filtering cells by $obs_name ..."
        dense = KernelDensity.kde(dat)
        # Find peaks and check distribution
        peak_idx = findpeaks(dense.density,min_height=0.1)
        if length(peak_idx) > 1
            @warn "Distribution of data might be abnormal!?"
        end
        # Set start, end for linear regression and keep 'all' smaller/larger cells
        if operator == "<"
            s = maximum(peak_idx)
            e = length(dense.x)
        else
            s = 1
            e = minimum(peak_idx)
        end
        # Linear regression for one cutoff
        @info "Linear regression for threshold ..."
        dif = DataFrame(axis_x=[dense.x[s],dense.x[e]],
                        axis_y=[dense.density[s],dense.density[e]]) |> 
            x -> lm(@formula(axis_y ~ axis_x),x) |> 
            x -> GLM.predict(x,DataFrame(axis_x=dense.x[s:e])) |>
            x -> x - dense.density[s:e]
        cp = findmax(dif)[2] + s - 1
        cutoff = dense.x[cp] + 0.5 # 0.5 is an empiric value.
        @info "The cutoff is about $(round(cutoff;digits=2))"
        f = Meta.parse("x -> x " * operator * " $cutoff") |> 
            x -> @RuntimeGeneratedFunction(x)
        idx = findall(f,obj.obs[:,obs_name])
        # dat
        obj.dat[ptr * "_dat"] = obj.dat[ptr * "_dat"][:,idx]
        # obs
        obj.obs = obj.obs[idx,:]
    elseif !isnothing(obs_name)
        @warn "Nothing to do! 'obs_name' is Not existed!"
    end

    # Filter cells by the cluster of UMI/Feature
    # ------
    if filter_UmiFeature
        @info "Filtering cells by the cluster of UMI/Feature"
        Random.seed!(1984)
        dat = hcat(log10.(obj.obs.cell_counts),
                   log10.(obj.obs.cell_features))'
        clusters = dbscan(dat,5e-2)
        big_one = findmax(clusters.counts)[2]
        idx = clusters.clusters[big_one].core_indices
        # dat
        obj.dat[ptr * "_dat"] = obj.dat[ptr * "_dat"][:,idx]
        # obs
        obj.obs = obj.obs[idx,:]
    end

    # log
    push!(obj.log,String("AutoFilter!(WsObj;obs_name=$obs_name,ptr=$ptr)"))

    return "Finished!"
end

# Convenient plotting function for quality checking
function QcPlot(obj::WsObj;
        obs_name::Union{Nothing,AbstractString} = nothing,
        ptr::AbstractString = "raw")
    # Data checking
    if isnothing(obs_name)
        if "cell_counts" in names(obj.obs) && "cell_features" in names(obj.obs)
            x = obj.obs.cell_counts
            y = obj.obs.cell_features
            # combine 2 plots
            xticks = obj.obs.cell_counts |> maximum |> 
                x -> round(x / 10000) |> x -> [ y * 10000 for y in 1:x ]
            p = plot(x,y,st=[:scatter,:histogram2d],xticks=xticks,ma=0.1,ms=5,
                     markershape=:hexagon,msc="white",legend=false,layout=2,
                     bin=75)
            p.series_list[1].plotattributes[:markercolor] = "red"
            p.series_list[2].plotattributes[:fillcolor] = ColorSchemes.Reds_3[:]
            p
        else
            throw(DomainError(obs_name,
                              String("The 'cell_counts' and/or " * 
                                     "'cell_features' are/is not existed in " * 
                                     "the obs keys!")))
        end
    else
        if obs_name in names(obj.obs)
            x = obj.obs[:,obs_name]
            plot(x,st=:density,c="red",legend=false,linewidth=4,linealpha=0.5)
        else
            @warn "No $obs_name in the keys of 'obs'!"
        end
    end
end
