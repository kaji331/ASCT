"""
    FeatureHeat(obj)

Draw a heatmap of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `features::Union{Integer,Vector{<: AbstractString}} = 10`: top DE number of 
  each cell group or the vector of feature names.
- `colorscheme::Symbol = :Hiroshige`: cells colorscheme from ColorSchemes 
  package.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `group_colorscheme::Symbol = :tol_rainbow`: groups colorscheme from 
  ColorSchemes package.
- `gene_label_size::Union{Nothing,Integer} = nothing`: set the size of gene 
  labels or switch off gene labels.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 30`: size of title font.
- `width::Integer = 600`: figure width.
- `height::Integer = 300`: figure height.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureHeat(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 10,
        colorscheme::Symbol = :Hiroshige,
        group_name::AbstractString = "clusters_latest",
        group_colorscheme::Symbol = :tol_rainbow,
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 30,
        width::Integer = 640,
        height::Integer = 480,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    grid_layout = figure[1,1] = GridLayout()

    # Data preparation
    if typeof(features) <: Integer
        gene_id = sort(obj.meta[group_name * "_DE"],[:group,:log2fc];
                       rev=[false,true],lt=natural) |> 
            x -> groupby(x,:group) |>
            x -> DataFrames.combine(x,:gene => y -> first(y,features)) |> 
            x -> indexin(unique(x.gene_function),obj.var.name)
    else
        gene_id = indexin(unique(features),obj.var.name)
    end
    levels = unique(obj.obs[!,group_name]) |> x -> sort(x;lt=natural)
    cell_id = vcat([ findall(x -> x == level,obj.obs[!,group_name]) 
                    for level in levels ]...)
    groups = vcat([ repeat([level],count(obj.obs[!,group_name] .== level)) 
                   for level in levels ]...)
    data = obj.dat["scale_dat"][gene_id,cell_id]

    # Group bar
    axis1,heat1 = heatmap(grid_layout[1,1],groups' |> Matrix |> transpose,
                          colormap=colorschemes[group_colorscheme])
    hidedecorations!(axis1,grid=false)

    # Heatmap
    cs = colorschemes[colorscheme]
    idx = range(1,length(cs);length=8) |> x -> round.(Integer,x) |> reverse
    cs = cs[idx]
    axis2,heat2 = heatmap(grid_layout[2,1],reverse(data',dims=2),
                          colormap=cs)

    # Legends
    Colorbar(grid_layout[2,2],heat2,height=round(Integer,height / 3 * 2))
    axis3,heat3 = heatmap(grid_layout[2,3],reverse(levels)' |> Matrix,
                          colormap=colorschemes[group_colorscheme])

    # Details
    axis1.height = 20
    axis1.width = width

    axis2.width = width
    axis2.height = height
    if isnothing(gene_label_size)
        axis2.yticklabelsvisible = false
    else
        axis2.yticks = (eachindex(gene_id) |> reverse,obj.var.name[gene_id])
        axis2.yticklabelsize = gene_label_size
        axis2.yticklabelrotation = π / 10
    end
    axis2.yticksvisible = false
    axis2.xticklabelsvisible = false
    axis2.xticksvisible = false

    axis3.height = round(Integer,height / 3 * 2)
    axis3.width = 20
    axis3.yticks = (eachindex(levels) |> reverse,string.(levels))
    axis3.yaxisposition = :right
    axis3.yticksvisible = false
    axis3.xticklabelsvisible = false
    axis3.xticksvisible = false
    rowgap!(grid_layout,5)
    colgap!(grid_layout,20)
    Label(grid_layout[0,:],text=title,fontsize=title_size,font=:bold)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    DimensionPoints(obj)

Draw a scatter figure of cells.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `dimension_name::AbstractString = "umap"`: name of dimension data in meta.
- `group_name::Union{Nothing,AbstractString} = "clusters_latest"`: a column name 
  in obs for different colors.
- `split_name::Union{Nothing,AbstractString} = nothing`: a column name in obs 
  for different facets.
- `colorscheme::Symbol = :tableau_20`: group colors from ColorSchemes package.
- `dim_prefix::Union{Nothing,AbstractString} = nothing`: prefix of axis labels. 
  default is the dimension_name.
- `point_size:Integer = 7`: change the size of points.
- `alpha::AbstractFloat = 0.7`: transparence of points from 0 to 1.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 30`: size of title font.
- `width::Integer = 400`: figure width.
- `height::Integer = 300`: figure height.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function DimensionPoints(obj::WsObj;
        dimension_name::AbstractString = "umap",
        group_name::Union{Nothing,AbstractString} = "clusters_latest",
        split_name::Union{Nothing,AbstractString} = nothing,
        colorscheme::Symbol = :tableau_20,
        dim_prefix::Union{Nothing,AbstractString} = nothing,
        point_size::Integer = 5,
        alpha::AbstractFloat = 0.5,
        title::AbstractString = "",
        title_size::Integer = 30,
        width::Integer = 640,
        height::Integer = 480,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    grid_layout_main = figure[1,1] = GridLayout()
    grid_layout_legend = figure[1,2] = GridLayout()

    # Data preparation
    if !isnothing(group_name)
        groups = obj.obs[!,group_name]
        levels = unique(groups) |> x -> sort(x;lt=natural)
    else
        groups = 1
        levels = 1
    end
    if !isnothing(split_name)
        facets = obj.obs[!,split_name]
        facet_levels = unique(facets) |> x -> sort(x;lt=natural)
        facet_rows = facet_cols = length(facet_levels) |> sqrt |> ceil |> Int
        if length(facet_levels) / facet_cols < facet_rows
            facet_rows = mod(length(facet_levels),facet_cols) == 0 ? 
                facet_rows - 1 : 
                facet_rows
        end
    else
        facets = 1
        facet_levels = 1
        facet_rows = facet_cols = 1
    end
    data = obj.meta[dimension_name][:,1:2]
    dim_prefix = isnothing(dim_prefix) ? 
        uppercase(dimension_name) : dim_prefix
    xrange = [minimum(data[:,1]),maximum(data[:,1])]
    xmargin = (xrange[2] - xrange[1]) / 10
    xlim = [ i == 1 ? v - xmargin : v + xmargin for (i,v) in enumerate(xrange) ]
    yrange = [minimum(data[:,2]),maximum(data[:,2])]
    ymargin = (yrange[2] - yrange[1]) / 10
    ylim = [ i == 1 ? v - ymargin : v + ymargin for (i,v) in enumerate(yrange) ]
    
    # Colorscheme
    if colorscheme in keys(colorschemes)
        cs = colorschemes[colorscheme]
        if 1 < length(levels) <= length(cs)
            cs = cs[round.(Integer,range(1,length(cs);length=length(levels)))] 
        else
            @warn "Only 1 group or groups number is more than colors number!"
            cs = cs[round.(Integer,range(1,length(cs);length=length(levels)))] 
        end
    else
        cs = cgrad(colorscheme,length(levels);categorical=true)
    end

    # Plots
    counter = 1
    for row in 1:facet_rows
        for col in 1:facet_cols
            # Facets
            if counter > length(facet_levels)
                break
            end

            # Axis
            axis = Axis(grid_layout_main[row,col])
            axis.width = width / (facet_cols + 0.3)
            axis.height = height / (facet_rows + 0.3)
            axis.limits = (xlim,ylim)
            if facet_levels == 1
                hidespines!(axis,:t,:r,:l,:b)
            end
            hidedecorations!(axis,grid=false)

            # Scatters
            for (i,level) in enumerate(levels)
                idx = (groups .== level) .& (facets .== facet_levels[counter])
                scatter!(axis,data[:,1][idx],data[:,2][idx],alpha=alpha,
                         label=string(level),markersize=point_size,
                         color=cs[i])
            end

            # Arrows
            if row == facet_rows && col == 1 && length(facet_levels) == 1
                end_len1 = (range(xlim[1],xlim[2],length=10)[2] - xlim[1]) * 
                    facet_cols
                end_len2 = (range(ylim[1],ylim[2],length=10)[2] - ylim[1]) * 
                    facet_rows
                factor = 0.2
                arrows!(axis,
                        [xlim[1] + factor * xmargin],
                        [ylim[1] + factor * ymargin],
                        [end_len1],
                        [0])
                arrows!(axis,
                       [xlim[1] + factor * xmargin],
                       [ylim[1] + factor * ymargin],
                       [0],
                       [end_len2])
            end

            counter += 1
        end
    end
    Label(grid_layout_main[0,:],text=title,fontsize=title_size,font=:bold)
    if length(facet_levels) == 1
        Label(grid_layout_main[end + 1,:],text=dim_prefix * "_1",
              fontsize=title_size / 2,halign=:left)
        Label(grid_layout_main[end - 1,0],text=dim_prefix * "_2",
              fontsize=title_size / 2,rotation=1.6,valign=:bottom)
    end
    # Make the fixed size of legend elements 
    legend_elements = [ MarkerElement(color=c,marker=:circle,makersize=15) 
                       for c in cs ]
    # legend_elements = [ MarkerElement(color=c,marker='●',makersize=15) 
    #                    for c in figure.scene.theme.palette.color[] ]
    Legend(grid_layout_legend[1,1],legend_elements,
           [ string(level) for level in levels ])
    colgap!(grid_layout_main,10)
    rowgap!(grid_layout_main,10)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    DrawQC(obj)

Draw quality checking figures.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `obs_name::Union{Nothing,AbstractString} = nothing: draw the distribution of 
  a column in obs.
- `scale_fun::Function = log10`: a scale function for counts/features figure.
- `hex_colorscheme::Symbol = :thermal`: heat color map of counts/features from 
  ColorSchemes package.
- `hex_number::Union{Integer,Tuple{Integer,Integer}} = 75`: number of hex bins 
  for x and y direction.
- `density_color::Tuple{Symbol,AbstractFloat} = (:orangered,0.3)`: color for 
  density figure of obs column.
- `height::Integer = 300`: figure height.
- `width::Union{Symbol,Integer} = :auto`: figure width.
- `title::Union{Symbol,AbstractString} = :auto`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function DrawQC(obj::WsObj;
        obs_name::Union{Nothing,AbstractString} = nothing,
        scale_fun::Function = log10,
        hex_colorscheme::Symbol = :viridis,
        hex_number::Union{Integer,Tuple{Integer,Integer}} = 75,
        density_color::Tuple{Symbol,AbstractFloat} = (:orangered,0.3),
        height::Integer = 480,
        width::Union{Symbol,Integer} = :auto,
        title::Union{Symbol,AbstractString} = :auto,
        title_size::Integer = 15,
        path::Union{Nothing,AbstractString} = nothing)
    
    # Data
    if isnothing(obs_name)
        if "cell_counts" in names(obj.obs) && "cell_features" in names(obj.obs)
            @info "Drawing counts/features hex plot..."
            x = obj.obs.cell_counts
            y = obj.obs.cell_features

            figure = Figure(fonts=(;regular="Nimbus Roman",
                                   bold="Nimbus Roman Bold"))
            if title == :auto
                title = "Distribution of Counts/Features"
            end
            axis = Axis(figure[1,1],xlabel="Cell Counts",ylabel="Cell Features",
                        title=title,titlefont=:bold,xtickformat="{:d}",
                        ytickformat="{:d}",titlesize=title_size,
                        xticklabelrotation=π / 4)
            axis.limits = ([0,maximum(x) + maximum(x) * 0.1],
                           [0,maximum(y) + maximum(y) * 0.1])
            xrange = maximum(x) - minimum(x)
            yrange = maximum(y) - minimum(y)
            ratio = yrange / xrange
            if width == :auto
                xcenter = height * ratio
                xleft = height * minimum(x) / xrange
                xright = height * maximum(x) / xrange * 0.1
                axis.width,axis.height = (xcenter + xleft + xright) * 1.5,height
            else
                axis.width,axis.height = width,height
            end
            hexbin!(axis,x,y;bins=hex_number,colorscale=scale_fun,
                    colormap=hex_colorscheme)
        else
            @warn "The 'cell_counts' and/or 'cell_features' are/is not " * 
                "existed in the obs keys!"
            return "Nothing to do!"
        end
    else
        if obs_name in names(obj.obs)
            x = obj.obs[!,obs_name]

            figure = Figure(fonts=(;regular="Nimbus Roman",
                                   bold="Nimbus Roman Bold"))
            if title == :auto
                title = "Distribution of '$obs_name'"
            end
            axis = Axis(figure[1,1],title=title,titlefont=:bold,
                        width=width == :auto ? height * 1.3 : width,
                        height=height)
            density!(axis,x;color=density_color,strokecolor=density_color[1],
                     strokewidth=3)
        else
            @warn "No $obs_name in the keys of 'obs'!"
            return "Nothing to do!"
        end
    end

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    FeatureVariances(obj)

Draw the variances of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `title::Union{Symbol,AbstractString} = :auto`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 400`: figure width.
- `height::Integer = 300`: figure height.
- `point_size:Integer = 5`: change the size of points.
- `alpha::AbstractFloat = 0.5`: transparence of points from 0 to 1.
- `hvg_color:Symbol = :orangered`: color of highly variable features.
- `base_color:Symbol = :gray`: color of base points.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureVariances(obj::WsObj;
        title::Union{Symbol,AbstractString} = :auto,
        title_size::Integer = 20,
        width::Integer = 400,
        height::Integer = 300,
        point_size::Integer = 5,
        alpha::AbstractFloat = 0.5,
        hvg_color::Symbol = :orangered,
        base_color::Symbol = :gray,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    if title == :auto
        title = "Standard variances of features"
    end
    axis = Axis(figure[1,1],title=title,titlesize=title_size,titlefont=:bold,
                xlabel="Means",ylabel="Standard variances",width=width,
                height=height)
    index_red = obj.meta["hvg_index"]
    index_black = eachindex(obj.var.name) |> 
        x -> x .∉ Ref(intersect(x,index_red))
    scatter!(axis,obj.var.hvg_mean[index_black],
             obj.var.hvg_var_std[index_black];
             color=base_color,alpha=alpha,size=point_size)
    scatter!(axis,obj.var.hvg_mean[index_red],
             obj.var.hvg_var_std[index_red];
             color=hvg_color,alpha=alpha,size=point_size)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    ElbowPCA(obj)

Draw the elbow scatter-line of PC variances.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `title::Union{Symbol,AbstractString} = :auto`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 400`: figure width.
- `height::Integer = 300`: figure height.
- `point_color::Union{Symbol,Tuple{Symbol,AbstractFloat}} = :gray90`: change 
  the filled color of points.
- `point_size::Integer = 12`: change the size of points.
- `line_color::Symbol = :black`: color of line.
- `line_width::Real = 0.8`: width of line.
- `stroke_color::Symbol = :gray`: color of point edges.
- `stroke_width::Integer = 1`: width of point edges.
- `threshold_line::Bool = true`: add the threshold line of PCs cut.
- `threshold_line_color::Symbol = :gray`: the color of threshold line.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function ElbowPCA(obj::WsObj;
        title::Union{Symbol,AbstractString} = :auto,
        title_size::Integer = 20,
        width::Integer = 400,
        height::Integer = 300,
        point_color::Union{Symbol,Tuple{Symbol,AbstractFloat}} = :gray90,
        point_size::Integer = 12,
        line_color::Symbol = :black,
        line_width::Real = 0.8,
        stroke_color::Symbol = :gray,
        stroke_width::Integer = 1,
        threshold_line::Bool = true,
        threshold_line_color::Symbol = :gray,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    if title == :auto
        title = "Elbow figure of PCA"
    end
    axis = Axis(figure[1,1],title=title,titlesize=title_size,titlefont=:bold,
                xlabel="PCs",ylabel="Proportions of variances",width=width,
                height=height)
    scatterlines!(axis,1:size(obj.meta["pca"],2),obj.meta["pca_var"];
                  markersize=point_size,color=line_color,markercolor=point_color,
                  strokecolor=stroke_color,strokewidth=stroke_width,
                  linewidth=line_width)
    if threshold_line
        vlines!(axis,obj.meta["pca_cut"];ymin=0.01,ymax=0.95,linestyle=:dash,
                linewidth=2,color=threshold_line_color)
    end

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    FeatureViolin(obj)

Draw a violin figure of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `features::Union{Integer,Vector{<: AbstractString}} = 2`: top DE number of 
  each cell group or the vector of feature names.
- `colorscheme::Symbol = :Hiroshige`: group colorscheme from ColorSchemes 
  package.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `gene_label_size::Union{Nothing,Integer} = nothing`: set the size of gene 
  labels or switch off gene labels.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 400`: figure width.
- `height::Integer = 300`: figure height.
- `font_size::Integer = 10`: font size factor.
- `y_axis_pos::Symbol = :right`: set the position of gene labels.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureViolin(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        colorscheme::Symbol = :tableau_20,
        group_name::AbstractString = "clusters_latest",
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Integer = 580,
        height::Integer = 480,
        font_size::Integer = 10,
        y_axis_pos::Symbol = :right,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))

    # Data preparation
    if typeof(features) <: Integer
        gene_id = sort(obj.meta[group_name * "_DE"],[:group,:log2fc];
                       rev=[false,true],lt=natural) |> 
            x -> groupby(x,:group) |> 
            x -> DataFrames.combine(x,:gene => y -> first(y,features)) |> 
            x -> indexin(unique(x.gene_function),obj.var.name)
    else
        gene_id = indexin(unique(features),obj.var.name)
    end
    levels = unique(obj.obs[!,group_name]) |> x -> sort(x;lt=natural)
    cell_id = vcat([ findall(x -> x == level,obj.obs[!,group_name]) 
                    for level in levels ]...)
    groups = vcat([ repeat([level],count(obj.obs[!,group_name] .== level)) 
                   for level in levels ]...)
    data = obj.dat["norm_dat"][gene_id,cell_id]

    # Colorscheme
    if colorscheme in keys(colorschemes)
        cs = colorschemes[colorscheme]
        cs = length(levels) <= length(cs) ? 
            cs[round.(Integer,range(1,length(cs);length=length(levels)))] : 
            nothing
    else
        cs = repeat([colorscheme],length(levels))
    end

    ll = length(levels)
    row_height = height / length(gene_id)
    # Violin for each gene
    for gid_ind in eachindex(gene_id)
        axis = Axis(figure[gid_ind,1],xlabel=group_name,
                    xticks=(eachindex(levels),string.(levels)),
                    height=row_height,width=width,ytickformat="{:2.1f}")
        # Violin for each group
        for (i,level) in enumerate(levels)
            idx = groups .== level
            value = data[gid_ind,idx] |> Vector
            if isnothing(cs)
                violin!(axis,repeat([i],count(idx)),value;
                        datalimits=(0,Inf))
            else
                violin!(axis,repeat([i],count(idx)),value;
                        datalimits=(0,Inf),color=cs[i])
            end
        end
        axis.yticks = [maximum(data[gid_ind,cell_id])]
        axis.yticklabelsize = 10
        axis.yticksize = font_size * 0.3
        axis.spinewidth = 0.2
        axis.yaxisposition = y_axis_pos
        axis.ylabel = obj.var.name[gene_id[gid_ind]]
        axis.ylabelsize = isnothing(gene_label_size) ? 
            font_size * 1.2 : 
            gene_label_size
        axis.ylabelrotation = π * 2
        if gid_ind != length(gene_id)
            hidexdecorations!(axis)
        else
            hidexdecorations!(axis,ticklabels=false,label=false)
            axis.xticklabelrotation = π / 4
            axis.xticklabelsize = font_size * 1.2
        end
        hideydecorations!(axis,ticks=false,ticklabels=false,label=false)
        hidespines!(axis,:t,:r,:l)
    end
    rowgap!(figure.layout,0)
    Label(figure.layout[0,:],text=title,fontsize=title_size,font=:bold)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    FeatureFracDots(obj)

Draw a bubble figure of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `features::Union{Integer,Vector{<: AbstractString}} = 2`: top DE number of 
  each cell group or the vector of feature names.
- `colorscheme::Union{Symbol,Vector{Symbol}} = [:gray90,:orangered,:red]`: 
  colorscheme for expression levels.
- `color_length::Integer = 100`: gradient length of colorscheme.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `gene_label_size::Union{Nothing,Integer} = nothing`: set the size of gene 
  labels or switch off gene labels.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 580`: figure width.
- `height::Integer = 480`: figure height.
- `font_size::Integer = 10`: font size factor.
- `circle_size::Integer = 30`: circle size factor.
- `y_axis_pos::Symbol = :right`: set the position of gene labels.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureFracDots(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        colorscheme::Union{Symbol,Vector{Symbol}} = [:gray90,:orangered,:red],
        color_length::Integer = 100,
        group_name::AbstractString = "clusters_latest",
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Integer = 580,
        height::Integer = 480,
        font_size::Integer = 10,
        circle_size::Integer = 30,
        y_axis_pos::Symbol = :right,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))

    # Data preparation
    if typeof(features) <: Integer
        gene_id = sort(obj.meta[group_name * "_DE"],[:group,:log2fc];
                       rev=[false,true],lt=natural) |> 
            x -> groupby(x,:group) |> 
            x -> DataFrames.combine(x,:gene => y -> first(y,features)) |> 
            x -> indexin(unique(x.gene_function),obj.var.name)
    else
        gene_id = indexin(unique(features),obj.var.name)
    end
    levels = unique(obj.obs[!,group_name]) |> x -> sort(x;lt=natural)
    cell_id = vcat([ findall(x -> x == level,obj.obs[!,group_name]) 
                    for level in levels ]...)
    groups = vcat([ repeat([level],count(obj.obs[!,group_name] .== level)) 
                   for level in levels ]...)
    data = obj.dat["norm_dat"][gene_id,cell_id]

    # Colorscheme
    if typeof(colorscheme) <: Symbol
        cs = colorschemes[colorscheme]
    else
        cs = cgrad(colorscheme,color_length)
    end
    color_start = range(minimum(data);stop=maximum(data),length=length(cs))

    ll = length(levels)
    lg = length(gene_id)
    row_height = height / length(gene_id)
    # Dots for each gene
    for gid_ind in eachindex(gene_id)
        axis = Axis(figure[gid_ind,1],xlabel=group_name,
                    xticks=(eachindex(levels),string.(levels)),
                    height=row_height,width=width)
        # Dots for each group
        for (i,level) in enumerate(levels)
            idx = groups .== level
            value = data[gid_ind,idx]
            value_mean = mean(value)
            percentage = count(value .!= 0) / count(idx)
            color = cs[findlast(x -> value_mean >= x,color_start)]
            scatter!(axis,i,5;datalimits=(0,Inf),color=color,strokewidth=0.5,
                     strokecolor=:black,
                     markersize=(percentage + 0.05) * circle_size)
        end
        axis.yaxisposition = y_axis_pos
        axis.ylabel = obj.var.name[gene_id[gid_ind]]
        axis.ylabelsize = isnothing(gene_label_size) ? 
            font_size * 1.2 : 
            gene_label_size
        axis.ylabelrotation = π * 2
        if gid_ind != lg
            hidexdecorations!(axis)
        else
            hidexdecorations!(axis,ticklabels=false,label=false)
            axis.xticklabelrotation = π / 4
            axis.xticklabelsize = font_size * 1.2
        end
        hideydecorations!(axis,label=false)
        hidespines!(axis,:t,:b,:l,:r)
    end
    rowgap!(figure.layout,0)
    Label(figure.layout[0,:],text=title,fontsize=title_size,font=:bold)
    Colorbar(figure.layout[round(Integer,lg * 0.45):round(Integer,lg * 0.55),
                           end + 1],
             limits=(0,maximum(data)),colormap=cs,tickformat="{:2.1f}",
             ticks=[0,maximum(data)],ticklabelsize=font_size * 1.2)
    legend_elements = [ MarkerElement(color=:grey80,strokewidth=0.5,
                                      strokecolor=:black,marker=:circle,
                                      markersize=size)
                       for size in 
                       range(0;stop=circle_size + 0.05,length=6)[2:6]]
    Legend(figure.layout[round(Integer,lg * 0.45):round(Integer,lg * 0.55),
                         end + 1],
           legend_elements,["0.2","0.4","0.6","0.8","1.0"],labelsize=12)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    FeatureDimension(obj)

Draw a dimension scatter figure of features expression.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `dimension_name::AbstractString = "umap"`: name of dimension data in meta.
- `features::Union{Integer,Vector{<: AbstractString}} = 2`: an integer or a 
  gene name vector. The integer means top number of each group for DE result.
- `de_name::Union{Nothing,AbstractString} = "clusters_latest_DE"`: the DE name 
  in the meta.
- `color::Union{Symbol,Vector{Symbol}} = 
  [:steelblue,:gray90,:orange,:orangered,:red]`: a symbol of colormap or a 
  vector of color symbols.
- `color_length::Integer = 100`: the number of color gradient.
- `dim_prefix::Union{Nothing,AbstractString} = nothing`: prefix of axis labels. 
  default is the dimension_name.
- `point_size:Integer = 5`: change the size of points.
- `alpha::AbstractFloat = 0.5`: transparence of points from 0 to 1.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 30`: size of title font.
- `width::Integer = 640`: figure width.
- `height::Integer = 480`: figure height.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureDimension(obj::WsObj;
        dimension_name::AbstractString = "umap",
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        de_name::Union{Nothing,AbstractString} = "clusters_latest_DE",
        color::Union{Symbol,Vector{Symbol}} = 
            [:steelblue,:gray90,:orange,:orangered,:red],
        color_length::Integer = 100,
        dim_prefix::Union{Nothing,AbstractString} = nothing,
        point_size::Integer = 5,
        alpha::AbstractFloat = 0.5,
        title::AbstractString = "",
        title_size::Integer = 30,
        width::Integer = 640,
        height::Integer = 480,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    grid_layout_main = figure[1,1] = GridLayout()
    grid_layout_legend = figure[1,2] = GridLayout()

    # Data preparation
    if typeof(features) <: Integer && typeof(de_name) <: AbstractString
        features = groupby(obj.meta[de_name],:group) |> 
            x -> DataFrames.combine(x) do y
                first(y,features)
            end |> x -> x.gene |> unique
    end
    if length(features) > 1
        facet_rows = facet_cols = length(features) |> sqrt |> ceil |> Int
        if length(features) / facet_cols < facet_rows
            facet_rows = mod(length(features),facet_cols) == 0 ? 
                facet_rows - 1 : 
                facet_rows
        end
    else
        facet_rows = facet_cols = 1
    end
    data = obj.meta[dimension_name][:,1:2]
    value = obj.dat["norm_dat"][indexin(features,obj.var.name),:]
    dim_prefix = isnothing(dim_prefix) ? 
        uppercase(dimension_name) : dim_prefix
    xrange = [minimum(data[:,1]),maximum(data[:,1])]
    xmargin = (xrange[2] - xrange[1]) / 10
    xlim = [ i == 1 ? v - xmargin : v + xmargin for (i,v) in enumerate(xrange) ]
    yrange = [minimum(data[:,2]),maximum(data[:,2])]
    ymargin = (yrange[2] - yrange[1]) / 10
    ylim = [ i == 1 ? v - ymargin : v + ymargin for (i,v) in enumerate(yrange) ]
    
    # Colorscheme
    if typeof(color) <: Symbol
        if color in keys(colorschemes)
            cs = cgrad(color,color_length)
        else
            return "Nothing to do! 'color' is not an existed colorscheme!"
        end
    else
        cs = cgrad(color,color_length)
    end
    value_windows = range(minimum(obj.dat["norm_dat"]),
                          maximum(obj.dat["norm_dat"]);length=color_length)

    # Plots
    counter = 1
    for row in 1:facet_rows
        for col in 1:facet_cols
            # Facets
            if counter > length(features)
                break
            end

            # Axis
            axis = Axis(grid_layout_main[row,col];title=features[counter],
                        titlesize=16 / sqrt(sqrt(length(features))))
            axis.width = width / (facet_cols + 0.3)
            axis.height = height / (facet_rows + 0.3)
            axis.limits = (xlim,ylim)
            if length(features) == 1
                hidespines!(axis,:t,:r,:l,:b)
            end
            hidedecorations!(axis,grid=false)

            # Scatters
            for cell_id in 1:value.n
                color_id = findfirst(x -> x >= value[counter,cell_id],
                                     value_windows)
                scatter!(axis,data[:,1][cell_id],data[:,2][cell_id],alpha=alpha,
                         markersize=point_size,color=cs[color_id])
            end

            # Arrows
            if row == facet_rows && col == 1 && length(features) == 1
                end_len1 = (range(xlim[1],xlim[2],length=10)[2] - xlim[1]) * 
                    facet_cols
                end_len2 = (range(ylim[1],ylim[2],length=10)[2] - ylim[1]) * 
                    facet_rows
                factor = 0.2
                arrows!(axis,
                        [xlim[1] + factor * xmargin],
                        [ylim[1] + factor * ymargin],
                        [end_len1],
                        [0])
                arrows!(axis,
                       [xlim[1] + factor * xmargin],
                       [ylim[1] + factor * ymargin],
                       [0],
                       [end_len2])
            end

            counter += 1
        end
    end
    Label(grid_layout_main[0,:],text=title,fontsize=title_size,font=:bold)
    if length(features) == 1
        Label(grid_layout_main[end + 1,:],text=dim_prefix * "_1",
              fontsize=title_size / 2,halign=:left)
        Label(grid_layout_main[end - 1,0],text=dim_prefix * "_2",
              fontsize=title_size / 2,rotation=1.6,valign=:bottom)
    end
    # Make the fixed size of legend elements 
    Colorbar(grid_layout_legend[1,1],colormap=cs,
             limits=(value_windows[1],value_windows[end]),
             ticks=[value_windows[1],value_windows[end] / 2,value_windows[end]],
             tickformat="{:2.1f}",height=height / 3,tellheight=false)
    colgap!(grid_layout_main,10)
    rowgap!(grid_layout_main,10)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    FeatureJitters(obj)

Draw a jitter figure of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `features::Union{Integer,Vector{<: AbstractString}} = 2`: top DE number of 
  each cell group or the vector of feature names.
- `colorscheme::Symbol = :tableau_20`: group colorscheme from ColorSchemes 
  package.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `gene_label_size::Union{Nothing,Integer} = nothing`: set the size of gene 
  labels or switch off gene labels.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 640`: figure width.
- `height::Integer = 480`: figure height.
- `font_size::Integer = 10`: font size factor.
- `point_size::Integer = 5`: point size.
- `alpha::Integer = 0.5`: point alpha.
- `jitter_width::Real = 0.75`: width among jitter groups.
- `y_axis_pos::Symbol = :right`: set the position of gene labels.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureJitters(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        colorscheme::Symbol = :tableau_20,
        group_name::AbstractString = "clusters_latest",
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Integer = 640,
        height::Integer = 480,
        font_size::Integer = 10,
        point_size::Integer = 5,
        alpha::AbstractFloat = 0.5,
        jitter_width::Real = 0.75,
        y_axis_pos::Symbol = :right,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))

    # Data preparation
    if typeof(features) <: Integer
        gene_id = sort(obj.meta[group_name * "_DE"],[:group,:log2fc];
                       rev=[false,true],lt=natural) |> 
            x -> groupby(x,:group) |> 
            x -> DataFrames.combine(x,:gene => y -> first(y,features)) |> 
            x -> indexin(unique(x.gene_function),obj.var.name)
    else
        gene_id = indexin(unique(features),obj.var.name)
    end
    levels = unique(obj.obs[!,group_name]) |> x -> sort(x;lt=natural)
    cell_id = vcat([ findall(x -> x == level,obj.obs[!,group_name]) 
                    for level in levels ]...)
    groups = vcat([ repeat([level],count(obj.obs[!,group_name] .== level)) 
                   for level in levels ]...)
    data = obj.dat["norm_dat"][gene_id,cell_id]

    # Colorscheme
    if colorscheme in keys(colorschemes)
        cs = colorschemes[colorscheme]
        cs = length(levels) <= length(cs) ? 
            cs[round.(Integer,range(1,length(cs);length=length(levels)))] : 
            nothing
        if !isnothing(cs)
            cs = cs[indexin(groups,levels)]
        end
    else
        cs = colorscheme
    end

    ll = length(levels)
    row_height = height / length(gene_id)
    # RainClouds for each gene
    for gid_ind in eachindex(gene_id)
        axis = Axis(figure[gid_ind,1],xlabel=group_name,
                    xticks=(eachindex(levels),string.(levels)),
                    height=row_height,width=width,ytickformat="{:2.1f}")
        value = data[gid_ind,:]
        if isnothing(cs)
            rainclouds!(axis,groups,value;clouds=nothing,plot_boxplots=false,
                        side_nudge=0,markersize=point_size,alpha=alpha,
                        jitter_width=jitter_width)
        else
            rainclouds!(axis,groups,value;clouds=nothing,plot_boxplots=false,
                        side_nudge=0,color=cs,markersize=point_size,alpha=alpha,
                        jitter_width=jitter_width)
        end
        axis.yticks = [maximum(data[gid_ind,cell_id])]
        axis.yticklabelsize = font_size
        axis.yticksize = font_size * 0.3
        axis.spinewidth = 0.2
        axis.yaxisposition = y_axis_pos
        axis.ylabel = obj.var.name[gene_id[gid_ind]]
        axis.ylabelsize = isnothing(gene_label_size) ? 
            font_size * 1.2 : 
            gene_label_size
        axis.ylabelrotation = π * 2
        if gid_ind != length(gene_id)
            hidexdecorations!(axis)
        else
            hidexdecorations!(axis,ticklabels=false,label=false)
            axis.xticklabelrotation = π / 4
            axis.xticklabelsize = font_size * 1.2
        end
        hideydecorations!(axis,ticks=false,ticklabels=false,label=false)
        hidespines!(axis,:t,:r,:l)
    end
    rowgap!(figure.layout,0)
    Label(figure.layout[0,:],text=title,fontsize=title_size,font=:bold)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end

"""
    PropBar(obj)

Draw a bar plot of cell proportions.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `type::Symbol = :dodged`: figure type. :dodged or :stacked.
- `colorscheme::Symbol = :tableau_20`: group colorscheme from ColorSchemes 
  package.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `library_name::Union{Nothing,AbstractString} = nothing`: a column name in obs 
  for samples batch or nothing for one sample.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Integer = 640`: figure width.
- `height::Integer = 480`: figure height.
- `label_size::Integer = 12`: font size of labels.
- `label_color::Symbol = :black`: font color of labels.
- `alpha::AbstractFloat = 0.8`: bar alpha.
- `label_rotation::Real = -π / 2.5`: label rotation.
- `font_size::Integer = 16`: font size factor.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function PropBar(obj::WsObj;
        type::Symbol = :dodged,
        colorscheme::Symbol = :tableau_20,
        group_name::AbstractString = "clusters_latest",
        library_name::Union{Nothing,AbstractString} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Integer = 640,
        height::Integer = 480,
        label_size::Integer = 14,
        label_color::Symbol = :black,
        label_digits::Integer = 3,
        alpha::AbstractFloat = 0.8,
        label_rotation::Real = -π / 2.5,
        font_size::Integer = 16,
        path::Union{Nothing,AbstractString} = nothing)

    # Start a figure
    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))

    # Data preparation
    if isnothing(library_name)
        s = obj.obs
        s = countmap(s[!,group_name]) |> i -> sort(i;lt=natural) |> 
            j -> [ (key,j[key] / size(s,1)) for key in keys(j) ] |> 
            DataFrame
        rename!(s,["Groups","Props"])
        s[!,"library_name"] .= 1
        libraries = s.library_name
    else
        s = groupby(obj.obs,library_name) |> 
            x -> DataFrames.combine(x) do y
                ncell = size(y,1)
                z = countmap(y[!,group_name]) |> i -> sort(i;lt=natural)
                [ (key,z[key] / ncell) for key in keys(z) ]
            end
        s.Groups = [ x[1] for x in s.x1 ]
        s.Props = [ x[2] for x in s.x1 ]
        s = s[!,Not(:x1)]
        libraries = s[!,library_name] |> x -> Int.(indexin(x,unique(x)))
    end
    groups = s.Groups |> x -> Int.(indexin(x,unique(x)))

    # Colorscheme
    if colorscheme in keys(colorschemes)
        cs = colorschemes[colorscheme]
        cn = findmax(groups)[1]
        cs = cn <= length(cs) ? cgrad(colorscheme,cn;categorical=true) : nothing
        if isnothing(cs)
            return "Nothing to do! Color number is less than the group number!"
        end
    else
        cs = colorscheme
    end

    # Bar
    axis = Axis(figure[1,1],height=height,width=width,title=title,
                titlesize=title_size)
    if type == :dodged
        if isnothing(library_name)
            axis.xticklabelsvisible = false
            axis.xticksvisible = false
        else
            axis.xticks = (unique(libraries),unique(s[!,library_name]))
        end
        barplot!(axis,libraries,s.Props;dodge=s.Groups,color=s.Groups,
                 colormap=cs,strokecolor=:gray25,strokewidth=0.25,
                 bar_labels=:y,label_rotation=label_rotation,
                 label_color=label_color,label_size=label_size,alpha=alpha,
                 label_formatter=
                    l -> (l = round(l;digits=label_digits);"$l"))
        ylims!(axis,0,maximum(s.Props) + 0.05)
    elseif type == :stacked
        labels = round.(s.Props;digits=label_digits) |> x -> string.(x)
        if isnothing(library_name)
            cum_props = [ sum([ s.Props[j] for j in 1:i ]) 
                         for i in 1:size(s,1) ]
            axis.xticklabelsvisible = false
            axis.xticksvisible = false
        else
            cum_props = groupby(s,library_name) |> 
                x -> DataFrames.combine(x) do y
                    [ sum([ y.Props[j] for j in 1:i ]) for i in 1:size(y,1) ]
                end |> x -> x.x1
            axis.xticks = (unique(libraries),unique(obj.obs[!,library_name]))
        end
        barplot!(axis,libraries,s.Props;stack=s.Groups,color=s.Groups,
                 colormap=cs,strokecolor=:gray25,strokewidth=0.25,alpha=alpha)
        text!(ax,labels;
              position=Point2f.((libraries .- 1) + 
                                (repeat(range(0.55,1.4;
                                              length=length(unique(s.Groups))),
                                        length(unique(libraries)))) 
                                .^ 3 ./ 3 .* 0.8 .+ 0.5 .+ 
                                (abs(label_rotation) / 20),
                                cum_props .* 0.995),
              fontsize=label_size,color=label_color,rotation=label_rotation)
        axis.limits = ([0.5,length(unique(libraries)) + 0.5],[0,1.05])
        axis.yticks = ([0.,1.],["0","1.0"])
    else
        return "Nothing to do! Only :dodged and :stacked types are supported!"
    end
    axis.xgridvisible = false
    axis.xticklabelsize = font_size
    axis.yticklabelsize = font_size
    legend_elements = [ MarkerElement(color=c,marker=:rect,
                                      markersize=font_size)
                       for c in cs ]
    Legend(figure[1,2],legend_elements,string.(unique(s.Groups)))
    colgap!(figure.layout,10)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end
