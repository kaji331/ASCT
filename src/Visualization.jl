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
        width::Integer = 600,
        height::Integer = 300,
        path::Union{Nothing,AbstractString} = nothing)

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
    axis2,heat2 = heatmap(grid_layout[2,1],reverse(data',dims=2),
                          colormap=colorschemes[colorscheme][8:-1:1])

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

Draw a scatter figure of features.

# Arguments
- `obj::WsObj`: a single-cell WsObj struct.

# Keyword Arguments
- `dimension_name::AbstractString = "pca"`: name of dimension data in meta.
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
        dimension_name::AbstractString = "pca",
        group_name::Union{Nothing,AbstractString} = "clusters_latest",
        split_name::Union{Nothing,AbstractString} = nothing,
        colorscheme::Symbol = :tableau_20,
        dim_prefix::Union{Nothing,AbstractString} = nothing,
        point_size::Integer = 7,
        alpha::AbstractFloat = 0.7,
        title::AbstractString = "",
        title_size::Integer = 30,
        width::Integer = 400,
        height::Integer = 300,
        path::Union{Nothing,AbstractString} = nothing)

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
    data = obj.meta[dimension_name]
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
        cs = [ (colorscheme,alpha) 
              for alpha in range(0.1,1;length=length(levels)) ]
    end
    figure.scene.theme.palette.color = cs

    # Plots
    final_axis = nothing
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
            hidespines!(axis,:t,:r,:l,:b)
            hidedecorations!(axis,grid=false)

            # Scatters
            for level in levels
                idx = (groups .== level) .& (facets .== facet_levels[counter])
                scatter!(axis,data[:,1][idx],data[:,2][idx],alpha=alpha,
                         label=string(level),markersize=point_size)
            end

            # Arrows
            if row == facet_rows && col == 1
                end_len1 = (range(xlim[1],xlim[2],length=10)[2] - xlim[1]) * 
                    facet_cols
                end_len2 = (range(ylim[1],ylim[2],length=10)[2] - ylim[1]) * 
                    facet_rows
                factor = facet_rows == facet_cols == 1 ? 0.2 : 0.5
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

            final_axis = axis
            counter += 1
        end
    end
    Label(grid_layout_main[0,:],text=title,fontsize=title_size,font=:bold)
    Label(grid_layout_main[end + 1,:],text=dim_prefix * "_1",
          fontsize=title_size / 2,halign=:left)
    Label(grid_layout_main[end - 1,0],text=dim_prefix * "_2",
          fontsize=title_size / 2,rotation=1.6,valign=:bottom)
    # Make the fixed size of legend elements 
    legend_elements = [ MarkerElement(color=c,marker=:circle,makersize=15) 
                       for c in figure.scene.theme.palette.color[] ]
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
        hex_colorscheme::Symbol = :thermal,
        hex_number::Union{Integer,Tuple{Integer,Integer}} = 75,
        density_color::Tuple{Symbol,AbstractFloat} = (:orangered,0.3),
        height::Integer = 300,
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

    figure = Figure(fonts=(;regular="Nimbus Roman",bold="Nimbus Roman Bold"))
    if title == :auto
        title = "Elbow figure of PCA"
    end
    axis = Axis(figure[1,1],title=title,titlesize=title_size,titlefont=:bold,
                xlabel="PCs",ylabel="Percentages of variances",width=width,
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
- `width::Union{Symbol,Integer} = :auto`: figure width.
- `row_height::Integer = 30`: row height for each gene.
- `y_axis_pos::Symbol = :right`: set the position of gene labels.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureViolin(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        colorscheme::Symbol = :Hiroshige,
        group_name::AbstractString = "clusters_latest",
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Union{Symbol,Integer} = :auto,
        row_height::Integer = 30,
        y_axis_pos::Symbol = :right,
        path::Union{Nothing,AbstractString} = nothing)

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
        cs = 1 < length(levels) <= length(cs) ? 
            cs[round.(Integer,range(1,length(cs);length=length(levels)))] : 
            nothing
    else
        cs = repeat([colorscheme],length(levels))
    end

    ll = length(levels)
    # Violin for each gene
    for gid_ind in eachindex(gene_id)
        axis = Axis(figure[gid_ind,1],xlabel=group_name,
                    xticks=(eachindex(levels),string.(levels)),
                    height=row_height,
                    width=width == :auto ? row_height * ll * 1.5 : width,
                    ytickformat="{:2.1f}")
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
        axis.yticklabelsize = row_height * 0.3
        axis.yticksize = row_height * 0.1
        axis.spinewidth = 0.2
        axis.yaxisposition = y_axis_pos
        axis.ylabel = obj.var.name[gene_id[gid_ind]]
        axis.ylabelsize = isnothing(gene_label_size) ? 
            row_height * 0.4 : 
            gene_label_size
        axis.ylabelrotation = π * 2
        if gid_ind != length(gene_id)
            hidexdecorations!(axis)
        else
            hidexdecorations!(axis,ticklabels=false,label=false)
            axis.xticklabelrotation = π / 4
            axis.xticklabelsize = row_height * 0.4
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
- `colorscheme::Symbol = :red`: just one color for expression levels.
- `group_name::AbstractString = "clusters_latest"`: a column name in obs for 
  cell groups.
- `gene_label_size::Union{Nothing,Integer} = nothing`: set the size of gene 
  labels or switch off gene labels.
- `title::AbstractString = ""`: title of figure.
- `title_size::Integer = 20`: size of title font.
- `width::Union{Symbol,Integer} = :auto`: figure width.
- `row_height::Integer = 30`: row height for each gene.
- `y_axis_pos::Symbol = :right`: set the position of gene labels.
- `path::Union{Nothing,AbstractString} = nothing`: figure path.
"""
function FeatureFracDots(obj::WsObj;
        features::Union{Integer,Vector{<: AbstractString}} = 2,
        colorscheme::Symbol = :red,
        group_name::AbstractString = "clusters_latest",
        gene_label_size::Union{Nothing,Integer} = nothing,
        title::AbstractString = "",
        title_size::Integer = 20,
        width::Union{Symbol,Integer} = :auto,
        row_height::Integer = 30,
        y_axis_pos::Symbol = :right,
        path::Union{Nothing,AbstractString} = nothing)

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
    else
        # Single color with different alpha
        cs = range(0;stop=1,length=256)
    end
    color_start = range(minimum(data);stop=maximum(data),length=length(cs))

    ll = length(levels)
    lg = length(gene_id)
    # Dots for each gene
    for gid_ind in eachindex(gene_id)
        axis = Axis(figure[gid_ind,1],xlabel=group_name,
                    xticks=(eachindex(levels),string.(levels)),
                    height=row_height,
                    width=width == :auto ? row_height * ll * 1.5 : width)
        # Dots for each group
        for (i,level) in enumerate(levels)
            idx = groups .== level
            value = data[gid_ind,idx]
            value_mean = mean(value)
            percentage = count(value .!= 0) / count(idx)
            color = typeof(cs) <: AbstractRange ? 
                (colorscheme,cs[findlast(x -> value_mean >= x,color_start)]) : 
                cs[findlast(x -> value_mean >= x,color_start)]
            scatter!(axis,i,5;datalimits=(0,Inf),color=color,strokewidth=0.5,
                     strokecolor=:black,markersize=(percentage + 0.05) * 30)
        end
        axis.yaxisposition = y_axis_pos
        axis.ylabel = obj.var.name[gene_id[gid_ind]]
        axis.ylabelsize = isnothing(gene_label_size) ? 
            row_height * 0.4 : 
            gene_label_size
        axis.ylabelrotation = π * 2
        if gid_ind != lg
            hidexdecorations!(axis)
        else
            hidexdecorations!(axis,ticklabels=false,label=false)
            axis.xticklabelrotation = π / 4
            axis.xticklabelsize = row_height * 0.4
        end
        hideydecorations!(axis,label=false)
        hidespines!(axis,:t,:b,:l,:r)
    end
    rowgap!(figure.layout,0)
    Label(figure.layout[0,:],text=title,fontsize=title_size,font=:bold)
    colormap = typeof(cs) <: AbstractRange ? 
        cgrad([colorscheme],256;alpha=(0,1)) : 
        cs
    Colorbar(figure.layout[round(Integer,lg * 0.75):lg,
                           end + 1],
             limits=(0,maximum(data)),colormap=colormap,tickformat="{:2.1f}",
             ticks=[0,maximum(data)],ticklabelsize=row_height * 0.4)
    legend_elements = [ MarkerElement(color=:grey80,strokewidth=0.5,
                                      strokecolor=:black,marker=:circle,
                                      markersize=size)
                       for size in 
                       range(0;stop=row_height,length=6)[2:6]]
    Legend(figure.layout[round(Integer,lg * 0.45):round(Integer,lg * 0.7),end],
           legend_elements,["0.2","0.4","0.6","0.8","1.0"],labelsize=12)

    # Save figure
    if !isnothing(path)
        save(path,figure;pt_per_unit=2,px_per_unit=2)
    end

    return figure
end
