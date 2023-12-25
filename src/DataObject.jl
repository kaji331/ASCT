mutable struct WsObj
    dat::Dict{AbstractString,Union{AbstractMatrix,AbstractSparseMatrix}}
    obs::DataFrame
    var::DataFrame
    log::Vector{AbstractString}
    meta::Union{Dict{AbstractString,
                     Union{AbstractVector,AbstractMatrix,AbstractString,Real,
                           DataFrame}},
                Nothing}
end
