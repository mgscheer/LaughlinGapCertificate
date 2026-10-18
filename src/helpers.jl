# Iteration through a sparse matrix.
function nonzero_entries(A::SparseMatrixCSC; check_zeros::Bool=true)
    rows = rowvals(A)
    vals = nonzeros(A)
    if check_zeros
        return ((rows[j], col, vals[j]) for col in axes(A, 2) for j in nzrange(A, col) if !iszero(vals[j]))
    else
        return ((rows[j], col, vals[j]) for col in axes(A, 2) for j in nzrange(A, col))
    end
end

# Iteration through a sparse vector.
function nonzero_entries(v::SparseVector; check_zeros::Bool=true)
    inds = rowvals(v)
    vals = nonzeros(v)
    if check_zeros
        return ((inds[j], vals[j]) for j in eachindex(vals) if !iszero(vals[j]))
    else
        return ((inds[j], vals[j]) for j in eachindex(vals))
    end
end
