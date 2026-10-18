#=
Representation of real semidefinite programs in the form useful for solvers.

Primal problem:
maximize_{Y,y} b⋅y subject to
(By)_i + ∑_j tr(A_{i,j} * Y_j) = c_i for each i, Y_j ⪰ 0 for each j, y_k ≥ 0 for n_eq < k.
A is a matrix of symmetric real matrices, B is a real matrix, b is a real vector,
c is a real vector, Y is a vector of symmetric real matrices, and y is a real vector.

Dual problem:
minimize_x c⋅x subject to ∑_i x_i * A_{i,j} ⪰ 0 for each j,
(Bᵀx)_k = b_k for 1 ≤ k ≤ n_eq, (Bᵀx)_k ≥ b_k for n_eq < k.
x is a real vector.
=#

export SolverSDP, apply_preconditioner, ConicSDP, symmetric_from_lower_triangle

struct SolverSDP{T}
    A::Vector{SparseMatrixCSC{T, Int}} # A[j][i, :] is vectorization of lower triangular part of A_{i,j}
    B::SparseMatrixCSC{T, Int}
    b::SparseVector{T, Int}
    c::SparseVector{T, Int}
    n_eq::Int
    function SolverSDP{T}(A::Vector{SparseMatrixCSC{T, Int}}, B::SparseMatrixCSC{T, Int},
        b::SparseVector{T, Int}, c::SparseVector{T, Int}, n_eq::Int) where T
        len_b, len_c = length(b), length(c)
        @assert size(B, 2) == len_b >= n_eq >= 0
        @assert size(B, 1) == len_c >= 1
        for a in A
            @assert size(a, 1) == len_c
        end
        return new{T}(A, B, b, c, n_eq)
    end
end

function SolverSDP(A::AbstractVector{<:SparseMatrixCSC{T, Int}}, B::SparseMatrixCSC{T, Int},
    b::SparseVector{T, Int}, c::SparseVector{T, Int}; n_eq::Integer=length(b)) where T
    return SolverSDP{T}(collect(A), B, b, c, Int(n_eq))
end

# Order of indices in the vectorization of the lower triangular part of an n x n matrix.
lower_triangular_inds(n::Integer) = [(ind1, ind2) for ind2 in 1:n for ind1 in ind2:n]

function triangular_dim_from_length(l::Integer)
    n = fld(isqrt(1 + 8 * l) - 1, 2)
    @assert n * (n + 1) == 2 * l
    return n
end

# Produce a symmetric matrix from a lower triangular vectorization.
function symmetric_from_lower_triangle(v::AbstractVector{T}) where T
    n = triangular_dim_from_length(length(v))
    mat = zeros(T, n, n)
    v_ind, row_ind, col_ind = 1, 1, 1
    while v_ind <= length(v)
        mat[row_ind, col_ind] = v[v_ind]
        mat[col_ind, row_ind] = v[v_ind]
        v_ind += 1
        if row_ind < n
            row_ind += 1
        else
            col_ind += 1
            row_ind = col_ind
        end
    end
    return mat
end

# Write ∑_j tr(A_{i,j} * Y_j) as a matrix vector product using the lower triangular parts of A_{i,j} and Y_j.
function linearize_A(A::AbstractVector{<:SparseMatrixCSC{T}}, num_rows::Integer) where T
    m = sum(nnz, A; init=0)
    row_inds, col_inds, vals = sizehint!(Int[], m), sizehint!(Int[], m), sizehint!(T[], m)
    col_offset = 0
    for a in A
        @assert size(a, 1) == num_rows
        n = triangular_dim_from_length(size(a, 2))
        lt_inds = lower_triangular_inds(n)
        for (row_ind, col_ind, val) in nonzero_entries(a)
            push!(row_inds, row_ind)
            push!(col_inds, col_offset + col_ind)
            r, c = lt_inds[col_ind]
            push!(vals, val * (r == c ? 1 : 2))
        end
        col_offset += size(a, 2)
    end
    return sparse(row_inds, col_inds, vals, num_rows, col_offset)
end

# Map a function on a SolverSDP.
function map(f, ssdp::SolverSDP)
    B0 = map(f, ssdp.B)
    b0 = map(f, ssdp.b)
    c0 = map(f, ssdp.c)
    T = eltype(B0) # f can change the data type
    A0 = sizehint!(SparseMatrixCSC{T, Int}[], length(ssdp.A))
    for a in ssdp.A
        push!(A0, map(f, a))
    end
    return SolverSDP(A0, B0, b0, c0; n_eq=ssdp.n_eq)
end

# Apply a scaling to a SolverSDP.
function apply_preconditioner(ssdp::SolverSDP{T}, scales::Vector{<:Real}) where T
    if !all(x -> x > 0, scales)
        error("All scales must be positive.")
    end
    m = nnz(ssdp.c)
    inds, vals = sizehint!(Int[], m), sizehint!(T[], m)
    s_offset = 0
    for (i, ci) in nonzero_entries(ssdp.c)
        si = scales[s_offset + i]
        push!(inds, i); push!(vals, si * ci)
    end
    c = sparsevec(inds, vals, length(ssdp.c))
    m = nnz(ssdp.B)
    row_inds, col_inds, vals = sizehint!(Int[], m), sizehint!(Int[], m), sizehint!(T[], m)
    t_offset = s_offset + length(ssdp.c)
    for (i, k, Bik) in nonzero_entries(ssdp.B)
        si = scales[s_offset + i]
        tk = scales[t_offset + k]
        push!(row_inds, i); push!(col_inds, k); push!(vals, si * Bik * tk)
    end
    B = sparse(row_inds, col_inds, vals, size(ssdp.B)...)
    u_offset = t_offset + size(ssdp.B, 2)
    A = sizehint!(SparseMatrixCSC{T, Int}[], length(ssdp.A))
    for a in ssdp.A
        n = triangular_dim_from_length(size(a, 2))
        lt_inds = lower_triangular_inds(n)
        m = nnz(a)
        row_inds, col_inds, vals = sizehint!(Int[], m), sizehint!(Int[], m), sizehint!(T[], m)
        for (i, ind, Aijxy) in nonzero_entries(a)
            x, y = lt_inds[ind]
            si = scales[s_offset + i]
            ujx = scales[u_offset + x]
            ujy = scales[u_offset + y]
            push!(row_inds, i); push!(col_inds, ind); push!(vals, si * ujx * ujy * Aijxy)
        end
        u_offset += n
        push!(A, sparse(row_inds, col_inds, vals, size(a)...))
    end
    if length(scales) != u_offset + 1
        error("scales has the wrong length.")
    end
    v = scales[u_offset + 1]
    m = nnz(ssdp.b)
    inds, vals = sizehint!(Int[], m), sizehint!(T[], m)
    for (k, bk) in nonzero_entries(ssdp.b)
        tk = scales[t_offset + k]
        push!(inds, k); push!(vals, v * bk * tk)
    end
    b = sparsevec(inds, vals, size(ssdp.b)...)
    return SolverSDP(A, B, b, c; n_eq=ssdp.n_eq), v
end

#=
Representation of real semidefinite programs in the form useful for primal iterative refinement and exact verification.

Primal problem:
maximize_v d⋅v subject to v ∈ C, R v = r.
C is a product of free cones (z is real), nonnegative cones (z ≥ 0), and PSD cones (Z ⪰ 0), in that order.
R is a real matrix and d and r are real vectors.
=#

struct ConicSDP{T}
    R::SparseMatrixCSC{T, Int}
    r::SparseVector{T, Int}
    d::SparseVector{T, Int}
    num_free::Int
    num_nonneg::Int
    psd_dims::Vector{Int}
    function ConicSDP{T}(R::SparseMatrixCSC{T, Int}, r::SparseVector{T, Int}, d::SparseVector{T, Int},
        num_free::Int, num_nonneg::Int, psd_dims::Vector{Int}) where T
        @assert num_free >= 0
        @assert num_nonneg >= 0
        @assert all(psd_dims .> 0)
        num_psd = sum(div(n * (n + 1), 2) for n in psd_dims; init=0)
        @assert size(R, 2) == length(d) == num_free + num_nonneg + num_psd
        @assert size(R, 1) == length(r)
        return new{T}(R, r, d, num_free, num_nonneg, psd_dims)
    end
end

function ConicSDP(R::SparseMatrixCSC{T, Int}, r::SparseVector{T, Int}, d::SparseVector{T, Int},
        num_free::Int, num_nonneg::Int, psd_dims::Vector{Int}) where T
    return ConicSDP{T}(R, r, d, num_free, num_nonneg, psd_dims)
end

function ConicSDP(ssdp::SolverSDP{T}) where T
    lin_A = linearize_A(ssdp.A, size(ssdp.B, 1))
    R = hcat(ssdp.B, lin_A)
    r = ssdp.c
    d = vcat(ssdp.b, spzeros(T, size(lin_A, 2)))
    num_free = ssdp.n_eq
    num_nonneg = length(ssdp.b) - num_free
    psd_dims = [triangular_dim_from_length(size(a, 2)) for a in ssdp.A]
    return ConicSDP{T}(R, r, d, num_free, num_nonneg, psd_dims)
end
