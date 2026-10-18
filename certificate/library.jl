using LaughlinGapCertificate, Combinatorics, BigCombinatorics, UniqueVectors, Distributed
using ProgressMeter, SparseArrays, LinearAlgebra, IntervalArithmetic
import LaughlinGapCertificate: scaled_fermion_anticommutator

K = Rational{BigInt} # The rational numbers are the coefficient field
R = OneVarPolynomial{K} # Polynomial ring
t = generator_power(K, 1) # Polynomial ring generator
F = ScaledFermion{R, Int} # Scaled fermion operators

#=
Creation of the SDP
=#

# ξ_m polynomial.
function scaled_fermion_anticommutator(m::Integer)
    ξ = Factorial(m) * one(R)
    for n in 1:(m-1)
        ξ *= 1 - n * t
    end
    return ξ
end

# f and f^† variables.
function create_variables(r::Integer)
    ms = 0:r
    fs = Dict(m => F(m, false) for m in ms)
    variables = [fs[m] for m in ms]
    append!(variables, [fs[m]' for m in ms])
    return fs, variables
end

# q_M operator.
function create_q(M::Integer, fs::Dict{Int, F})
    q = zero(Polynomial{R, F})
    M0 = BigInt(M)
    for m in 0:M
        m0 = BigInt(m)
        coef = (2 * m0 - M0) // (2 * Factorial(m0) * Factorial(M0 - m0))
        q += coef * fs[M - m] * fs[m]
    end
    return q
end

# The polynomial η_{M_num} / η_{M_den}.
function η_ratio(M_num::Integer, M_den::Integer)
    @assert M_num >= M_den >= 1
    p = one(R)
    for l in 2:M_num
        if !(2 <= l <= M_den)
            p *= 2 - l * t
        end
    end
    return p
end

# We represent (p' + p)/2 by one of p' and p, which we choose canonically using isless.
function orbit_rep(m::Monomial{R, F})
    m_dag = m'
    if isless(m, m_dag)
        return Polynomial(m)
    else
        return Polynomial(m_dag)
    end
end

# Apply orbit_rep to each term in a polynomial.
function orbit_rep(p::Polynomial{R, F})
    if iszero(p)
        return zero(Polynomial{R, F})
    else
        return sum([coef * orbit_rep(m) for (m, coef) in p.terms])
    end
end

# \tilde{η} * W after application of orbit_rep.
function create_W(r::Integer, qs::Dict{Int, Polynomial{R, F}}, fs::Dict{Int, F})
    η_ratios = [η_ratio(r, M) for M in 1:r]
    W = zero(Polynomial{R, F})
    for M in 1:r
        for N in 1:M
            coef = (M == N ? 1 : 2) * η_ratios[M] * η_ratios[N]
            W += coef * orbit_rep(qs[M]' * fs[N]' * fs[M] * qs[N])
        end
    end
    return W * (1 + t) / 2^(2 * BigInt(r) - 1)
end

# Compute the coefficient α in the equation \tilde{η} * \tilde{Q}^\dagger_M \tilde{Q}_M = α * q^†_M q_M.
function create_Q_dag_Q_coef(M::Integer, r::Integer)
    return (Factorial(M - 1) // 2^(2 * BigInt(r) - 1)) * (2 - t) * η_ratio(r, M) * η_ratio(r, 1)
end

# Compute the coefficient β in the equation
# \tilde{η} * \tilde{Q}^\dagger_{M} \tilde{Q}^\dagger_N \tilde{Q}_N \tilde{Q}_{M} = β * q^†_{M} q^†_N q_N q_{M}.
function create_QQ_dag_QQ_coef(M::Integer, N::Integer, r::Integer)
    return (Factorial(M - 1) * Factorial(N - 1) // 2^(2 * BigInt(r) - 1)) * (2 - t) * η_ratio(r, M) * η_ratio(r, N)
end

# Compute [J_+, v] where v is either f^†_m or f_m for some m.
# Treats [J_+, f^†_r] = 0 since we will use J_plus_comm_in_monomials to be sure not to encounter this case.
function J_plus_comm(v::F, fs::Dict{Int, F})
    m = v.label
    @assert m >= 0
    if v.dagger && (m + 1 in keys(fs))
        return Polynomial{R, F}(fs[m + 1]')
    elseif !v.dagger && (m >= 1) && (m - 1 in keys(fs))
        m = BigInt(m)
        return -m * (1 - (m - 1) * t) * fs[m - 1]
    else
        return zero(Polynomial{R, F})
    end
end

# For m = O_1 ⋯ O_N where each O is an f^† or f, compute [J_+, m] = ∑_j O_1 ⋯ [J_+, O_j] ⋯ O_N.
function J_plus_comm(m::Monomial{R, F}, fs::Dict{Int, F})
    output = zero(Polynomial{R, F})
    for j in eachindex(m.variables)
        output += Monomial(m.variables[1:(j-1)]) * J_plus_comm(m.variables[j], fs) * Monomial(m.variables[(j+1):end])
    end
    return output
end

# Determine efficiently if orbit_rep([J_+, m]) will be in the span of the given monomials.
function J_plus_comm_in_monomials(m::Monomial{R, F}, monomials::AbstractVector{Monomial{R, F}})
    function is_normal_ordered(u1::F, u2::F)
        if u1.dagger != u2.dagger
            return u1.dagger
        elseif u1.dagger
            return u1.label < u2.label
        else
            return u1.label > u2.label
        end
    end

    for j in eachindex(m.variables)
        v = m.variables[j]
        new_label = v.label + (v.dagger ? 1 : -1)
        if new_label < 0
            # [J_+, O_j] = 0, so this term vanishes
            continue
        end
        new_v = F(new_label, v.dagger)
        # now find where [J_+, O_j] goes in the product without tracking coefficients
        new_vars = F[]
        other_inds = filter(!isequal(j), eachindex(m.variables))
        k = 1
        while (k <= length(other_inds)) && is_normal_ordered(m.variables[other_inds[k]], new_v)
            push!(new_vars, m.variables[other_inds[k]])
            k += 1
        end
        if (k <= length(other_inds)) && (m.variables[other_inds[k]] == new_v)
            # this term has a repeated f^† or f and therefore vanishes
            continue
        end
        push!(new_vars, new_v)
        for other_v in m.variables[other_inds[k:end]]
            @assert is_normal_ordered(new_v, other_v)
            push!(new_vars, other_v)
        end
        new_m = Monomial{R, F}(new_vars)
        if !(new_m in monomials) && !(new_m' in monomials)
            # found a term not in monomials
            return false
        end
    end
    # all terms are in monomials
    return true
end

# Particle number and J_z quantum numbers (J_z is angular momentum relative to the south pole).
function get_quantum_numbers(m::Monomial{R, F})
    n, j = 0, 0
    for v in m.variables
        if v.dagger
            n += 1
            j += v.label
        else
            n -= 1
            j -= v.label
        end
    end
    return n, j
end

# Add all normal ordered monomials of a given degree to the monomials dict, which maps
# quantum number tuples to sets of monomials.
function add_monomials!(monomials::Dict{Tuple{Int, Int}, Set{Monomial{R, F}}},
    vars::AbstractVector{F}, degree::Integer)
    for var_combo in combinations(vars, degree)
        poly = normal_order(var_combo)
        for m in keys(poly.terms)
            if length(m.variables) == degree
                quantum_numbers = get_quantum_numbers(m)
                if !(quantum_numbers in keys(monomials))
                    monomials[quantum_numbers] = Set{Monomial{R, F}}()
                end
                push!(monomials[quantum_numbers], m)
            end
        end
    end
    return monomials
end

# Add a polynomial to the polynomials dict, which maps quantum number tuples to vectors of polynomials.
function add_polynomial!(polynomials::Dict{Tuple{Int, Int}, Vector{Polynomial{R, F}}}, p::Polynomial{R, F})
    qns = get_quantum_numbers.(keys(p.terms))
    if !isempty(qns)
        qn = qns[1]
        @assert all([x == qn for x in qns])
        if !(qn in keys(polynomials))
            polynomials[qn] = Polynomial{R, F}[]
        end
        push!(polynomials[qn], p)
    end
    return polynomials
end

function setup_SDP(r_prime::Integer, r::Integer)
    fs, variables = create_variables(r)
    qs = Dict(M => create_q(M, fs) for M in 1:r)
    W = create_W(r, qs, fs)
    normalization = create_Q_dag_Q_coef(1, r) * orbit_rep(qs[1]' * qs[1])

    # Create a dict mapping symmetry sectors to polynomials.
    polynomials = Dict{Tuple{Int, Int}, Vector{Polynomial{R, F}}}()
    for M in 1:r_prime
        add_polynomial!(polynomials, qs[M])
        for m in 0:r_prime
            for n in 0:r_prime
                if 1 <= -m + n + M <= r_prime
                    add_polynomial!(polynomials, fs[m]' * fs[n] * qs[M])
                end
            end
            for n in (m+1):r_prime
                if m + n + M <= 2 * r_prime + 4
                    add_polynomial!(polynomials, fs[m] * fs[n] * qs[M])
                end
            end
        end
    end
    for M in 1:r
        add_polynomial!(polynomials, fs[M]' * qs[M])
        add_polynomial!(polynomials, fs[M] * qs[M])
    end
    P_keys = sort(collect(keys(polynomials)))
    P_arrs = [collect(polynomials[tup]) for tup in P_keys]
    println("length(P_arrs): $(length(P_arrs))")
    println("sum(length.(P_arrs)): $(sum(length.(P_arrs)))")
    flush(stdout)

    # Perform all distinct products of monomials needed.
    P_mono_arrs = [UniqueVector(sort(collect(union([keys(p.terms) for p in P_arr]...)))) for P_arr in P_arrs]
    prod_keys = vcat([[(P_ind, ind1, ind2) for ind1 in eachindex(P_mono) for ind2 in eachindex(P_mono)]
        for (P_ind, P_mono) in enumerate(P_mono_arrs)]...)
    println("length(prod_keys): $(length(prod_keys))")
    flush(stdout)

    pool = CachingPool(workers())
    P_mono_prods = @showprogress desc="Multiplying: " pmap(pool, prod_keys; batch_size=1000) do (P_ind, ind1, ind2)
        orbit_rep(P_mono_arrs[P_ind][ind1]' * P_mono_arrs[P_ind][ind2])
    end
    clear!(pool)
    monomial_set = Set{Monomial{R, F}}()
    for p in P_mono_prods
        union!(monomial_set, keys(p.terms))
    end
    monomials = UniqueVector(sort(collect(monomial_set)))
    println("length(monomials): $(length(monomials))")
    flush(stdout)

    # Expand each polynomial into monomials and uses the products above
    # to compute the structure constants.
    cutoff_inds_P = [0; cumsum([length(P_mono)^2 for P_mono in P_mono_arrs])]
    A = @showprogress desc="Expanding p1' * p2 in monomials: " pmap(pool, eachindex(P_arrs)) do P_ind
        C = expand_in_basis(P_arrs[P_ind], P_mono_arrs[P_ind])
        D = expand_in_basis(P_mono_prods[(cutoff_inds_P[P_ind] + 1):(cutoff_inds_P[P_ind+1])], monomials)
        M = D * kron(C, C)
        N = length(P_arrs[P_ind])
        inds = [k1 + N * (k0-1) for k1 in 1:N for k0 in k1:N]
        M[:, inds]
    end
    clear!(pool)

    # Construct a vector of monomials which will be used for SU(2) symmetry.
    P_plus = UniqueVector{Monomial{R, F}}()
    for (r0, deg_max) in [(r_prime, 4), (r, 3)]
        f_dag_vars = filter(v -> v.dagger && (v.label <= r0), variables)
        f_vars = filter(v -> (!v.dagger) && (v.label <= r0), variables)
        f_dag_monomials = Dict{Tuple{Int, Int}, Set{Monomial{R, F}}}()
        f_monomials = Dict{Tuple{Int, Int}, Set{Monomial{R, F}}}()
        for degree in 1:deg_max
            add_monomials!(f_dag_monomials, f_dag_vars, degree)
            add_monomials!(f_monomials, f_vars, degree)
        end
        for (n, j) in keys(f_dag_monomials)
            if (-n, -1 - j) in keys(f_monomials)
                for op_f_dag in f_dag_monomials[(n, j)]
                    for op_f in f_monomials[(-n, -1 - j)]
                        p = op_f_dag * op_f
                        m = only(keys(p.terms))
                        if !(m in P_plus)
                            if J_plus_comm_in_monomials(m, monomials)
                                push!(P_plus, m)
                            end
                        end
                    end
                end
            end
        end
    end
    P_plus = UniqueVector(sort(collect(P_plus)))
    println("length(P_plus): $(length(P_plus))")
    flush(stdout)

    equation_polys = @showprogress desc="Computing commutators with J_+: " pmap(pool, P_plus; batch_size=200) do m
        orbit_rep(J_plus_comm(m, fs))
    end
    clear!(pool)

    QHQ = sum([create_QQ_dag_QQ_coef(r_prime, N, r) * orbit_rep(
        qs[r_prime]' * qs[N]' * qs[N] * qs[r_prime]) for N in 1:r_prime])
    objective = QHQ - W + normalization

    # Check that the objective and normalization polynomials are in the span of the monomials.
    @assert issubset(union([keys(p.terms) for p in [objective, normalization]]...), Set(monomials))

    B = expand_in_basis([normalization; equation_polys], monomials)
    b = spzeros(R, size(B, 2))
    b[1] = 1
    c = expand_in_basis(objective, monomials)
    return SolverSDP(A, B, b, c)
end

#=
SDP certification helpers
=#

# Check transpose(lhs) * x == rhs using efficient sparse matrix traversal and early failure.
function fast_equality_check(lhs::SparseMatrixCSC{T, Int}, x::AbstractVector{T}, rhs::AbstractVector{T}; rhs_scale::T=one(T)) where T
    if size(lhs) != (length(x), length(rhs))
        return false
    end
    rows, vals = rowvals(lhs), nonzeros(lhs)
    @showprogress desc="Checking equality: " for j in axes(lhs, 2)
        val = zero(T)
        for k in nzrange(lhs, j)
            val += vals[k] * x[rows[k]]
        end
        if val != rhs_scale * rhs[j]
            return false
        end
    end
    return true
end

# Exact positivity check with vals, vecs certificate.
function positivity_bound(Y::Matrix{Rational{BigInt}}, vals::Vector{Rational{BigInt}},
    vecs::Matrix{Rational{BigInt}})
    d0 = minimum(vals)
    if d0 <= 0
        return nothing
    end
    σ = maximum(sum(abs, Y - vecs * Diagonal(vals) * transpose(vecs); dims=2))
    τ = maximum(sum(abs, vecs * transpose(vecs) - I; dims=2))
    bound = d0 * (1 - τ) - σ
    if bound <= 0
        return nothing
    end
    return bound
end

# Exact positivity check for all PSD variables in v.
function verify_positivity(csdp::ConicSDP{Rational{BigInt}}, v::Vector{Rational{BigInt}},
    certificate::Vector{Tuple{Vector{Rational{BigInt}}, Matrix{Rational{BigInt}}}})
    if length(certificate) != length(csdp.psd_dims)
        return false
    end
    for ind in (csdp.num_free + 1):(csdp.num_free + csdp.num_nonneg)
        if v[ind] <= 0
            return false
        end
    end
    v_ranges = Vector{Rational{BigInt}}[]
    ind = csdp.num_free + csdp.num_nonneg
    for n in csdp.psd_dims
        inc = div(n * (n + 1), 2)
        push!(v_ranges, v[(ind + 1):(ind + inc)])
        ind += inc
    end
    if ind != length(v)
        return false
    end
    order = sortperm(csdp.psd_dims; rev=true) # the biggest problems will start first
    bounds = @showprogress desc="Checking positivity: " pmap(zip(certificate[order], v_ranges[order])) do ((vals, vecs), v_range)
        positivity_bound(symmetric_from_lower_triangle(v_range), vals, vecs)
    end
    return all(!isnothing, bounds)
end

#=
Interval arithmetic tail bound helpers
=#

int(n) = interval(Float64, n)

# ∑_{n in a:b} sqrt(n) * x^(n+1).
function partial_χ(x, a, b)
    χ = int(0)
    for n in a:b
        χ += sqrt(int(n)) * pown(x, n + 1)
    end
    return χ
end

# Interval arithmetic calculation of omega_tilde_2.
function ω_tilde_2_interval(r; N0=200)
    if !(1 <= r <= N0)
        error("we require 1 <= r <= N0.")
    end
    x0 = inv(sqrt(int(2)))
    χ_r = partial_χ(x0, 1, r)
    finite_tail = partial_χ(x0, r + 1, N0)
    tail_upper = pown(x0, N0 + 2) * (int(N0 + 1) - int(N0) * x0) / pown(int(1) - x0, 2)
    return int(4) * χ_r * (finite_tail + tail_upper)
end
