#=
General noncommutative algebra system over a ring R.
The ring R must have methods of zero, one, +, -, *, conj, ==, isequal, hash
along with conversion and promotion methods as desired. Division by a scalar z is interpreted as
multiplication by inv(convert(R, z)) except with R = OneVarPolynomial{K} in which case it is
interpreted as multiplication by OneVarPolynomial{K}(inv(convert(K, z))).
=#

export AbstractTerm, AbstractVariable, Monomial, Polynomial, normal_order_generic, normal_order, expand_in_basis

#=
AbstractTerm
=#

# Supertype for AbstractVariable, Monomial, and Polynomial.
abstract type AbstractTerm{R} end

# We treat AbstractTerm subtypes as immutable types
copy(t::AbstractTerm) = t

function hash(t::U, h::UInt) where U <: AbstractTerm
    output = hash(nameof(U), h)
    for f in fieldnames(U)
        output = hash(getfield(t, f), output)
    end
    return output
end

function isequal(t1::U, t2::U) where U <: AbstractTerm
    for f in fieldnames(U)
        if !isequal(getfield(t1, f), getfield(t2, f))
            return false
        end
    end
    return true
end

isequal(::AbstractTerm, ::AbstractTerm) = false
isequal(::AbstractTerm{R}, ::Union{R, Number}) where R = false
isequal(::Union{R, Number}, ::AbstractTerm{R}) where R = false

function ==(t1::U, t2::U) where U <: AbstractTerm
    for f in fieldnames(U)
        if getfield(t1, f) != getfield(t2, f)
            return false
        end
    end
    return true
end

#=
AbstractVariable
=#

# Type for the variables in a Monomial or Polynomial.
# Any U <: AbstractVariable{R} must implement relation(u1::U, u2::U) -> Union{Polynomial{R, U}, Nothing}
# which encodes the relations to be used during normal ordering. If U is not always Hermitian
# or daggers to a polynomial, adjoint should be overridden. All fields in the type U must be of types that implement
# ==, isless, and hash.
abstract type AbstractVariable{R} <: AbstractTerm{R} end

# lexicographic ordering
function isless(u1::U, u2::U) where U <: AbstractVariable
    for f in fieldnames(U)
        f1, f2 = getfield(u1, f), getfield(u2, f)
        if isless(f1, f2)
            return true
        elseif isless(f2, f1)
            return false
        end
    end
    return false
end

# By default, AbstractVariables are Hermitian
adjoint(u::AbstractVariable) = u
*(u1::U, u2::V) where {R, U <: AbstractVariable{R}, V <: AbstractVariable{R}} = Monomial{R, U}(u1) * Monomial{R, V}(u2)

function prod(us::AbstractVector{U}) where {R, U <: AbstractVariable{R}}
    if isempty(us)
        return one(Polynomial{R, U})
    end
    return normal_order(us)
end

zero(::U) where {R, U <: AbstractVariable{R}} = zero(Polynomial{R, U})
zero(::Type{U}) where {R, U <: AbstractVariable{R}} = zero(Polynomial{R, U})
one(::U) where {R, U <: AbstractVariable{R}} = one(Monomial{R, U})
one(::Type{U}) where {R, U <: AbstractVariable{R}} = one(Monomial{R, U})
+(u1::U, u2::V) where {R, U <: AbstractVariable{R}, V <: AbstractVariable{R}} = Polynomial{R, U}(u1) + Polynomial{R, V}(u2)
sum(us::AbstractVector{U}) where {R, U <: AbstractVariable{R}} = sum(Polynomial{R, U}.(us))
*(z::Union{R, Number}, u::U) where {R, U <: AbstractVariable{R}} = convert(R, z) * Monomial{R, U}(u)
*(u::U, z::Union{R, Number}) where {R, U <: AbstractVariable{R}} = convert(R, z) * Monomial{R, U}(u)

#=
Monomial
=#

# A product of variables.
struct Monomial{R, U <: AbstractVariable{R}} <: AbstractTerm{R}
    variables::Vector{U} # should not be mutated after construction

    # Strict inner constructor for finer control of types.
    function Monomial{R, U}(variables::Vector{U}) where {R, U <: AbstractVariable{R}}
        return new{R, U}(variables)
    end
end

Monomial(m::Monomial) = m
Monomial(variables::Vector{U}) where {R, U <: AbstractVariable{R}} = Monomial{R, U}(variables)
Monomial(variables::AbstractVector{U}) where {R, U <: AbstractVariable{R}} = Monomial{R, U}(Vector{U}(variables))
Monomial{R, U}(variables::AbstractVector{U}) where {R, U <: AbstractVariable{R}} = Monomial{R, U}(Vector{U}(variables))
Monomial{R, U}(u::U) where {R, U <: AbstractVariable{R}} = Monomial{R, U}(Vector{U}([u]))
Monomial(u::U) where {R, U <: AbstractVariable{R}} = Monomial{R, U}(u)

function isless(m1::M, m2::M) where M <: Monomial
    l1, l2 = length(m1.variables), length(m2.variables)
    if l1 < l2
        return true
    elseif l1 > l2
        return false
    else # lexicographic ordering
        for ind in eachindex(m1.variables)
            v1, v2 = m1.variables[ind], m2.variables[ind]
            if isless(v1, v2)
                return true
            elseif isless(v2, v1)
                return false
            end
        end
        return false
    end
end

zero(::Monomial{R, U}) where {R, U} = zero(Polynomial{R, U})
zero(::Type{Monomial{R, U}}) where {R, U} = zero(Polynomial{R, U})
one(::Monomial{R, U}) where {R, U} = Monomial{R, U}(U[])
one(::Type{Monomial{R, U}}) where {R, U} = Monomial{R, U}(U[])
isone(m::Monomial) = isempty(m.variables)
# By default, adjoint of each variable is another variable
adjoint(m::Monomial{R, U}) where {R, U} = normal_order(Vector{U}(reverse(adjoint.(m.variables))))
*(m1::Monomial{R, U}, m2::Monomial{R, U}) where {R, U} = normal_order(Vector{U}(vcat(m1.variables, m2.variables)))

function prod(ms::AbstractVector{Monomial{R, U}}) where {R, U}
    if isempty(ms)
        return one(Polynomial{R, U})
    end
    return normal_order(Vector{U}(vcat([m.variables for m in ms]...)))
end

+(m1::Monomial{R, U}, m2::Monomial{R, U}) where {R, U} = Polynomial{R, U}(m1) + Polynomial{R, U}(m2)
sum(ms::AbstractVector{Monomial{R, U}}) where {R, U} = sum(Polynomial{R, U}.(ms))

function *(z::Union{R, Number}, m::Monomial{R, U}) where {R, U}
    z = convert(R, z)
    if iszero(z)
        return zero(Polynomial{R, U})
    end
    return Polynomial{R, U}(Dict(m => z))
end

*(m::Monomial{R, U}, z::Union{R, Number}) where {R, U} = z * m

#=
Polynomial
=#

# A linear combination of Monomials.
struct Polynomial{R, U <: AbstractVariable{R}} <: AbstractTerm{R}
    terms::Dict{Monomial{R, U}, R} # should not include zero coefficients

    # Strict inner constructor for finer control of types.
    function Polynomial{R, U}(terms::Dict{Monomial{R, U}, R}) where {R, U <: AbstractVariable{R}}
        return new{R, U}(terms)
    end
end

Polynomial(p::Polynomial) = p
Polynomial(terms::Dict{Monomial{R, U}, R}) where {R, U} = Polynomial{R, U}(terms)
Polynomial(terms::AbstractDict{Monomial{R, U}, S}) where {R, U, S} = Polynomial{R, U}(convert(Dict{Monomial{R, U}, R}, terms))
Polynomial{R, U}(m::Monomial{R, U}) where {R, U} = Polynomial{R, U}(Dict(m => one(R)))
Polynomial(m::Monomial{R, U}) where {R, U} = Polynomial{R, U}(m)
Polynomial{R, U}(u::U) where {R, U <: AbstractVariable{R}} = Polynomial{R, U}(Monomial{R, U}(u))
Polynomial(u::U) where {R, U <: AbstractVariable{R}} = Polynomial{R, U}(u)

zero(::Polynomial{R, U}) where {R, U} = Polynomial{R, U}(Dict{Monomial{R, U}, R}())
zero(::Type{Polynomial{R, U}}) where {R, U} = Polynomial{R, U}(Dict{Monomial{R, U}, R}())
iszero(p::Polynomial) = isempty(p.terms)
one(::Polynomial{R, U}) where {R, U} = Polynomial{R, U}(one(Monomial{R, U}))
one(::Type{Polynomial{R, U}}) where {R, U} = Polynomial{R, U}(one(Monomial{R, U}))

function adjoint(p::Polynomial{R, U}) where {R, U}
    if iszero(p)
        return zero(Polynomial{R, U})
    end
    return sum(Polynomial{R, U}[conj(coef) * m' for (m, coef) in p.terms])
end

function +(p1::P, p2::P) where P <: Polynomial
    return P(Dict(k => v for (k, v) in mergewith(+, p1.terms, p2.terms) if !iszero(v)))
end

function sum(ps::AbstractVector{P}) where P <: Polynomial
    if isempty(ps)
        return zero(P)
    end
    terms = empty(first(ps).terms)
    for p in ps
        mergewith!(+, terms, p.terms)
    end
    return P(Dict(k => v for (k, v) in terms if !iszero(v)))
end

function *(p1::Polynomial{R, U}, p2::Polynomial{R, U}) where {R, U}
    terms = Dict{Monomial{R, U}, R}()
    for (m1, coef1) in p1.terms
        for (m2, coef2) in p2.terms
            p = m1 * m2
            for (m, coef) in p.terms
                val = get(terms, m, zero(R)) + coef1 * coef2 * coef
                if iszero(val)
                    delete!(terms, m)
                else
                    terms[m] = val
                end
            end
        end
    end
    return Polynomial{R, U}(terms)
end

function prod(ps::AbstractVector{P}) where P <: Polynomial
    output = one(P)
    for p in ps
        output = output * p
    end
    return output
end

function *(z::Union{R, Number}, p::Polynomial{R, U}) where {R, U}
    z = convert(R, z)
    if iszero(z)
        return zero(Polynomial{R, U})
    end
    terms = sizehint!(Dict{Monomial{R, U}, R}(), length(p.terms))
    for (k, v) in p.terms
        w = z * v
        # z could be a zero-divisor
        if !iszero(w)
            terms[k] = w
        end
    end
    return Polynomial{R, U}(terms)
end

*(p::Polynomial{R, U}, z::Union{R, Number}) where {R, U} = z * p

#=
Normal ordering
=#

# generic normal ordering with no assumptions
function normal_order_generic(vars::AbstractVector{U}) where {R, U <: AbstractVariable{R}}
    for ind in eachindex(vars)[1:end-1]
        v_rel = relation(vars[ind], vars[ind+1])
        if !isnothing(v_rel)
            if iszero(v_rel)
                return zero(Polynomial{R, U})
            end
            vl, vr = vars[1:ind-1], vars[ind+2:end]
            return sum(Polynomial{R, U}[coef * normal_order_generic(Vector{U}(vcat(vl, m.variables, vr))) for (m, coef) in v_rel.terms])
        end
    end
    return Polynomial{R, U}(Monomial{R, U}(vars))
end

# by default, use normal_order_generic
normal_order(vars::AbstractVector{U}) where {R, U <: AbstractVariable{R}} = normal_order_generic(vars)

#=
Promotion, conversion, additional operations
=#

promote_rule(::Type{Polynomial{R, U}}, ::Type{Polynomial{R, U}}) where {R, U} = Polynomial{R, U}
promote_rule(::Type{Polynomial{R, U}}, ::Type{Polynomial{R, V}}) where {R, U, V} = Polynomial{R, Union{U, V}}
promote_rule(::Type{Polynomial{R, U}}, ::Type{Monomial{R, U}}) where {R, U} = Polynomial{R, U}
promote_rule(::Type{Polynomial{R, U}}, ::Type{Monomial{R, V}}) where {R, U, V} = Polynomial{R, Union{U, V}}
promote_rule(::Type{Polynomial{R, U}}, ::Type{U}) where {R, U} = Polynomial{R, U}
promote_rule(::Type{Polynomial{R, U}}, ::Type{V}) where {R, U, V <: AbstractVariable{R}} = Polynomial{R, Union{U, V}}
promote_rule(::Type{Polynomial{R, U}}, ::Type{R}) where {R, U} = Polynomial{R, U}
promote_rule(::Type{Monomial{R, U}}, ::Type{Monomial{R, U}}) where {R, U} = Monomial{R, U}
promote_rule(::Type{Monomial{R, U}}, ::Type{Monomial{R, V}}) where {R, U, V} = Monomial{R, Union{U, V}}
promote_rule(::Type{Monomial{R, U}}, ::Type{U}) where {R, U} = Monomial{R, U}
promote_rule(::Type{Monomial{R, U}}, ::Type{V}) where {R, U, V <: AbstractVariable{R}} = Monomial{R, Union{U, V}}
promote_rule(::Type{Monomial{R, U}}, ::Type{R}) where {U, R} = Polynomial{R, U}
promote_rule(::Type{U}, ::Type{R}) where {R, U <: AbstractVariable{R}} = Polynomial{R, U}

convert(::Type{Polynomial{R, U}}, p::Polynomial{R, U}) where {R, U} = p
convert(::Type{Polynomial{R, U}}, p::Polynomial{R, V}) where {R, U, V} = Polynomial{R, U}(convert(Dict{Monomial{R, U}, R}, p.terms))
convert(::Type{Polynomial{R, U}}, m::Monomial{R, U}) where {R, U} = Polynomial{R, U}(m)
convert(::Type{Polynomial{R, U}}, m::Monomial{R, V}) where {R, U, V} = Polynomial{R, U}(convert(Monomial{R, U}, m))
convert(::Type{Polynomial{R, U}}, u::U) where {R, U} = Polynomial{R, U}(u)
convert(::Type{Polynomial{R, U}}, u::V) where {R, U, V <: AbstractVariable{R}} = Polynomial{R, U}(convert(U, u))
convert(::Type{Polynomial{R, U}}, z::Union{R, Number}) where {R, U} = convert(R, z) * one(Monomial{R, U})
convert(::Type{Monomial{R, U}}, m::Monomial{R, U}) where {R, U} = m
convert(::Type{Monomial{R, U}}, m::Monomial{R, V}) where {R, U, V} = Monomial{R, U}(U[convert(U, v) for v in m.variables])
convert(::Type{Monomial{R, U}}, u::U) where {R, U} = Monomial{R, U}(u)
convert(::Type{Monomial{R, U}}, u::V) where {R, U, V <: AbstractVariable{R}} = Monomial{R, U}(convert(U, u))
convert(::Type{U}, u::U) where U <: AbstractVariable = u

function ==(t1::AbstractTerm, t2::AbstractTerm)
    if (t1 isa AbstractVariable) && (t2 isa AbstractVariable)
        return false
    end
    return ==(promote(t1, t2)...)
end

+(t1::AbstractTerm, t2::AbstractTerm) = +(promote(t1, t2)...)
*(t1::AbstractTerm, t2::AbstractTerm) = *(promote(t1, t2)...)
-(t::AbstractTerm{R}) where R = -one(R) * t
-(t1::AbstractTerm, t2::AbstractTerm) = t1 + (-t2)

==(t::AbstractTerm{R}, z::Union{R, Number}) where R = ==(promote(t, convert(R, z))...)
==(z::Union{R, Number}, t::AbstractTerm{R}) where R = ==(promote(convert(R, z), t)...)
+(t::AbstractTerm{R}, z::Union{R, Number}) where R = +(promote(t, convert(R, z))...)
+(z::Union{R, Number}, t::AbstractTerm{R}) where R = +(promote(convert(R, z), t)...)
-(t::AbstractTerm{R}, z::Union{R, Number}) where R = t + (-convert(R, z))
-(z::Union{R, Number}, t::AbstractTerm{R}) where R = z + (-t)

/(t::AbstractTerm{OneVarPolynomial{K}}, z::Number) where K = t * OneVarPolynomial{K}(inv(convert(K, z)))
/(t::AbstractTerm{R}, z::Number) where R = t * inv(convert(R, z))

function ^(t::AbstractTerm, num::Integer)
    @assert num >= 0
    output = Polynomial(one(t))
    while true
        if isodd(num)
            output *= t
        end
        num = fld(num, 2)
        if num > 0
            t *= t
        else
            return output
        end
    end
end

# Find the expansion of a polynomial in the given monomials or given basis.
function expand_in_basis(poly::Polynomial{R, U}, monomials::UniqueVector{Monomial{R, U}};
    monomial_to_basis::Union{AbstractMatrix, <:Nothing}=nothing) where {R, U <: AbstractVariable{R}}
    l = length(poly.terms)
    inds, vals = sizehint!(Int[], l), sizehint!(R[], l)
    for (m, coef) in poly.terms
        push!(inds, findfirst(isequal(m), monomials))
        push!(vals, coef)
    end
    v = sparsevec(inds, vals, length(monomials))
    return isnothing(monomial_to_basis) ? v : monomial_to_basis * v
end

function expand_in_basis(polys::AbstractVector{Polynomial{R, U}}, monomials::UniqueVector{Monomial{R, U}};
    monomial_to_basis::Union{AbstractMatrix, <:Nothing}=nothing) where {R, U <: AbstractVariable{R}}
    row_inds, col_inds, vals = Int[], Int[], R[]
    for (col_ind, poly) in enumerate(polys)
        v = expand_in_basis(poly, monomials, monomial_to_basis=monomial_to_basis)
        for (row_ind, val) in nonzero_entries(v)
            push!(row_inds, row_ind)
            push!(col_inds, col_ind)
            push!(vals, val)
        end
    end
    num_rows = isnothing(monomial_to_basis) ? length(monomials) : size(monomial_to_basis, 1)
    return sparse(row_inds, col_inds, vals, num_rows, length(polys))
end
