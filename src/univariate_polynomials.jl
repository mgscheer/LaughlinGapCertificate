#=
Polynomials of a single variable over a field K.
The field K must have methods of +, -, *, inv, zero, iszero, one, ==, hash
along with conversion and promotion methods as desired.
By default, we assume coefficients are real and are invariant under conj.
If that is not true, then conj methods must be overridden.
=#

export OneVarPolynomial, degree, leading_coefficient, get_coefficient, get_value, generator_power

struct OneVarPolynomial{K}
    coefficients::Vector{K} # should not include trailing zeros or be mutated after construction
end

# We treat OneVarPolynomial as an immutable type
copy(p::OneVarPolynomial) = p

# By default, the variable is considered real and invariant under conj
conj(p::OneVarPolynomial) = p
adjoint(p::OneVarPolynomial) = conj(p)
transpose(p::OneVarPolynomial) = p

function hash(p::U, h::UInt) where U <: OneVarPolynomial
    output = hash(nameof(U), h)
    output = hash(p.coefficients, output)
    return output
end

function isequal(p1::U, p2::U) where U <: OneVarPolynomial
    return isequal(p1.coefficients, p2.coefficients)
end

isequal(::OneVarPolynomial{K}, ::Union{K, Number}) where K = false
isequal(::Union{K, Number}, ::OneVarPolynomial{K}) where K = false

function ==(p1::U, p2::U) where U <: OneVarPolynomial
    return p1.coefficients == p2.coefficients
end

zero(::OneVarPolynomial{K}) where K = OneVarPolynomial{K}(K[])
zero(::Type{OneVarPolynomial{K}}) where K = OneVarPolynomial{K}(K[])
iszero(p::OneVarPolynomial) = isempty(p.coefficients)
one(::OneVarPolynomial{K}) where K = OneVarPolynomial{K}(K[one(K)])
one(::Type{OneVarPolynomial{K}}) where K = OneVarPolynomial{K}(K[one(K)])
OneVarPolynomial{K}(p::OneVarPolynomial{K}) where K = p
OneVarPolynomial(p::OneVarPolynomial) = p
OneVarPolynomial{K}(z::K) where K = iszero(z) ? zero(OneVarPolynomial{K}) : OneVarPolynomial{K}(K[z])
OneVarPolynomial(z::K) where K = OneVarPolynomial{K}(z)

function trim_zeros!(coefficients::Vector)
    while !isempty(coefficients) && iszero(coefficients[end])
        pop!(coefficients)
    end
    return coefficients
end

function +(p1::OneVarPolynomial{K}, p2::OneVarPolynomial{K}) where K
    coefficients = zeros(K, max(length(p1.coefficients), length(p2.coefficients)))
    for coefs in (p1.coefficients, p2.coefficients)
        for j in eachindex(coefs)
            coefficients[j] += coefs[j]
        end
    end
    return OneVarPolynomial{K}(trim_zeros!(coefficients))
end

function *(p1::OneVarPolynomial{K}, p2::OneVarPolynomial{K}) where K
    if iszero(p1) || iszero(p2)
        return zero(OneVarPolynomial{K})
    end
    coefficients = zeros(K, length(p1.coefficients) + length(p2.coefficients) - 1)
    for (j1, coef1) in enumerate(p1.coefficients)
        for (j2, coef2) in enumerate(p2.coefficients)
            coefficients[j1 + j2 - 1] += coef1 * coef2
        end
    end
    return OneVarPolynomial{K}(trim_zeros!(coefficients))
end

function *(z::Union{K, Number}, p::OneVarPolynomial{K}) where K
    z = convert(K, z)
    if iszero(z)
        return zero(OneVarPolynomial{K})
    else
        return OneVarPolynomial{K}([z * c for c in p.coefficients])
    end
end

*(p::OneVarPolynomial{K}, z::Union{K, Number}) where K = convert(K, z) * p
/(p::OneVarPolynomial{K}, z::Union{K, Number}) where K = inv(convert(K, z)) * p

function ^(p::OneVarPolynomial{K}, num::Integer) where K
    @assert num >= 0
    output = one(p)
    while true
        if isodd(num)
            output *= p
        end
        num = fld(num, 2)
        if num > 0
            p *= p
        else
            return output
        end
    end
end

degree(p::OneVarPolynomial) = length(p.coefficients) - 1

function leading_coefficient(p::OneVarPolynomial{K}) where K
    if degree(p) < 0
        return zero(K)
    else
        return p.coefficients[end]
    end
end

function get_coefficient(p::OneVarPolynomial{K}, d::Integer) where K
    @assert d >= 0
    if d > degree(p)
        return zero(K)
    else
        return p.coefficients[d + 1]
    end
end

get_coefficient(d::Integer) = p -> get_coefficient(p, d)

# Evaluate polynomial at a given value x.
function get_value(p::OneVarPolynomial{K}, x::K) where K
    if iszero(x)
        return get_coefficient(p, 0)
    end
    output = zero(K)
    for coef in Iterators.reverse(p.coefficients)
        output = x * output + coef
    end
    return output
end

get_value(x::K) where K = p -> get_value(p, x)

function generator_power(::Type{K}, n::Integer) where K
    @assert n >= 0
    coefficients = zeros(K, n + 1)
    coefficients[end] = one(K)
    return OneVarPolynomial{K}(coefficients)
end

#=
Promotion and conversion
=#

function promote_rule(::Type{OneVarPolynomial{K}}, ::Type{OneVarPolynomial{L}}) where {K, L}
    return OneVarPolynomial{promote_type(K, L)}
end
promote_rule(::Type{OneVarPolynomial{K}}, ::Type{K}) where K = OneVarPolynomial{K}

convert(::Type{OneVarPolynomial{K}}, p::OneVarPolynomial{K}) where K = p
function convert(::Type{OneVarPolynomial{K}}, p::OneVarPolynomial{L}) where {K, L}
    return OneVarPolynomial{K}(trim_zeros!(convert.(K, p.coefficients)))
end
convert(::Type{OneVarPolynomial{K}}, z::Union{K, Number}) where K = OneVarPolynomial{K}(convert(K, z))

==(p1::OneVarPolynomial, p2::OneVarPolynomial) = ==(promote(p1, p2)...)
==(p::OneVarPolynomial{K}, z::Union{K, Number}) where K = ==(promote(p, convert(K, z))...)
==(z::Union{K, Number}, p::OneVarPolynomial{K}) where K = ==(promote(convert(K, z), p)...)

+(p1::OneVarPolynomial, p2::OneVarPolynomial) = +(promote(p1, p2)...)
+(p::OneVarPolynomial{K}, z::Union{K, Number}) where K = +(promote(p, convert(K, z))...)
+(z::Union{K, Number}, p::OneVarPolynomial{K}) where K = +(promote(convert(K, z), p)...)

*(p1::OneVarPolynomial, p2::OneVarPolynomial) = *(promote(p1, p2)...)

-(p::OneVarPolynomial{K}) where K = (-one(K)) * p
-(p1::OneVarPolynomial, p2::OneVarPolynomial) = p1 + (-p2)
-(p::OneVarPolynomial{K}, z::Union{K, Number}) where K = p + (-convert(K, z))
-(z::Union{K, Number}, p::OneVarPolynomial{K}) where K = z + (-p)
