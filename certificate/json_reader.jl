using JSON

function parse_rational(s::AbstractString)
    p = split(s, "/")
    if length(p) != 2
        error("not a valid rational number string: $(s).")
    end
    return parse(BigInt, p[1]) // parse(BigInt, p[2])
end

function matrix_from_rows(rows::AbstractVector)
    return permutedims(reduce(hcat, [parse_rational.(row) for row in rows]))
end

function read_certificate(path::AbstractString)
    d = JSON.parsefile(path)
    certificate = [(parse_rational.(c["vals"]), matrix_from_rows(c["vecs"])) for c in d["certificate"]]
    rs = d["rank_stability"]
    return (r_prime = d["r_prime"], r = d["r"], scales = parse_rational.(d["scales"]), v0 = parse_rational.(d["v0"]),
        certificate = certificate, indep_row_inds = Int.(rs["indep_row_inds"]), indep_col_inds = Int.(rs["indep_col_inds"]),
        invertible_prime = Int(rs["invertible_prime"]))
end
