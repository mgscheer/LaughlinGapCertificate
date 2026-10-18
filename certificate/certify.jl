include("library.jl")
include("json_reader.jl")

r_prime = 15
r = 22
num_workers = 3
println("(r', r): ($(r_prime), $(r))")
println("num_workers: $(num_workers)")
results = read_certificate("$(@__DIR__)/certificates_$(r_prime)_$(r).json");
if (results.r_prime != r_prime) || (results.r != r)
    error("incorrect certificate file.")
end
println("read certificate file")

if num_workers > 0
    worker_ids = addprocs(num_workers; exeflags="--threads=1")
    @everywhere worker_ids include("$(@__DIR__)/library.jl")
end
println("nprocs(): $(nprocs())")

# Build SDP as a function of t
ssdp_unscaled = setup_SDP(r_prime, r);

# Specialize to t = 0 and apply a scaling
ssdp_unscaled_0 = map(get_coefficient(0), ssdp_unscaled);
ssdp, objective_scale = apply_preconditioner(ssdp_unscaled_0, results.scales);
csdp = ConicSDP(ssdp);

# Check for exact equality: csdp.R * v == csdp.r
RT = sparse(transpose(csdp.R));
equality_check = fast_equality_check(RT, results.v0, csdp.r)
println("Exact equality: $(equality_check)")

# Check for exact positivity
positivity_check = verify_positivity(csdp, results.v0, results.certificate)
println("Exact positivity: $(positivity_check)")

# Objective value
μ1 = dot(csdp.d, results.v0) / objective_scale

# Tail bound with interval arithmetic
δ = (r + 2)//2^BigInt(r)
ω_tilde_2 = ω_tilde_2_interval(r)
ϵ = int(μ1 - δ) - ω_tilde_2
lower_bound_check = ϵ > int(1//100)
interval_guarantee_check = isguaranteed(ϵ)

# Numerical values
println("Float64(μ1) = $(Float64(μ1))")
println("Float64(δ) = $(Float64(δ))")
println("ω_tilde_2 < $(sup(ω_tilde_2))")
println("ϵ > $(inf(ϵ))")

# Exact results
println("t = 0 SDP and tail bound exactly certify ϵ > 1/100: $(equality_check && positivity_check && lower_bound_check && interval_guarantee_check)")
