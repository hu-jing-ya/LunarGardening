module MoonNeutronCapture

using QuadGK, LsqFit
using Random, Statistics, Distributions
using CSV, DataFrames

#= Irradiation =#

const nA = 6.022e23
const yrs = 365.25 * 86400
const bilyrs = 1.0e9 * yrs

const mavgSm = 150.25
const aSm = Int64[144, 147, 148, 149, 150, 152, 154]
const nmSm = ["144Sm", "147Sm", "148Sm", "149Sm", "150Sm", "152Sm", "154Sm"]
const mSm = Float64[143.911999, 146.9148979, 147.9148227, 148.9171847, 149.9172755, 151.9197324, 153.9222093]
const abunSm0 = Float64[0.0307, 0.1499, 0.1124, 0.1382, 0.0738, 0.2675, 0.2275]
const r152Sm0 = abunSm0 ./ abunSm0[6]
const r150_149Sm0 = abunSm0[5] / abunSm0[4]
const r150_152Sm0 = abunSm0[5] / abunSm0[6]

const mavgGd = 157.5202
const aGd = Int64[152, 154, 155, 156, 157, 158, 160]
const nmGd = ["152Gd", "154Gd", "155Gd", "156Gd", "157Gd", "158Gd", "160Gd"]
const mGd = Float64[151.919791, 153.9208656, 154.922622, 155.922622, 156.92396, 157.9244039, 159.9270541]
const abunGd0 = Float64[0.00202914, 0.0218015, 0.148011, 0.204715, 0.156511, 0.248318, 0.218616]
const r160Gd0 = abunGd0 ./ abunGd0[7]
const r158_157Gd0 = abunGd0[6] / abunGd0[5]
const r158_160Gd0 = abunGd0[6] / abunGd0[7]

const normPlanet = 2.0

const sigmaGd157 = 255059 / 1e6# *1e-18 m2
const sigmaSm149 = 40697.5 / 1e6# *1e-18 m2

function init(input)
    rownb = Int64[]
    for i in eachindex(input[!, 1])
        input[i, 1] == "*" && push!(rownb, i)
    end

    rnglist = [rownb[i]+1:rownb[i+1]-1 for i in 1:length(rownb)-1]
    input_list = DataFrame[]
    for rng in rnglist
        push!(input_list, parse.(Float64, input[rng, :]))
    end

    for df in input_list
        df."depth" = (df."from" .+ df."to") / 200.0
        insert!.(eachcol(df), 1, zeros(Float64, 10))
    end

    return input_list
end

fb(x, p) = x * (p[1] + x * (p[2] + x * (p[3] + x * (p[4] + x * (p[5] + x * p[6])))))
modelb(x, p) = x .* (p[1] .+ x .* (p[2] .+ x .* (p[3] .+ x .* (p[4] .+ x .* (p[5] .+ x * p[6])))))
fit_body(depth, val, p0) = curve_fit(modelb, depth, val, p0).param

ft(x, p) = p[1] * exp(p[2] * x)
modelt(x, p) = p[1] * exp.(p[2] * x)
function fit_tail(depth, val, p0, adj)
    st = round(Int64, length(depth) * adj)
    curve_fit(modelt, depth[st:end-1], val[st:end-1], p0).param
end

function f(x, npr, rng)
    x < rng && return fb(x, npr.body)

    return ft(x, npr.tail)
end

struct NeutronProductionRate
    body::Vector{Float64}
    tail::Vector{Float64}
end

NeutronProductionRate(depth, val, pb0, pt0, adj) = NeutronProductionRate(fit_body(depth, val, pb0), fit_tail(depth, val, pt0, adj))
npr(depth, val) = NeutronProductionRate(depth, val, ones(Float64, 6), ones(Float64, 2), 0.8)

function npr(df::DataFrame)
    nprlist = NeutronProductionRate[]
    for col in 3:9
        push!(nprlist, npr(df[!, "depth"], df[!, col]))
    end

    return nprlist
end

function renorm!(npr::NeutronProductionRate, factor)
    for i in eachindex(npr.body)
        npr.body[i] *= factor
    end
    npr.tail[1] *= factor

    return nothing
end

function genatom(npr, depth, rng)
    depth > 15.0 && return zero(depth)
    depth > rng && return ft(depth, npr.tail)
    return fb(depth, npr.body)
end

new(npr, intrvl, depth, rng) = map(x -> genatom(x, depth, rng) * intrvl * bilyrs, npr)

function new!(new_abun, npr, intrvl, depth, rng)
    for i in eachindex(new_abun)
        new_abun[i] = genatom(npr[i], depth, rng) * intrvl * bilyrs
    end

    return nothing
end

function addSm!(dest, old, new)
    for i in eachindex(new)
        new[i] *= old[i]
        dest[i] = old[i] - new[i]
    end

    dest[3] += new[2]
    dest[4] += new[3]
    dest[5] += new[4]

    return nothing
end

function addGd!(dest, old, new)
    for i in eachindex(new)
        new[i] *= old[i]
        dest[i] = old[i] - new[i]
    end

    dest[3] += new[2]
    dest[4] += new[3]
    dest[5] += new[4]
    dest[6] += new[5]

    return nothing
end

function norm!(arr)
    ttl = sum(arr)
    for i in eachindex(arr)
        arr[i] /= ttl
    end

    return nothing
end

function comp(hist, npr, abun0, rng, add!)
    intrvlnb = length(hist[1]) - 1
    abun = [zeros(Float64, 7) for _ in 1:intrvlnb+1]
    new_abun = zeros(Float64, 7)
    copy!(abun[1], abun0)

    t, h = hist[1], hist[2]

    for j in 1:intrvlnb
        new!(new_abun, npr, t[j+1] - t[j], h[j], rng)
        add!(abun[j+1], abun[j], new_abun)
        norm!(abun[j+1])
    end

    return abun
end

compSm(hist, npr, abun0, rng) = comp(hist, npr, abun0, rng, addSm!)
compGd(hist, npr, abun0, rng) = comp(hist, npr, abun0, rng, addGd!)

function fast_comp(hist, npr, abun0, rng, ind, add!)
    indnb = length(ind)
    ind_abun = [zeros(Float64, 7) for _ in 1:indnb]

    prev_abun, next_abun, new_abun = deepcopy(abun0), zeros(Float64, 7), zeros(Float64, 7)

    t, h = hist[1], hist[2]

    for i in 2:ind[1]
        new!(new_abun, npr, t[i] - t[i-1], h[i-1], rng)
        add!(next_abun, prev_abun, new_abun)
        norm!(next_abun)
        copy!(prev_abun, next_abun)
    end
    copy!(ind_abun[1], prev_abun)

    for j in 2:indnb
        for i in ind[j-1]+1:ind[j]
            new!(new_abun, npr, t[i] - t[i-1], h[i-1], rng)
            add!(next_abun, prev_abun, new_abun)
            norm!(next_abun)
            copy!(prev_abun, next_abun)
        end
        copy!(ind_abun[j], prev_abun)
    end

    return ind_abun
end

fast_compSm(hist, npr, abun0, rng, ind) = fast_comp(hist, npr, abun0, rng, ind, addSm!)
fast_compGd(hist, npr, abun0, rng, ind) = fast_comp(hist, npr, abun0, rng, ind, addGd!)

function convert_r!(comp, nb)
    denum = zero(Float64)
    for j in eachindex(comp)
        denum = comp[j][nb]
        for i in eachindex(comp[j])
            comp[j][i] /= denum
        end
    end

    return nothing
end

convert_rSm!(comp) = convert_r!(comp, 6)
convert_rGd!(comp) = convert_r!(comp, 7)

function convert_d!(comp, r0)
    for j in eachindex(comp)
        for i in eachindex(r0)
            comp[j][i] = 1000(comp[j][i] / r0[i] - one(Float64))
        end
    end

    return nothing
end

convert_dSm!(r) = convert_d!(r, r152Sm0)
convert_dGd!(r) = convert_d!(r, r160Gd0)

eps(r, r0) = (r - r0) / (one(Float64) + r)
epsSm(r) = eps(r, r150_149Sm0)
epsGd(r) = eps(r, r158_157Gd0)

delta(r, r0) = 1000((r / r0) - one(Float64))
deltaSm(r) = delta(r, r150_152Sm0)
deltaGd(r) = delta(r, r158_160Gd0)

function convert_df(table, nm)
    npr_df = DataFrame(Isotope=repeat(nm, Int(size(table, 1) / length(nm))))

    npr_df[!, :b1] = table[:, 1]
    npr_df[!, :b2] = table[:, 2]
    npr_df[!, :b3] = table[:, 3]
    npr_df[!, :b4] = table[:, 4]
    npr_df[!, :b5] = table[:, 5]
    npr_df[!, :b6] = table[:, 6]
    npr_df[!, :f1] = table[:, 7]
    npr_df[!, :f2] = table[:, 8]

    return npr_df
end

function read_npr_tab(table)
    npr_all = reshape(mapslices(row -> NeutronProductionRate(row[1:6], row[7:8]), table; dims=2), 7, Int(size(table, 1) / 7))

    return [npr_all[:, i] for i in axes(npr_all, 2)]
end

#= Initializing parameters for fitting neutron production rate

sm_input = CSV.File("Sm input.csv") |> DataFrame
gd_input = CSV.File("Gd input.csv") |> DataFrame

sm_input_list = init(sm_input)
gd_input_list = init(gd_input)

sm_npr_list = npr.(sm_input_list)
gd_npr_list = npr.(gd_input_list)

sm_table = hcat([vcat(sm_npr_list[i][j].body, sm_npr_list[i][j].tail) for i in eachindex(sm_npr_list) for j in eachindex(sm_npr_list[1])]...) |> transpose
gd_table = hcat([vcat(gd_npr_list[i][j].body, gd_npr_list[i][j].tail) for i in eachindex(gd_npr_list) for j in eachindex(gd_npr_list[1])]...) |> transpose

sm_npr_df = convert_df(sm_table, nmSm)
gd_npr_df = convert_df(gd_table, nmGd)

CSV.write("Sm NPR.csv", sm_npr_df)
CSV.write("Gd NPR.csv", gd_npr_df)

=#

sm_table = Float64.((CSV.File("Sm NPR.csv")|>Tables.matrix)[:, 2:9])
gd_table = Float64.((CSV.File("Gd NPR.csv")|>Tables.matrix)[:, 2:9])

sm_npr_list = read_npr_tab(sm_table)
gd_npr_list = read_npr_tab(gd_table)

foreach(npr -> renorm!.(npr, mSm * normPlanet / nA), sm_npr_list)
foreach(npr -> renorm!.(npr, mGd * normPlanet / nA), gd_npr_list)

#= Gardening process =#

const niu = 0.4
const miu = 0.41
const k1 = 1.03
const y = 0.01e6
const rhom = 2500.0
const rhot = 1500.0

const k = 0.3
const c1 = 0.55

const rMoon = 1737400.0
const g = 1.62
const vescMoon = sqrt(2g * rMoon)

const vesc = vescMoon

# Considering primary impacts only
#const b0 = 2.59
#const a0 = 7.75e-8

# Considering secondary impacts
const b0 = 2.97
const a0 = 1.67e-7

const phi = 10.0

# Max lunar impact velocity in m
const vmax = 40000.0

# Lunar escape velocity in m
vnMoon() = 1000(40 - sqrt(1000rand(Float64))) * sqrt(rand(Float64))
randvn = vnMoon

const t1 = pi * a0 * b0 * phi^2 * k1^2 / 4
const t2 = (g / 2) * (rhom / rhot)^(2niu / miu)
const t3 = (rhom / rhot)^(niu * (miu + 2) / miu) * (y / rhot)^((miu + 2) / 2)
const t4 = -miu / (miu + 2)

dc(v, dm) = k1 * dm * ((t2 * dm + t3 / v^miu) / v^2)^t4

dcmax(dm) = dc(vmax, dm)

dlmbd(v, dm) = t1 * ((t2 * dm + t3 / v^miu) / v^2)^2t4 * dm^(one(Float64) - b0)
lmbd(v, d1, d2) = quadgk(dm -> dlmbd(v, dm), d1, d2, rtol=1e-5)[1]
lmbd(d1, d2) = lmbd(vmax, d1, d2)

const p1 = 3miu
const p2 = (k / 8) * c1^p1 * (rhom / rhot)^3niu

const q1 = phi / (5phi - 4)
const q2 = one(Float64) / phi - one(Float64)
const q3 = k / 4pi * c1^p1 * (rhom / rhot)^3niu

vole(rc, fHe) = -2pi * fHe * rc^3 * q2

volr(rc, fHe) = pi * fHe^2 * rc^3 / (8 / 15 + 2fHe)

volesc(vn, vesc, dm) = p2 * (vn / vesc)^p1 * dm^3

fheabs(vn, vesc, dm, rc) = q1 * (4 / 15 - (q3 * (vn / vesc)^p1 - 1 / 3) * (dm / rc)^3)

function fhe(vn, vesc, dm, rc)
    fHe = fheabs(vn, vesc, dm, rc)
    fHe < zero(fHe) ? zero(fHe) : fHe
end

zc(r, rc, fHe) = rc * ((4 / 15 + fHe) * (r / rc)^2 - 4 / 15)

zb(r, rc) = (8 / 15) * rc * ((r / rc)^2 - one(Float64))

cover(r, rc, fHe) = fHe * rc * (rc / r)^3

mix(rc, fHe) = rc * rand(Float64) * sqrt((fHe^2 - 16 / 225) * rand(Float64) + 16 / 225)
function mix(rc, fHe, upb)
    dpmx = mix(rc, fHe)

    return (dpmx, dpmx - upb)
end

function eject(rc, fHe)
    m1 = rand(Float64)
    m2 = rc * fHe * (one(Float64) + q2 * rand(Float64))^3

    return ((one(Float64) - m1) * m2, m1 * m2)
end

function impact(z, z0, r, dm, vesc)
    vn = randvn()
    rc = dc(vn, dm) / 2
    fHe = fhe(vn, vesc, dm, rc)

    r > phi * rc && return (-z, -z0)
    r > rc && return (cover(r, rc, fHe) - z, -z0)
    z < zb(r, rc) && return (zc(r, rc, fHe) - z, -z0)

    upb = zc(r, rc, fHe)
    z < upb && return mix(rc, fHe, upb)

    ve = vole(rc, fHe)
    vr = volr(rc, fHe)
    vtot = ve + vr + volesc(vn, vesc, dm)

    rnd = rand(Float64)
    rnd < ve / vtot && return eject(rc, fHe)
    rnd < (ve + vr) / vtot && return mix(rc, fHe, upb)

    return (10000.0, 10000.0)
end

impact(z, z0, r, dm) = impact(z, z0, r, dm, vesc)

# Number of time intervals
const binnb = 1200
# Max impactor diameter considered
const dmmax = 100.0
const binprm = 125.0

# Set of impactor diameter range for sampling
const dmbinval = reverse!([dmmax * exp(-n / binprm) for n in 0:binnb])
const dmbinset = [(dmbinval[i], dmbinval[i+1]) for i in 1:length(dmbinval)-1]
const dmset = mean.(dmbinset)
const rmaxset = (phi / 2) .* (dcmax.(dmset))
const scale_factor = map(x -> 1 / 1000lmbd(x...), dmbinset) # scale factor in unit of per billion years
const expdist = Exponential.(scale_factor)

function randt!(t, dist)
    for i in eachindex(t)
        rand!(dist[i], t[i])
    end

    return nothing
end

function accumt!(output, input)
    for i in eachindex(input)
        accumulate!(+, output[i], input[i])
    end

    return nothing
end

function end_of!(dest, itemlist, item_end)
    for i in eachindex(itemlist)
        dest[i] = findfirst(x -> x > item_end, itemlist[i]) - 1
    end

    return nothing
end

function end_of(itemlist, item_end)
    dest = length(itemlist)
    end_of!(dest, itemlist, item_end)

    return dest
end

end_of(df::DataFrame, item_end) = findfirst(x -> x > item_end, df[!, 1]) - 1

rand_rad(r) = r * sqrt(rand(Float64))

function histloop(rlist, dmlist, dp0, impnb)
    dplist = zeros(Float64, impnb + 2)
    dplist_abs = zeros(Float64, impnb + 2)
    dplist[1] = dp0
    dplist_abs[1] = dp0
    for i in 2:impnb+1
        dplist[i], dplist_abs[i] = impact(-dplist[i-1], -dplist_abs[i-1], rlist[i], dmlist[i])
    end

    dplist[impnb+2] = dplist[impnb+1]
    dplist_abs[impnb+2] = dplist_abs[impnb+1]

    return dplist, dplist_abs
end

function histgen!(dp0, t, expdist, randt, trunt, agenb)
    randt!(randt, expdist)
    accumt!(trunt, randt)

    end_of!(agenb, trunt, t)
    impnb = sum(agenb)

    agelist, rlist, dmlist = (zeros(Float64, impnb + 1) for _ = 1:3)

    n = 1
    m = n

    for j in eachindex(agenb)
        iszero(agenb[j]) && continue

        for i in 1:agenb[j]
            m = n + i

            agelist[m] = trunt[j][i]
            dmlist[m] = dmset[j]
            rlist[m] = rand_rad(rmaxset[j])
        end

        n += agenb[j]
    end

    seq = sortperm(agelist)
    agelist_reorder = agelist[seq]
    push!(agelist_reorder, t)

    return (agelist_reorder, histloop(view(rlist, seq), view(dmlist, seq), dp0, impnb)...)
end

function histgen(dp0, t, expdist, simnb)
    randt = [zeros(Float64, simnb) for _ in 1:binnb]
    trunt = [zeros(Float64, simnb) for _ in 1:binnb]
    agenb = zeros(Int64, binnb)

    return histgen!(dp0, t, expdist, randt, trunt, agenb)
end

# the last age searched must be larger than the largest age in the agelist
function findage(agelist, age...)
    ind = zeros(Int64, length(age))
    i = 1
    for j in eachindex(agelist)
        agelist[j] >= age[i] && (ind[i] = j - 1; i = i + 1)
    end

    ind[end] = length(agelist)

    return ind
end

function fillnc!(output, j, hist, ind_depth, ind_depth_abs, abunSm0, abunGd0, nprSm, nprGd, rng, ind)

    ind_compSm = fast_compSm(hist, nprSm, abunSm0, rng, ind)
    ind_compGd = fast_compGd(hist, nprGd, abunGd0, rng, ind)

    eSm = map(comp -> epsSm(comp[5] / comp[4]), ind_compSm)
    eGd = map(comp -> epsGd(comp[6] / comp[5]), ind_compGd)

    dSm = map(comp -> deltaSm(comp[5] / comp[6]), ind_compSm)
    dGd = map(comp -> deltaGd(comp[6] / comp[7]), ind_compGd)

    for i in eachindex(ind)
        n = 6i
        output[n-5, j] = ind_depth[i]
        output[n-4, j] = ind_depth_abs[i]
        output[n-3, j] = eSm[i]
        output[n-2, j] = eGd[i]
        output[n-1, j] = dSm[i]
        output[n, j] = dGd[i]
    end

    return nothing
end

function n_capture(dp, t, abunSm0, abunGd0, nprSm, nprGd, rng, simnb, disp_age)
    prtnb = length(dp)
    smpnb = length(rng)
    indnb = length(disp_age)
    randt = [zeros(Float64, simnb) for _ in 1:binnb]
    trunt = [zeros(Float64, simnb) for _ in 1:binnb]
    agenb = zeros(Int64, binnb)

    output = [zeros(Float64, 6indnb, prtnb) for _ in 1:smpnb]

    for j in eachindex(dp)
        hist = histgen!(dp[j], t, expdist, randt, trunt, agenb)
        ind = findage(hist[1], disp_age...)
        ind_depth = hist[2][ind]
        ind_depth_abs = hist[3][ind]
        for k in eachindex(rng)
            fillnc!(output[k], j, hist, ind_depth, ind_depth_abs, abunSm0, abunGd0, nprSm[k], nprGd[k], rng[k], ind)
        end
    end

    return DataFrame.(transpose.(output), :auto)
end

function n_capture_static(dp, t, abunSm0, abunGd0, nprSm, nprGd, rng, disp_age)
    prtnb = length(dp)
    smpnb = length(rng)
    indnb = length(disp_age)

    output = [zeros(Float64, 6indnb, prtnb) for _ in 1:smpnb]

    for j in eachindex(dp)
        hist = ([(t / 1000.0) * (i - 1) for i in 1:1000], fill(dp[j], 1000), fill(dp[j], 1000))

        ind = findage(hist[1], disp_age...)
        ind_depth = hist[2][ind]
        ind_depth_abs = hist[3][ind]
        for k in eachindex(rng)
            fillnc!(output[k], j, hist, ind_depth, ind_depth_abs, abunSm0, abunGd0, nprSm[k], nprGd[k], rng[k], ind)
        end
    end

    return DataFrame.(transpose.(output), :auto)
end

end