include("NeutronCapture.jl")
#include("NeutronCapturePrimary.jl")

using .MoonNeutronCapture
MNC = MoonNeutronCapture

using Tables, CSV, DataFrames
using Plots, StatsPlots
using Statistics, Distributions

function save_nc(df, npr_name, prtnb, dp)
    dpcol = [1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32]
    df_dp = df[1][!, dpcol]
    df_dp[!, "dp0"] = dp
    CSV.write("Depth.csv", df_dp)
    for i in eachindex(df)
        nm = npr_name[i]
        CSV.write("$nm raw.csv", df[i])
        df_age = (df[i][!, 1:6], df[i][!, 7:12], df[i][!, 13:18], df[i][!, 19:24], df[i][!, 25:30], df[i][!, 31:36])
        sort!.(df_age)
        CSV.write("$nm sort.csv", hcat(df_age...)[1:Int64(prtnb * 0.20), :]) # 0.20 here truncate the maximum depth to 20 m
    end

    return nothing
end

function save_nc_static(df, npr_name)
    for i in eachindex(df)
        nm = npr_name[i]
        CSV.write("$nm.csv", df[i])
    end

    return nothing
end

function save_h(df, dp)
    dpcol = [1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32]
    df_dp = df[1][!, dpcol]
    df_dp[!, "dp0"] = dp
    CSV.write("Depth discrete.csv", df_dp)

    return nothing
end

const rgbset = [
    RGB(81 / 255, 81 / 255, 81 / 255),
    RGB(241 / 255, 64 / 255, 64 / 255),
    RGB(26 / 255, 111 / 255, 223 / 255),
    RGB(55 / 255, 173 / 255, 107 / 255),
    RGB(177 / 255, 119 / 255, 222 / 255),
    RGB(0 / 255, 203 / 255, 204 / 255),
    RGB(125 / 255, 78 / 255, 78 / 255),
    RGB(142 / 255, 142 / 255, 0 / 255),
    RGB(251 / 255, 101 / 255, 1 / 255),
    RGB(102 / 255, 153 / 255, 204 / 255),
    RGB(111 / 255, 184 / 255, 2 / 255)
]

function dp_pos(dplist, stp)
    h = stp
    j = 1
    pos = zeros(Int64, ceil(Int64, dplist[end] / stp))
    for i in eachindex(dplist)
        dplist[i] > h && (pos[j] = i - 1; j = j + 1; h = h + stp)
    end

    pos[end] = length(dplist)

    return pos
end

function mean_dp(df, indlist)
    df_avg = combine(df[1:indlist[1], :], All() .=> mean)

    for i in 2:length(indlist)-1
        push!(df_avg, combine(df[indlist[i]:indlist[i+1], :], All() .=> mean)[1, :])
    end

    return df_avg
end

mean_dp(df) = mean_dp(df, dp_pos(df[!, 1], 0.01))

function mean_list(df)
    dflist = [df[!, 1:6], df[!, 1+6:6+6], df[!, 1+12:6+12], df[!, 1+18:6+18], df[!, 1+24:6+24], df[!, 1+30:6+30]]

    return mean_dp.(dflist)
end

function save_mean(df, ln)
    sort!(df)
    df_mean = mean_list(df)

    return hcat(first.(df_mean, ln)...)
end

function add_logdepth!(df)
    df[!, :logDepth_1] = log10.(df[!, :Depth_1])
    df[!, :logDepth_2] = log10.(df[!, :Depth_2])
    df[!, :logDepth_3] = log10.(df[!, :Depth_3])
    df[!, :logDepth_4] = log10.(df[!, :Depth_4])
    df[!, :logDepth_5] = log10.(df[!, :Depth_5])
    df[!, :logDepth_6] = log10.(df[!, :Depth_6])

    return nothing
end

function add_flux!(df)
    df[!, :jSm_1] = df[!, :eSm_1] / MNC.sigmaSm149
    df[!, :jSm_2] = df[!, :eSm_2] / MNC.sigmaSm149
    df[!, :jSm_3] = df[!, :eSm_3] / MNC.sigmaSm149
    df[!, :jSm_4] = df[!, :eSm_4] / MNC.sigmaSm149
    df[!, :jSm_5] = df[!, :eSm_5] / MNC.sigmaSm149
    df[!, :jSm_6] = df[!, :eSm_6] / MNC.sigmaSm149

    df[!, :jGd_1] = df[!, :eGd_1] / MNC.sigmaGd157
    df[!, :jGd_2] = df[!, :eGd_2] / MNC.sigmaGd157
    df[!, :jGd_3] = df[!, :eGd_3] / MNC.sigmaGd157
    df[!, :jGd_4] = df[!, :eGd_4] / MNC.sigmaGd157
    df[!, :jGd_5] = df[!, :eGd_5] / MNC.sigmaGd157
    df[!, :jGd_6] = df[!, :eGd_6] / MNC.sigmaGd157

    df[!, :rflux_1] = df[!, :jGd_1] ./ df[!, :jSm_1]
    df[!, :rflux_2] = df[!, :jGd_2] ./ df[!, :jSm_2]
    df[!, :rflux_3] = df[!, :jGd_3] ./ df[!, :jSm_3]
    df[!, :rflux_4] = df[!, :jGd_4] ./ df[!, :jSm_4]
    df[!, :rflux_5] = df[!, :jGd_5] ./ df[!, :jSm_5]
    df[!, :rflux_6] = df[!, :jGd_6] ./ df[!, :jSm_6]

    return nothing
end

function plotCRE(df, df_static, colx, coly, leg, xlab, ylab)
    lwd = 5
    mksz = 4
    mkwd = 0.25

    p = plot(framestyle=:box, xlabel=xlab, ylabel=ylab, tickfontsize=18, xguidefontsize=24, yguidefontsize=24)

    scatter!(df[!, colx[1]], df[!, coly[1]], label=leg[1], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[2])
    scatter!(df[!, colx[2]], df[!, coly[2]], label=leg[2], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[3])
    scatter!(df[!, colx[3]], df[!, coly[3]], label=leg[3], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[4])
    scatter!(df[!, colx[4]], df[!, coly[4]], label=leg[4], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[5])
    scatter!(df[!, colx[5]], df[!, coly[5]], label=leg[5], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[6])
    scatter!(df[!, colx[6]], df[!, coly[6]], label=leg[6], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[7])

    plot!(df_static[!, colx[1]], df_static[!, coly[1]], label=leg[1], linewidth=lwd, linecolor=rgbset[2])
    plot!(df_static[!, colx[2]], df_static[!, coly[2]], label=leg[2], linewidth=lwd, linecolor=rgbset[3])
    plot!(df_static[!, colx[3]], df_static[!, coly[3]], label=leg[3], linewidth=lwd, linecolor=rgbset[4])
    plot!(df_static[!, colx[4]], df_static[!, coly[4]], label=leg[4], linewidth=lwd, linecolor=rgbset[5])
    plot!(df_static[!, colx[5]], df_static[!, coly[5]], label=leg[5], linewidth=lwd, linecolor=rgbset[6])
    plot!(df_static[!, colx[6]], df_static[!, coly[6]], label=leg[6], linewidth=lwd, linecolor=rgbset[7])

    return p
end

col_leg = ["0.3 Gy", "0.6 Gy", "0.9 Gy", "1.2 Gy", "1.5 Gy", "1.8 Gy"]
col_legoff = [false for _ in 1:6]
col_Depth = ["Depth_1", "Depth_2", "Depth_3", "Depth_4", "Depth_5", "Depth_6"]
col_logDepth = ["logDepth_1", "logDepth_2", "logDepth_3", "logDepth_4", "logDepth_5", "logDepth_6"]
col_dSm = ["dSm_1", "dSm_2", "dSm_3", "dSm_4", "dSm_5", "dSm_6"]
col_dSm = ["dSm_1", "dSm_2", "dSm_3", "dSm_4", "dSm_5", "dSm_6"]
col_dGd = ["dGd_1", "dGd_2", "dGd_3", "dGd_4", "dGd_5", "dGd_6"]
col_jSm = ["jSm_1", "jSm_2", "jSm_3", "jSm_4", "jSm_5", "jSm_6"]
col_jGd = ["jGd_1", "jGd_2", "jGd_3", "jGd_4", "jGd_5", "jGd_6"]

plotdSm(df, df_static, xlab, ylab) = plot(plotCRE(df, df_static, col_Depth, col_dSm, col_leg, xlab, ylab),xlim=(0,8))
plotdGd(df, df_static, xlab, ylab) = plot(plotCRE(df, df_static, col_Depth, col_dGd, col_legoff, xlab, ylab),xlim=(0,8))

plotlogdSm(df, df_static, xlab, ylab) = plotCRE(df, df_static, col_logDepth, col_dSm, col_legoff, xlab, ylab)
plotlogdGd(df, df_static, xlab, ylab) = plotCRE(df, df_static, col_logDepth, col_dGd, col_legoff, xlab, ylab)

plotjSm(df, df_static, xlab, ylab) = plot(plotCRE(df, df_static, col_Depth, col_jSm, col_legoff, xlab, ylab),xlim=(0,8))
plotjGd(df, df_static, xlab, ylab) = plot(plotCRE(df, df_static, col_Depth, col_jGd, col_legoff, xlab, ylab),xlim=(0,8))

plotlogjSm(df, df_static, xlab, ylab) = plotCRE(df, df_static, col_logDepth, col_jSm, col_legoff, xlab, ylab)
plotlogjGd(df, df_static, xlab, ylab) = plotCRE(df, df_static, col_logDepth, col_jGd, col_legoff, xlab, ylab)


function disp_d(df_p, df_s, df_static)

    p = (
        plotdSm(df_p, df_static, "Depth (m)", "δSm"), plotlogdSm(df_p, df_static, "lg(Depth)", "δSm"),
        plotdGd(df_p, df_static, "Depth (m)", "δGd"), plotlogdGd(df_p, df_static, "lg(Depth)", "δGd"),
        plotdSm(df_s, df_static, "Depth (m)", "δSm"), plotlogdSm(df_s, df_static, "lg(Depth)", "δSm"),
        plotdGd(df_s, df_static, "Depth (m)", "δGd"), plotlogdGd(df_s, df_static, "lg(Depth)", "δGd")
    )

    return plot(p..., layout=(4, 2), size=(1200, 1600), legendfontsize=12, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=-2Plots.mm)
end

function disp_j(df_p, df_s, df_static)

    p = (
        plotjSm(df_p, df_static, "Depth (m)", "JSm"), plotlogjSm(df_p, df_static, "lg(Depth)", "JSm"),
        plotjGd(df_p, df_static, "Depth (m)", "JGd"), plotlogjGd(df_p, df_static, "lg(Depth)", "JGd"),
        plotjSm(df_s, df_static, "Depth (m)", "JSm"), plotlogjSm(df_s, df_static, "lg(Depth)", "JSm"),
        plotjGd(df_s, df_static, "Depth (m)", "JGd"), plotlogjGd(df_s, df_static, "lg(Depth)", "JGd")
    )

    return plot(p..., layout=(4, 2), size=(1200, 1600), legendfontsize=12, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=-2Plots.mm)
end

function disp_dj(df, df_static)

    p = (
        plotdSm(df, df_static, " ", "δ150Sm"), plotlogdSm(df, df_static, " ", " "),
        plotdGd(df, df_static, " ", "δ158Gd"), plotlogdGd(df, df_static, " ", " "),
        plotjSm(df, df_static, " ", "J149Sm"), plotlogjSm(df, df_static, " ", " "),
        plotjGd(df, df_static, "Depth (m)", "J157Gd"), plotlogjGd(df, df_static, "lg(Depth)", " ")
    )

    return plot(p..., layout=(4, 2), size=(1200, 1600), legendfontsize=12, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=-2Plots.mm)
end

function plotcomp(ap4p, ap5p, ap6p, ap7p, ap8p, ap4s, ap5s, ap6s, ap7s, ap8s, ap4t, ap5t, ap6t, ap7t, ap8t, x, y, xlab, ylab)
    lwd = 2
    mksz = 1

    p = plot(framestyle=:box, xlabel=xlab, ylabel=ylab, tickfontsize=18, xguidefontsize=24, yguidefontsize=24)

    plot!(ap4p[!, x], ap4p[!, y], label="Lunar A prim", line=(:dash, lwd, rgbset[1]))
    plot!(ap5p[!, x], ap5p[!, y], label="Lunar B prim", line=(:dash, lwd, rgbset[2]))
    plot!(ap6p[!, x], ap6p[!, y], label="Lunar C prim", line=(:dash, lwd, rgbset[3]))
    plot!(ap7p[!, x], ap7p[!, y], label="Lunar D prim", line=(:dash, lwd, rgbset[4]))
    plot!(ap8p[!, x], ap8p[!, y], label="Luanr E prim", line=(:dash, lwd, rgbset[5]))

    plot!(ap4s[!, x], ap4s[!, y], label="Lunar A sec", line=(lwd+2, rgbset[1]))
    plot!(ap5s[!, x], ap5s[!, y], label="Lunar B sec", line=(lwd+2, rgbset[2]))
    plot!(ap6s[!, x], ap6s[!, y], label="Lunar C sec", line=(lwd+2, rgbset[3]))
    plot!(ap7s[!, x], ap7s[!, y], label="Lunar D sec", line=(lwd+2, rgbset[4]))
    plot!(ap8s[!, x], ap8s[!, y], label="Lunar E sec", line=(lwd+2, rgbset[5]))

    plot!(ap4t[!, x], ap4t[!, y], label="Lunar A stat", line=(:dashdot, lwd, rgbset[1]))
    plot!(ap5t[!, x], ap5t[!, y], label="Lunar B stat", line=(:dashdot, lwd, rgbset[2]))
    plot!(ap6t[!, x], ap6t[!, y], label="Lunar C stat", line=(:dashdot, lwd, rgbset[3]))
    plot!(ap7t[!, x], ap7t[!, y], label="Lunar D stat", line=(:dashdot, lwd, rgbset[4]))
    plot!(ap8t[!, x], ap8t[!, y], label="Lunar E stat", line=(:dashdot, lwd, rgbset[5]))

    return p
end

function fitcauchy(x, y)
    ymax = maximum(y)
    l, r = findfirst(t -> t >= ymax / 2, y), findlast(t -> t >= ymax / 2, y)

    return (x[l] + x[r]) / 2, 1 / (pi * ymax)
end

# Skew Cauchy distribution
fsc(x, miu, sigma, lambda) = 1 / (sigma * pi * ((x - miu)^2 / (sigma * (lambda * (x - miu) + 1))^2 + 1))
fsc(x, miu, sigma) = fsc(x, miu, sigma, 0)

function cauchyplot(input)
    mlist = input
    mplot = density(mlist, label="Sim")
    hm = mplot[1][1]
    mx, my = hm[:x], hm[:y]
    miu, sigma = fitcauchy(mx, my)
    plot!(x -> fsc(x, miu, sigma), -10, 10, label="fit")

    return mplot, (miu, sigma)
end