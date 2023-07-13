df_ap4 = CSV.File("Primary/Ap4 sum.csv") |> DataFrame
df_ap4_sec = CSV.File("Secondary/Ap4 sum.csv") |> DataFrame
df_ap4_static = CSV.File("Static/Ap4 static sum.csv") |> DataFrame

df_ap5 = CSV.File("Primary/Ap5 sum.csv") |> DataFrame
df_ap5_sec = CSV.File("Secondary/Ap5 sum.csv") |> DataFrame
df_ap5_static = CSV.File("Static/Ap5 static sum.csv") |> DataFrame

df_ap6 = CSV.File("Primary/Ap6 sum.csv") |> DataFrame
df_ap6_sec = CSV.File("Secondary/Ap6 sum.csv") |> DataFrame
df_ap6_static = CSV.File("Static/Ap6 static sum.csv") |> DataFrame

df_ap7 = CSV.File("Primary/Ap7 sum.csv") |> DataFrame
df_ap7_sec = CSV.File("Secondary/Ap7 sum.csv") |> DataFrame
df_ap7_static = CSV.File("Static/Ap7 static sum.csv") |> DataFrame

df_ap8 = CSV.File("Primary/Ap8 sum.csv") |> DataFrame
df_ap8_sec = CSV.File("Secondary/Ap8 sum.csv") |> DataFrame
df_ap8_static = CSV.File("Static/Ap8 static sum.csv") |> DataFrame

df_depth = CSV.File("Primary/Depth sum.csv") |> DataFrame
df_depth_sec = CSV.File("Secondary/Depth sum.csv") |> DataFrame

df_depth_disc = CSV.File("Primary/Depth discrete sum.csv") |> DataFrame
df_depth_disc_sec = CSV.File("Secondary/Depth discrete sum.csv") |> DataFrame

dp_list = [df_depth[50000i+1:50000i+50000, :] for i in 0:50]
dp_list_sec = [df_depth_sec[50000i+1:50000i+50000, :] for i in 0:50]

dp_disc_list = [df_depth_disc[50000i+1:50000i+50000, :] for i in 0:50]
dp_disc_list_sec = [df_depth_disc_sec[50000i+1:50000i+50000, :] for i in 0:50]

# Isotpic anomalies and neutron fluence plotted as functions of depths

disp_d(df_ap4, df_ap4_sec, df_ap4_static)
disp_d(df_ap5, df_ap5_sec, df_ap5_static)
disp_d(df_ap6, df_ap6_sec, df_ap6_static)
disp_d(df_ap7, df_ap7_sec, df_ap7_static)
disp_d(df_ap8, df_ap8_sec, df_ap8_static)

disp_j(df_ap4, df_ap4_sec, df_ap4_static)
disp_j(df_ap5, df_ap5_sec, df_ap5_static)
disp_j(df_ap6, df_ap6_sec, df_ap6_static)
disp_j(df_ap7, df_ap7_sec, df_ap7_static)
disp_j(df_ap8, df_ap8_sec, df_ap8_static)

plotcomp(x, y, xlab, ylab) = plotcomp(df_ap4, df_ap5, df_ap6, df_ap7, df_ap8, df_ap4_sec, df_ap5_sec, df_ap6_sec, df_ap7_sec, df_ap8_sec, df_ap4_static, df_ap5_static, df_ap6_static, df_ap7_static, df_ap8_static, x, y, xlab, ylab)

p = (
    plot(plotcomp("Depth_3", "dSm_3", "Depth", "δ150Sm"),xlim=(0,4)), 
    plot(plotcomp("logDepth_3", "dSm_3", "lg(Depth)", "δ150Sm"), legend=false),
    plot(plotcomp("Depth_3", "dGd_3", "Depth", "δ158Gd"),xlim=(0,4), legend=false), 
    plot(plotcomp("logDepth_3", "dGd_3", "lg(Depth)", "δ158Gd"), legend=false),
    plot(plotcomp("Depth_3", "jSm_3", "Depth", "J149Sm"),xlim=(0,4), legend=false), 
    plot(plotcomp("logDepth_3", "jSm_3", "lg(Depth)", "J149Sm"), legend=false),
    plot(plotcomp("Depth_3", "jGd_3", "Depth", "J157Gd"),xlim=(0,4), legend=false), 
    plot(plotcomp("logDepth_3", "jGd_3", "lg(Depth)", "J157Gd"), legend=false)
)

plot(p..., layout=(4, 2), size=(1200, 1600), legendfontsize=10, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=-2Plots.mm)

lwd = 4

p1=plot(xlabel="Depth (m)", ylabel=" ", framestyle=:box, xguidefontsize=28, yguidefontsize=28, tickfontsize=24, legendfontsize=24, legend=:bottomleft, bottom_margin=15Plots.mm, size=(1500, 800))
plot!(df_depth_sec[!, :Depth_1], seriestype=:hist, bin=0.05:0.1:120, label="0.3 Gy", linewidth=lwd, color=rgbset[2])
plot!(df_depth_sec[!, :Depth_2], seriestype=:hist, bin=0.05:0.1:120, label="0.6 Gy", linewidth=lwd, color=rgbset[3])
plot!(df_depth_sec[!, :Depth_3], seriestype=:hist, bin=0.05:0.1:120, label="0.9 Gy", linewidth=lwd, color=rgbset[4])
plot!(df_depth_sec[!, :Depth_4], seriestype=:hist, bin=0.05:0.1:120, label="1.2 Gy", linewidth=lwd, color=rgbset[5])
plot!(df_depth_sec[!, :Depth_5], seriestype=:hist, bin=0.05:0.1:120, label="1.5 Gy", linewidth=lwd, color=rgbset[6])
plot!(df_depth_sec[!, :Depth_6], seriestype=:hist, bin=0.05:0.1:120, label="1.8 Gy", linewidth=lwd, color=rgbset[7])

p2=plot(xlabel="Depth (m)", ylabel=" ", framestyle=:box, xguidefontsize=28, yguidefontsize=28, tickfontsize=24, legendfontsize=24, legend=:bottomleft, bottom_margin=15Plots.mm, size=(1500, 800))
plot!(df_depth[!, :Depth_1], seriestype=:hist, bin=0.05:0.1:120, label="0.3 Gy", linewidth=lwd, color=rgbset[2])
plot!(df_depth[!, :Depth_2], seriestype=:hist, bin=0.05:0.1:120, label="0.6 Gy", linewidth=lwd, color=rgbset[3])
plot!(df_depth[!, :Depth_3], seriestype=:hist, bin=0.05:0.1:120, label="0.9 Gy", linewidth=lwd, color=rgbset[4])
plot!(df_depth[!, :Depth_4], seriestype=:hist, bin=0.05:0.1:120, label="1.2 Gy", linewidth=lwd, color=rgbset[5])
plot!(df_depth[!, :Depth_5], seriestype=:hist, bin=0.05:0.1:120, label="1.5 Gy", linewidth=lwd, color=rgbset[6])
plot!(df_depth[!, :Depth_6], seriestype=:hist, bin=0.05:0.1:120, label="1.8 Gy", linewidth=lwd, color=rgbset[7])

plot(p1,p2, layout=(2, 1), size=(1600, 1200), legendfontsize=16, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=12Plots.mm)

# Relative depth varied as time

grpnb = 21
lwd = 0.0

p1=plot(xlabel="Relative depth (m)", framestyle=:box, xguidefont=22, tickfontsize=18, legendfontsize=15, xlim=(-4, 4), ylim=(0, 8e3), size=(1000, 700), bottom_margin=12Plots.mm)
plot!(dp_list[grpnb][!, :Net_1], seriestype=:hist, bin=-4:0.025:4, label="0.3 Gy", linewidth=lwd, alpha=0.75, color=rgbset[2])
plot!(dp_list[grpnb][!, :Net_2], seriestype=:hist, bin=-4:0.025:4, label="0.4 Gy", linewidth=lwd, alpha=0.75, color=rgbset[3])
plot!(dp_list[grpnb][!, :Net_3], seriestype=:hist, bin=-4:0.025:4, label="0.9 Gy", linewidth=lwd, alpha=0.75, color=rgbset[4])
plot!(dp_list[grpnb][!, :Net_4], seriestype=:hist, bin=-4:0.025:4, label="1.2 Gy", linewidth=lwd, alpha=0.75, color=rgbset[5])
plot!(dp_list[grpnb][!, :Net_5], seriestype=:hist, bin=-4:0.025:4, label="1.5 Gy", linewidth=lwd, alpha=0.75, color=rgbset[6])
plot!(dp_list[grpnb][!, :Net_6], seriestype=:hist, bin=-4:0.025:4, label="1.8 Gy", linewidth=lwd, alpha=0.75, color=rgbset[7])

p2=plot(xlabel="Relative depth (m)", framestyle=:box, xguidefont=22, tickfontsize=18, legendfontsize=15, xlim=(-4, 4), ylim=(0, 2e3), size=(1000, 700), bottom_margin=12Plots.mm)
plot!(dp_list_sec[grpnb][!, :Net_1], seriestype=:hist, bin=-4:0.025:4, label="0.3 Gy", linewidth=lwd, alpha=0.75, color=rgbset[2])
plot!(dp_list_sec[grpnb][!, :Net_2], seriestype=:hist, bin=-4:0.025:4, label="0.4 Gy", linewidth=lwd, alpha=0.75, color=rgbset[3])
plot!(dp_list_sec[grpnb][!, :Net_3], seriestype=:hist, bin=-4:0.025:4, label="0.9 Gy", linewidth=lwd, alpha=0.75, color=rgbset[4])
plot!(dp_list_sec[grpnb][!, :Net_4], seriestype=:hist, bin=-4:0.025:4, label="1.2 Gy", linewidth=lwd, alpha=0.75, color=rgbset[5])
plot!(dp_list_sec[grpnb][!, :Net_5], seriestype=:hist, bin=-4:0.025:4, label="1.5 Gy", linewidth=lwd, alpha=0.75, color=rgbset[6])
plot!(dp_list_sec[grpnb][!, :Net_6], seriestype=:hist, bin=-4:0.025:4, label="1.8 Gy", linewidth=lwd, alpha=0.75, color=rgbset[7])

plot(p1,p2, layout=(1, 2), size=(1600, 750), legendfontsize=16, left_margin=10Plots.mm, right_margin=0.05Plots.mm, bottom_margin=12Plots.mm)

# Depth mixing histgram

lwd = 0.0
pdep = plot(dp_disc_list[1][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="0 m", linewidth=lwd, alpha=0.8, color=rgbset[2])
plot!(dp_disc_list[6][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="2.5 m", linewidth=lwd, alpha=0.8, color=rgbset[3])
plot!(dp_disc_list[11][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="5 m", linewidth=lwd, alpha=0.8, color=rgbset[4])
plot!(dp_disc_list[16][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="7.5 m", linewidth=lwd, alpha=0.8, color=rgbset[5])
plot!(dp_disc_list[21][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="10 m", linewidth=lwd, alpha=0.8, color=rgbset[6])
plot!(dp_disc_list[26][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="12.5 m", linewidth=lwd, alpha=0.8, color=rgbset[7])
plot!(dp_disc_list[31][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="15 m", linewidth=lwd, alpha=0.8, color=rgbset[8])

pdep_abs = plot(dp_disc_list[1][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="0 m", linewidth=lwd, alpha=0.5, color=rgbset[2], yscale=:log10)
plot!(dp_disc_list[6][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="2.5 m", linewidth=lwd, alpha=0.5, color=rgbset[3], yscale=:log10)
plot!(dp_disc_list[11][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="5 m", linewidth=lwd, alpha=0.5, color=rgbset[4], yscale=:log10)
plot!(dp_disc_list[16][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="7.5 m", linewidth=lwd, alpha=0.5, color=rgbset[5], yscale=:log10)
plot!(dp_disc_list[21][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="10 m", linewidth=lwd, alpha=0.5, color=rgbset[6], yscale=:log10)
plot!(dp_disc_list[26][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="12.5 m", linewidth=lwd, alpha=0.5, color=rgbset[7], yscale=:log10)
plot!(dp_disc_list[31][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="15 m", linewidth=lwd, alpha=0.5, color=rgbset[8], yscale=:log10)

plot(pdep, pdep_abs, layout=2)
plot!(xlabel="Depth (m)", framestyle=:box, xguidefontsize=36, tickfontsize=32, legendfontsize=20, size=(2000, 1500), bottom_margin=18Plots.mm)

# Depth mixing histgram

lwd = 0.0
pdep = plot(dp_disc_list[1][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="0 m", linewidth=lwd, alpha=0.8, color=rgbset[2])
plot!(dp_disc_list_sec[6][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="2.5 m", linewidth=lwd, alpha=0.8, color=rgbset[3])
plot!(dp_disc_list_sec[11][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="5 m", linewidth=lwd, alpha=0.8, color=rgbset[4])
plot!(dp_disc_list_sec[16][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="7.5 m", linewidth=lwd, alpha=0.8, color=rgbset[5])
plot!(dp_disc_list_sec[21][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="10 m", linewidth=lwd, alpha=0.8, color=rgbset[6])
plot!(dp_disc_list_sec[26][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="12.5 m", linewidth=lwd, alpha=0.8, color=rgbset[7])
plot!(dp_disc_list_sec[31][!, :Depth_6], seriestype=:hist, bin=0:0.1:20, label="15 m", linewidth=lwd, alpha=0.8, color=rgbset[8])

pdep_abs = plot(dp_disc_list[1][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="0 m", linewidth=lwd, alpha=0.5, color=rgbset[2], yscale=:log10)
plot!(dp_disc_list_sec[6][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="2.5 m", linewidth=lwd, alpha=0.5, color=rgbset[3], yscale=:log10)
plot!(dp_disc_list_sec[11][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="5 m", linewidth=lwd, alpha=0.5, color=rgbset[4], yscale=:log10)
plot!(dp_disc_list_sec[16][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="7.5 m", linewidth=lwd, alpha=0.5, color=rgbset[5], yscale=:log10)
plot!(dp_disc_list_sec[21][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="10 m", linewidth=lwd, alpha=0.5, color=rgbset[6], yscale=:log10)
plot!(dp_disc_list_sec[26][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="12.5 m", linewidth=lwd, alpha=0.5, color=rgbset[7], yscale=:log10)
plot!(dp_disc_list_sec[31][!, :Depth_abs_6], seriestype=:hist, bin=0:0.1:20, label="15 m", linewidth=lwd, alpha=0.5, color=rgbset[8], yscale=:log10)

plot(pdep, pdep_abs, layout=2)
plot!(xlabel="Depth (m)", framestyle=:box, xguidefontsize=36, tickfontsize=32, legendfontsize=20, size=(2000, 1500), bottom_margin=18Plots.mm)

# Apollo fittings with secondary impacts

# 0.3, 0.4, 0.5, 0.6, 0.7, 0.8
df_ap51 = CSV.File("Secondary/Apollo fitting/Ap5 sum 1.csv") |> DataFrame
df_ap71 = CSV.File("Secondary/Apollo fitting/Ap7 sum 1.csv") |> DataFrame

# 0.9, 1.0, 1.1, 1.2, 1.3, 1.4
df_ap52 = CSV.File("Secondary/Apollo fitting/Ap5 sum 2.csv") |> DataFrame
df_ap72 = CSV.File("Secondary/Apollo fitting/Ap7 sum 2.csv") |> DataFrame

df_smp = CSV.File("Appolo samples.csv") |> DataFrame

smp_Sm = groupby(dropmissing(df_smp[!, Not(6, 8, 10)]), ["Type", "References"])
smp_Gd = groupby(dropmissing(df_smp[!, Not(5, 7, 9)]), ["Type", "References"])

mksz = 5
mkwd = 0.5
lwd = 2
leg = ["Hidaka et al. 2000", "Russ et al. 1972",
    "Hidaka et al. 2007", "Russ et al. 1973",
    "Hidaka et al. 2007", "Curtis and Wasserburg 1975"
]

pSm = plot(ylabel="J¹⁴⁹Sm (10¹⁸)", framestyle=:box, xlim=(-0.25, 3), ylim=(0, 0.21), legend=:topleft)
scatter!(smp_Sm[1][!, "Depth (m)"], smp_Sm[1][!, "J149Sm"], label=leg[1], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9], markershape=:rect)
scatter!(smp_Sm[2][!, "Depth (m)"], smp_Sm[2][!, "J149Sm"], label=leg[2], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9])
scatter!(smp_Sm[3][!, "Depth (m)"], smp_Sm[3][!, "J149Sm"], label=leg[3], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10], markershape=:rect)
scatter!(smp_Sm[4][!, "Depth (m)"], smp_Sm[4][!, "J149Sm"], label=leg[4], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10])
scatter!(smp_Sm[5][!, "Depth (m)"], smp_Sm[5][!, "J149Sm"], label=leg[5], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[11], markershape=:rect)
scatter!(smp_Sm[6][!, "Depth (m)"], smp_Sm[6][!, "J149Sm"], label=leg[6], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[11])

pGd = plot(xlabel="Depth (m)", ylabel="J¹⁵⁷Gd (10¹⁸)", framestyle=:box, xlim=(-0.25, 3), ylim=(0, 0.051), legend=:topleft)
scatter!(smp_Gd[1][!, "Depth (m)"], smp_Gd[1][!, "J157Gd"], label=leg[1], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9], markershape=:rect)
scatter!(smp_Gd[2][!, "Depth (m)"], smp_Gd[2][!, "J157Gd"], label=leg[2], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9])
scatter!(smp_Gd[3][!, "Depth (m)"], smp_Gd[3][!, "J157Gd"], label=leg[3], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10], markershape=:rect)
scatter!(smp_Gd[4][!, "Depth (m)"], smp_Gd[4][!, "J157Gd"], label=leg[4], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10])
scatter!(smp_Gd[5][!, "Depth (m)"], smp_Gd[5][!, "J157Gd"], label=leg[5], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[11], markershape=:rect)
scatter!(smp_Gd[6][!, "Depth (m)"], smp_Gd[6][!, "J157Gd"], label=leg[6], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[11])

plot(pSm, pGd, layout=(2, 1), size=(950, 800), xguidefont=20, yguidefont=20, tickfontsize=18, legendfontsize=10, left_margin=6Plots.mm)

pSm15 = plot(ylabel="J¹⁴⁹Sm (10¹⁸)", framestyle=:box, xlim=(-0.25, 4), ylim=(0, 0.21), legend=:topleft)
@df df_ap51 plot!(:Depth_5, :jSm_5, label="0.7 Gy mixed", linewidth=lwd, linecolor=rgbset[4])
@df df_ap51 plot!(:Depth_6, :jSm_6, label="0.8 Gy mixed", linewidth=lwd, linecolor=rgbset[5])
@df df_ap52 plot!(:Depth_1, :jSm_1, label="0.9 Gy mixed", linewidth=lwd, linecolor=rgbset[6])
@df df_ap52 plot!(:Depth_2, :jSm_2, label="1.0 Gy mixed", linewidth=lwd, linecolor=rgbset[7])
@df df_ap52 plot!(:Depth_3, :jSm_3, label="1.1 Gy mixed", linewidth=lwd, linecolor=rgbset[8])
@df df_ap52 plot!(:Depth_4, :jSm_4, label="1.2 Gy mixed", linewidth=lwd, linecolor=rgbset[9])
scatter!(smp_Sm[1][!, "Depth (m)"], smp_Sm[1][!, "J149Sm"], label=leg[1], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9], markershape=:rect)
scatter!(smp_Sm[2][!, "Depth (m)"], smp_Sm[2][!, "J149Sm"], label=leg[2], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9])

pSm16 = plot(ylabel="J¹⁴⁹Sm (10¹⁸)", framestyle=:box, xlim=(-0.25, 4), ylim=(0, 0.26), legend=:topleft)
@df df_ap71 plot!(:Depth_4, :jSm_4, label="0.6 Gy mixed", linewidth=lwd, linecolor=rgbset[4])
@df df_ap71 plot!(:Depth_5, :jSm_5, label="0.7 Gy mixed", linewidth=lwd, linecolor=rgbset[5])
@df df_ap71 plot!(:Depth_6, :jSm_6, label="0.8 Gy mixed", linewidth=lwd, linecolor=rgbset[6])
@df df_ap72 plot!(:Depth_1, :jSm_1, label="0.9 Gy mixed", linewidth=lwd, linecolor=rgbset[7])
@df df_ap72 plot!(:Depth_2, :jSm_2, label="1.0 Gy mixed", linewidth=lwd, linecolor=rgbset[8])
scatter!(smp_Sm[3][!, "Depth (m)"], smp_Sm[3][!, "J149Sm"], label=leg[3], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10], markershape=:rect)
scatter!(smp_Sm[4][!, "Depth (m)"], smp_Sm[4][!, "J149Sm"], label=leg[4], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10])

pGd15 = plot(xlabel="Depth (m)", ylabel="J¹⁵⁷Gd (10¹⁸)", framestyle=:box, xlim=(-0.25, 4), ylim=(0, 0.041), legend=:topleft)
@df df_ap51 plot!(:Depth_2, :jGd_2, label="0.4 Gy mixed", linewidth=lwd, linecolor=rgbset[4])
@df df_ap51 plot!(:Depth_3, :jGd_3, label="0.5 Gy mixed", linewidth=lwd, linecolor=rgbset[5])
@df df_ap51 plot!(:Depth_4, :jGd_4, label="0.6 Gy mixed", linewidth=lwd, linecolor=rgbset[6])
@df df_ap51 plot!(:Depth_5, :jGd_5, label="0.7 Gy mixed", linewidth=lwd, linecolor=rgbset[7])
scatter!(smp_Gd[1][!, "Depth (m)"], smp_Gd[1][!, "J157Gd"], label=leg[1], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9], markershape=:rect)
scatter!(smp_Gd[2][!, "Depth (m)"], smp_Gd[2][!, "J157Gd"], label=leg[2], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[9])

pGd16 = plot(xlabel="Depth (m)", ylabel="J¹⁵⁷Gd (10¹⁸)", framestyle=:box, xlim=(-0.25, 4), ylim=(0, 0.081), legend=:topleft)
@df df_ap71 plot!(:Depth_1, :jGd_1, label="0.3 Gy mixed", linewidth=lwd, linecolor=rgbset[4])
@df df_ap71 plot!(:Depth_2, :jGd_2, label="0.4 Gy mixed", linewidth=lwd, linecolor=rgbset[5])
@df df_ap71 plot!(:Depth_3, :jGd_3, label="0.5 Gy mixed", linewidth=lwd, linecolor=rgbset[6])
scatter!(smp_Gd[3][!, "Depth (m)"], smp_Gd[3][!, "J157Gd"], label=leg[3], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10], markershape=:rect)
scatter!(smp_Gd[4][!, "Depth (m)"], smp_Gd[4][!, "J157Gd"], label=leg[4], markersize=mksz, markerstrokewidth=mkwd, markercolor=rgbset[10])

plot(pSm15, pGd15, pSm16, pGd16, layout=(2, 2), size=(1600, 1200), xguidefont=20, yguidefont=20, tickfontsize=18, legendfontsize=10, left_margin=8Plots.mm, bottom_margin=6Plots.mm)

#=

#Fitting using Cauchy distribution 

cmiu = zeros(6)
csig = zeros(6)

cp1, (cmiu[1], csig[1]) = cauchyplot(sort(df_depth_disc[!, :Net_1])[50000:end-50000])
cp2, (cmiu[2], csig[2]) = cauchyplot(sort(df_depth_disc[!, :Net_2])[50000:end-50000])
cp3, (cmiu[3], csig[3]) = cauchyplot(sort(df_depth_disc[!, :Net_3])[50000:end-50000])
cp4, (cmiu[4], csig[4]) = cauchyplot(sort(df_depth_disc[!, :Net_4])[50000:end-50000])
cp5, (cmiu[5], csig[5]) = cauchyplot(sort(df_depth_disc[!, :Net_5])[50000:end-50000])
cp6, (cmiu[6], csig[6]) = cauchyplot(sort(df_depth_disc[!, :Net_6])[50000:end-50000])

plot(cp1, cp2, cp3, cp4, cp5, cp6, layout=(3, 2), size=(1000, 800))

disp_age = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0]

plot(framestyle=:box)
plot!(x -> 0.2162x - 0.0310, linewidth=3, xlim=(0, 4), label="miu fit")
scatter!((disp_age, cmiu), label="miu")
plot!(x -> 0.1664x - 0.0256, linewidth=3, xlim=(0, 4), label="sigma fit")
scatter!((disp_age, csig), label="sigma")

=#

#=
plot(Beta(0.5, 0.12), xlim=(0, 1), ylim=(0, 10))

rbdy = 1.8
mlist = dp_disc_list[2][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.0005:0.5005)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5005), ylim=(0, 10))
plot!(Arcsine(0, rbdy^2), xlim=(0, 0.5), ylim=(0, 10))

rbdy = 3.4
mlist = dp_disc_list[4][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, rbdy^2), xlim=(0, 0.5), ylim=(0, 6))

rbdy = 7
mlist = dp_disc_list[8][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, rbdy^2), xlim=(0, 0.5), ylim=(0, 3))

rbdy = 12
mlist = dp_disc_list[12][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, rbdy^2), xlim=(0, 0.5), ylim=(0, 1))

rbdy = 18
mlist = dp_disc_list[16][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, rbdy^2), xlim=(0, 0.5), ylim=(0, 1))

scatter([0.5, 1.5, 3.5, 5.5, 7.5], [1.8, 3.4, 7, 12, 18])

rbdy = 0.3
mlist = dp_disc_list[2][!, :Depth_abs_6]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 10))

rbdy = 0.2
mlist = dp_disc_list[2][!, :Depth_abs_5]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 1))

rbdy = 0.12
mlist = dp_disc_list[2][!, :Depth_abs_4]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 1))

rbdy = 0.05
mlist = dp_disc_list[2][!, :Depth_abs_3]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 1))

rbdy = 0.02
mlist = dp_disc_list[2][!, :Depth_abs_2]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 1))

rbdy = 0.004
mlist = dp_disc_list[2][!, :Depth_abs_1]
hm = histogram(mlist, bin=0.001:0.001:0.5)[1][1]
mx, my = hm[:x][1:6:end], hm[:y][1:6:end] / (0.001 * length(mlist))
plot((mx, my), label="Sim", xlim=(0, 0.5))
plot!(Arcsine(0, 1 / rbdy), xlim=(0, 0.5), ylim=(0, 0.5))

scatter(disp_age, [0.004, 0.02, 0.05, 0.12, 0.2, 0.3])

=#
