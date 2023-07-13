# Simulation considering primary impacts only

tot_dp = 100.0 #m
t = 1.8 # Gy
simnb = 400
disp_age = (0.3, 0.6, 0.9, 1.2, 1.5, 1.8)
prtnb = 5_000_000
npr_name = ["ap4", "ap5", "ap6", "ap7", "ap8"]
rng = [2.5, 3.6, 2.65, 2.2, 2.2]
dp = Float64[(tot_dp / prtnb) * i for i in 1:prtnb]

df = MNC.n_capture(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[4:8], MNC.gd_npr_list[4:8], rng, simnb, disp_age)
save_nc(df, npr_name, prtnb, dp)

prtnb = 50_000
dp = Float64[0.5i - 0.5 for i in 1:51 for _ in 1:prtnb]
df_h = MNC.n_capture(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[4:8], MNC.gd_npr_list[4:8], rng, simnb, disp_age)
save_h(df_h, dp)

headap = [:Depth_1, :Depth_abs_1, :eSm_1, :eGd_1, :dSm_1, :dGd_1, :Depth_2, :Depth_abs_2, :eSm_2, :eGd_2, :dSm_2, :dGd_2, :Depth_3, :Depth_abs_3, :eSm_3, :eGd_3, :dSm_3, :dGd_3, :Depth_4, :Depth_abs_4, :eSm_4, :eGd_4, :dSm_4, :dGd_4, :Depth_5, :Depth_abs_5, :eSm_5, :eGd_5, :dSm_5, :dGd_5, :Depth_6, :Depth_abs_6, :eSm_6, :eGd_6, :dSm_6, :dGd_6]

df = CSV.File("Primary/ap4 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap4 mean.csv", df_avg)

df = CSV.File("Primary/ap5 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap5 mean.csv", df_avg)

df = CSV.File("Primary/ap6 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap6 mean.csv", df_avg)

df = CSV.File("Primary/ap7 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap7 mean.csv", df_avg)

df = CSV.File("Primary/ap8 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap8 mean.csv", df_avg)

headdp = [:Depth_1, :Depth_abs_1, :Depth_2, :Depth_abs_2, :Depth_3, :Depth_abs_3, :Depth_4, :Depth_abs_4, :Depth_5, :Depth_abs_5, :Depth_6, :Depth_abs_6, :Depth_0]

dp = CSV.File("Primary/Depth.csv") |> DataFrame

rename!(dp, headdp)

dp[!, :Net_1] = dp[!, :Depth_1] - dp[!, :Depth_0]
dp[!, :Net_abs_1] = dp[!, :Depth_abs_1] - dp[!, :Depth_0]
dp[!, :Net_2] = dp[!, :Depth_2] - dp[!, :Depth_0]
dp[!, :Net_abs_2] = dp[!, :Depth_abs_2] - dp[!, :Depth_0]
dp[!, :Net_3] = dp[!, :Depth_3] - dp[!, :Depth_0]
dp[!, :Net_abs_3] = dp[!, :Depth_abs_3] - dp[!, :Depth_0]
dp[!, :Net_4] = dp[!, :Depth_4] - dp[!, :Depth_0]
dp[!, :Net_abs_4] = dp[!, :Depth_abs_4] - dp[!, :Depth_0]
dp[!, :Net_5] = dp[!, :Depth_5] - dp[!, :Depth_0]
dp[!, :Net_abs_5] = dp[!, :Depth_abs_5] - dp[!, :Depth_0]
dp[!, :Net_6] = dp[!, :Depth_6] - dp[!, :Depth_0]
dp[!, :Net_abs_6] = dp[!, :Depth_abs_6] - dp[!, :Depth_0]

CSV.write("Depth sum.csv", dp)

headdp = [:Depth_1, :Depth_abs_1, :Depth_2, :Depth_abs_2, :Depth_3, :Depth_abs_3, :Depth_4, :Depth_abs_4, :Depth_5, :Depth_abs_5, :Depth_6, :Depth_abs_6, :Depth_0]

dp = CSV.File("Primary/Depth discrete.csv") |> DataFrame

rename!(dp, headdp)

dp[!, :Net_1] = dp[!, :Depth_1] - dp[!, :Depth_0]
dp[!, :Net_abs_1] = dp[!, :Depth_abs_1] - dp[!, :Depth_0]
dp[!, :Net_2] = dp[!, :Depth_2] - dp[!, :Depth_0]
dp[!, :Net_abs_2] = dp[!, :Depth_abs_2] - dp[!, :Depth_0]
dp[!, :Net_3] = dp[!, :Depth_3] - dp[!, :Depth_0]
dp[!, :Net_abs_3] = dp[!, :Depth_abs_3] - dp[!, :Depth_0]
dp[!, :Net_4] = dp[!, :Depth_4] - dp[!, :Depth_0]
dp[!, :Net_abs_4] = dp[!, :Depth_abs_4] - dp[!, :Depth_0]
dp[!, :Net_5] = dp[!, :Depth_5] - dp[!, :Depth_0]
dp[!, :Net_abs_5] = dp[!, :Depth_abs_5] - dp[!, :Depth_0]
dp[!, :Net_6] = dp[!, :Depth_6] - dp[!, :Depth_0]
dp[!, :Net_abs_6] = dp[!, :Depth_abs_6] - dp[!, :Depth_0]

CSV.write("Depth discrete sum.csv", dp)

df_ap4 = CSV.File("Primary/Ap4 mean.csv") |> DataFrame

df_ap5 = CSV.File("Primary/Ap5 mean.csv") |> DataFrame

df_ap6 = CSV.File("Primary/Ap6 mean.csv") |> DataFrame

df_ap7 = CSV.File("Primary/Ap7 mean.csv") |> DataFrame

df_ap8 = CSV.File("Primary/Ap8 mean.csv") |> DataFrame

df = df_ap4
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap4 sum.csv", df)

df = df_ap5
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap5 sum.csv", df)

df = df_ap6
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap6 sum.csv", df)

df = df_ap7
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap7 sum.csv", df)

df = df_ap8
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap8 sum.csv", df)

# Simulation considering primary and secondary impacts

tot_dp = 100.0 #m
t = 1.8 # Gy
simnb = 400
disp_age = (0.3, 0.6, 0.9, 1.2, 1.5, 1.8)
prtnb = 5_000_000
npr_name = ["ap4", "ap5", "ap6", "ap7", "ap8"]
rng = [2.5, 3.6, 2.65, 2.2, 2.2]
dp = Float64[(tot_dp / prtnb) * i for i in 1:prtnb]

df = MNC.n_capture(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[4:8], MNC.gd_npr_list[4:8], rng, simnb, disp_age)
save_nc(df, npr_name, prtnb, dp)

prtnb = 50_000
dp = Float64[0.5i - 0.5 for i in 1:51 for _ in 1:prtnb]
df_h = MNC.n_capture(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[4:8], MNC.gd_npr_list[4:8], rng, simnb, disp_age)
save_h(df_h, dp)

headap = [:Depth_1, :Depth_abs_1, :eSm_1, :eGd_1, :dSm_1, :dGd_1, :Depth_2, :Depth_abs_2, :eSm_2, :eGd_2, :dSm_2, :dGd_2, :Depth_3, :Depth_abs_3, :eSm_3, :eGd_3, :dSm_3, :dGd_3, :Depth_4, :Depth_abs_4, :eSm_4, :eGd_4, :dSm_4, :dGd_4, :Depth_5, :Depth_abs_5, :eSm_5, :eGd_5, :dSm_5, :dGd_5, :Depth_6, :Depth_abs_6, :eSm_6, :eGd_6, :dSm_6, :dGd_6]

df = CSV.File("Secondary/ap4 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap4 mean.csv", df_avg)

df = CSV.File("Secondary/ap5 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap5 mean.csv", df_avg)

df = CSV.File("Secondary/ap6 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap6 mean.csv", df_avg)

df = CSV.File("Secondary/ap7 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap7 mean.csv", df_avg)

df = CSV.File("Secondary/ap8 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap8 mean.csv", df_avg)

headdp = [:Depth_1, :Depth_abs_1, :Depth_2, :Depth_abs_2, :Depth_3, :Depth_abs_3, :Depth_4, :Depth_abs_4, :Depth_5, :Depth_abs_5, :Depth_6, :Depth_abs_6, :Depth_0]

dp = CSV.File("Secondary/Depth.csv") |> DataFrame

rename!(dp, headdp)

dp[!, :Net_1] = dp[!, :Depth_1] - dp[!, :Depth_0]
dp[!, :Net_abs_1] = dp[!, :Depth_abs_1] - dp[!, :Depth_0]
dp[!, :Net_2] = dp[!, :Depth_2] - dp[!, :Depth_0]
dp[!, :Net_abs_2] = dp[!, :Depth_abs_2] - dp[!, :Depth_0]
dp[!, :Net_3] = dp[!, :Depth_3] - dp[!, :Depth_0]
dp[!, :Net_abs_3] = dp[!, :Depth_abs_3] - dp[!, :Depth_0]
dp[!, :Net_4] = dp[!, :Depth_4] - dp[!, :Depth_0]
dp[!, :Net_abs_4] = dp[!, :Depth_abs_4] - dp[!, :Depth_0]
dp[!, :Net_5] = dp[!, :Depth_5] - dp[!, :Depth_0]
dp[!, :Net_abs_5] = dp[!, :Depth_abs_5] - dp[!, :Depth_0]
dp[!, :Net_6] = dp[!, :Depth_6] - dp[!, :Depth_0]
dp[!, :Net_abs_6] = dp[!, :Depth_abs_6] - dp[!, :Depth_0]

CSV.write("Depth sum.csv", dp)

headdp = [:Depth_1, :Depth_abs_1, :Depth_2, :Depth_abs_2, :Depth_3, :Depth_abs_3, :Depth_4, :Depth_abs_4, :Depth_5, :Depth_abs_5, :Depth_6, :Depth_abs_6, :Depth_0]

dp = CSV.File("Secondary/Depth discrete.csv") |> DataFrame

rename!(dp, headdp)

dp[!, :Net_1] = dp[!, :Depth_1] - dp[!, :Depth_0]
dp[!, :Net_abs_1] = dp[!, :Depth_abs_1] - dp[!, :Depth_0]
dp[!, :Net_2] = dp[!, :Depth_2] - dp[!, :Depth_0]
dp[!, :Net_abs_2] = dp[!, :Depth_abs_2] - dp[!, :Depth_0]
dp[!, :Net_3] = dp[!, :Depth_3] - dp[!, :Depth_0]
dp[!, :Net_abs_3] = dp[!, :Depth_abs_3] - dp[!, :Depth_0]
dp[!, :Net_4] = dp[!, :Depth_4] - dp[!, :Depth_0]
dp[!, :Net_abs_4] = dp[!, :Depth_abs_4] - dp[!, :Depth_0]
dp[!, :Net_5] = dp[!, :Depth_5] - dp[!, :Depth_0]
dp[!, :Net_abs_5] = dp[!, :Depth_abs_5] - dp[!, :Depth_0]
dp[!, :Net_6] = dp[!, :Depth_6] - dp[!, :Depth_0]
dp[!, :Net_abs_6] = dp[!, :Depth_abs_6] - dp[!, :Depth_0]

CSV.write("Depth discrete sum.csv", dp)

df_ap4 = CSV.File("Secondary/Ap4 mean.csv") |> DataFrame

df_ap5 = CSV.File("Secondary/Ap5 mean.csv") |> DataFrame

df_ap6 = CSV.File("Secondary/Ap6 mean.csv") |> DataFrame

df_ap7 = CSV.File("Secondary/Ap7 mean.csv") |> DataFrame

df_ap8 = CSV.File("Secondary/Ap8 mean.csv") |> DataFrame

df = df_ap4
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap4 sum.csv", df)

df = df_ap5
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap5 sum.csv", df)

df = df_ap6
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap6 sum.csv", df)

df = df_ap7
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap7 sum.csv", df)

df = df_ap8
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap8 sum.csv", df)

# Simulation of static profiles

tot_dp = 15.0 #m
t = 1.8 # Gy
disp_age = (0.3, 0.6, 0.9, 1.2, 1.5, 1.8)
prtnb = 1_000
npr_name_static = ["ap4 static", "ap5 static", "ap6 static", "ap7 static", "ap8 static"]
rng = [2.5, 3.6, 2.65, 2.2, 2.2]
dp = Float64[(tot_dp/prtnb)*i for i in 1:prtnb]

df_static = MNC.n_capture_static(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[4:8], MNC.gd_npr_list[4:8], rng, disp_age)
save_nc_static(df_static, npr_name_static)

df_ap4_static = CSV.File("Static/Ap4 static.csv") |> DataFrame
df_ap5_static = CSV.File("Static/Ap5 static.csv") |> DataFrame
df_ap6_static = CSV.File("Static/Ap6 static.csv") |> DataFrame
df_ap7_static = CSV.File("Static/Ap7 static.csv") |> DataFrame
df_ap8_static = CSV.File("Static/Ap8 static.csv") |> DataFrame

headap = [:Depth_1, :Depth_abs_1, :eSm_1, :eGd_1, :dSm_1, :dGd_1, :Depth_2, :Depth_abs_2, :eSm_2, :eGd_2, :dSm_2, :dGd_2, :Depth_3, :Depth_abs_3, :eSm_3, :eGd_3, :dSm_3, :dGd_3, :Depth_4, :Depth_abs_4, :eSm_4, :eGd_4, :dSm_4, :dGd_4, :Depth_5, :Depth_abs_5, :eSm_5, :eGd_5, :dSm_5, :dGd_5, :Depth_6, :Depth_abs_6, :eSm_6, :eGd_6, :dSm_6, :dGd_6]

df = df_ap4_static
rename!(df, headap)
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap4 static sum.csv", df)

df = df_ap5_static
rename!(df, headap)
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap5 static sum.csv", df)

df = df_ap6_static
rename!(df, headap)
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap6 static sum.csv", df)

df = df_ap7_static
rename!(df, headap)
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap7 static sum.csv", df)

df = df_ap8_static
rename!(df, headap)
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap8 static sum.csv", df)

# Apollo measurements fitting

tot_dp = 100.0 #m
t = 1.4 # Gy
simnb = 400
#disp_age = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
disp_age = (0.9, 1.0, 1.1, 1.2, 1.3, 1.4)
prtnb = 5_000_000
npr_name = ["ap5", "ap7"]
rng = [3.6, 2.2]
dp = Float64[(tot_dp / prtnb) * i for i in 1:prtnb]

df = MNC.n_capture(dp, t, MNC.abunSm0, MNC.abunGd0, MNC.sm_npr_list[[5, 7]], MNC.gd_npr_list[[5, 7]], rng, simnb, disp_age)
save_nc(df, npr_name, prtnb, dp)

headap = [:Depth_1, :Depth_abs_1, :eSm_1, :eGd_1, :dSm_1, :dGd_1, :Depth_2, :Depth_abs_2, :eSm_2, :eGd_2, :dSm_2, :dGd_2, :Depth_3, :Depth_abs_3, :eSm_3, :eGd_3, :dSm_3, :dGd_3, :Depth_4, :Depth_abs_4, :eSm_4, :eGd_4, :dSm_4, :dGd_4, :Depth_5, :Depth_abs_5, :eSm_5, :eGd_5, :dSm_5, :dGd_5, :Depth_6, :Depth_abs_6, :eSm_6, :eGd_6, :dSm_6, :dGd_6]

df = CSV.File("ap5 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap5 mean.csv", df_avg)
Cosmochem91

df = CSV.File("ap7 sort.csv") |> DataFrame
df_avg = save_mean(df, 1500)
rename!(df_avg, headap)
CSV.write("Ap7 mean.csv", df_avg)

df_ap5 = CSV.File("Ap5 mean.csv") |> DataFrame

df_ap7 = CSV.File("Ap7 mean.csv") |> DataFrame

df = df_ap5
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap5 sum.csv", df)

df = df_ap7
add_logdepth!(df)
add_flux!(df)
CSV.write("Ap7 sum.csv", df)