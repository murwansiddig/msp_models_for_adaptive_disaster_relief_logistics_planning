@time objs_RH2SSP, RH2SSP_bar, RH2SSP_low, RH2SSP_high, elapsed_RH2SSP = rolling_horizon();


fname = "./output/benchmark/roliing2SPresults.csv"
df = CSV.read(fname,DataFrame);

results_RH2SSP = Matrix(df);
results_RH2SSP[inst,1] = 0
results_RH2SSP[inst,2] = RH2SSP_bar
results_RH2SSP[inst,3] = RH2SSP_bar-RH2SSP_low
results_RH2SSP[inst,4] = 0
results_RH2SSP[inst,5] = elapsed_RH2SSP
results_RH2SSP[inst,6] = 0

updf = DataFrame(results_RH2SSP, :auto);
CSV.write(fname,updf)


println("############################################################")
println("############################################################")