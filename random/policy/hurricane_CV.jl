#clairvoyant solution 

#define the clairvoyant model 
m_cv, x_cv, f_cv, y_cv, z_cv, v_cv, dCons_cv = deterministic_model();

#evaluate the clairvoyant solution
objs_cv, cv_bar, cv_low, cv_high, elapsed_cv = clairvoyant_eval();

sleep(30*inst)

fname = "./output/benchmark/CVresults.csv"
df = CSV.read(fname,DataFrame);
results_cv = Matrix(df);
results_cv[inst,1] = 0
results_cv[inst,2] = cv_bar
results_cv[inst,3] = cv_bar-cv_low
results_cv[inst,4] = cv_bar+cv_low
results_cv[inst,5] = elapsed_cv
results_cv[inst,6] = 0

updf = DataFrame(results_cv, :auto);
CSV.write(fname,updf)