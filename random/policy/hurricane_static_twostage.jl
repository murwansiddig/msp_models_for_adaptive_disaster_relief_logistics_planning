#max_iter = 100;
s = 1
t_roll = 1
nbstages1 = T
nbstages2 = nbstages1-1
x_init = x_0
ξ = zeros(Nj)
#define the the model.
master, x, f, y1, z1, v1, θ, subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons, feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons = RH_2SSP_define_models(t_roll,nbstages1,nbstages2,x_init,ξ)

#solve the model.
LB_st2SSP, UB_st2SSP, iter_st2SSP, xval_st2SSP, fval_st2SSP, θval_st2SSP, Elapsed_st2SSP = RH_2SSP_solve_roll(s,t_roll,nbstages1,nbstages2,master,subproblem,x,f,θ,y1,z1,v1,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons);

#eval the model.
objs_st2SSP, st2SSP_bar, st2SSP_low, st2SSP_high, elapsed_st2SSP = static_2SSP_evaluation();

fname = "./output/benchmark/static2SPresults.csv"
df = CSV.read(fname,DataFrame);

results_st2SSP = Matrix(df);
results_st2SSP[inst,1] = LB_st2SSP[end]
results_st2SSP[inst,2] = st2SSP_bar
results_st2SSP[inst,3] = st2SSP_bar-st2SSP_low
results_st2SSP[inst,4] = Elapsed_st2SSP
results_st2SSP[inst,5] = elapsed_st2SSP
results_st2SSP[inst,6] = iter_st2SSP

updf = DataFrame(results_st2SSP, :auto);
CSV.write(fname,updf)
