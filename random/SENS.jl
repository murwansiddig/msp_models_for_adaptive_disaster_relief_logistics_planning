osfname = "./data/OOS"*PARAMS[6]*".csv"
temp = CSV.read(osfname,DataFrame);
CSV.write("./data/OOS.csv",temp);


#Fully adaptive
#Defne the model 
m_fa, x_fa, f_fa, y_fa, z_fa, v_fa, ϴ_fa, dCons_fa, FB1Cons_fa, FB2Cons_fa = define_models();


#train the model
LB, train_time, iter = train_models_offline();

println("***********************************************")
println("***********************************************")
println("finished training!")
println("LB = ", LB[end])
println("training time = ", train_time)
println("number of iterations = ", iter)

#evaluate the model 
start=time();
OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
objs_fa = zeros(nbOS,T);
xval_fa = Array{Any,2}(undef,nbOS,T); fval_fa = Array{Any,2}(undef,nbOS,T);
yval_fa = Array{Any,2}(undef,nbOS,T); zval_fa = Array{Any,2}(undef,nbOS,T); vval_fa = Array{Any,2}(undef,nbOS,T);
procurmnt_amount = zeros(T); 
for s=1:nbOS
    xval = zeros(Ni,T);
    for t=1:T
        #the state is known in the first stage; if not sample a new state k 
        k_t = OS_paths[s,t];
        m = OS_M[s]; # realization from OS path corresponding to layer 2
        if k_t ∉ absorbing_states
            if t > 1
                MSP_fa_update_RHS(k_t,t,xval,m);
            end
            #solve the model
            optimize!(m_fa[t,k_t])

            #check the status 
            status = termination_status(m_fa[t,k_t]);
            if status != MOI.OPTIMAL
                println(" in evaluation")
                println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
                return 
            else
                #collect values
                xval_fa[s,t] = value.(x_fa[t,k_t]); xval[:,t] = xval_fa[s,t]
                fval_fa[s,t] = value.(f_fa[t,k_t]);                 
                yval_fa[s,t] = value.(y_fa[t,k_t]);
                zval_fa[s,t] = value.(z_fa[t,k_t]);
                vval_fa[s,t] = value.(v_fa[t,k_t]);
                objs_fa[s,t] = objective_value(m_fa[t,k_t])- value(ϴ_fa[t,k_t]);

                procurmnt_amount[t] += (sum(fval_fa[s,t][N0,i] for i=1:Ni))/nbOS 
            end
        end
    end        
end

sleep(inst*20)

fname = "./output/sensitivity.csv"
df = CSV.read(fname,DataFrame);
results_fa = Matrix(df);
results_fa[inst,1:T] = procurmnt_amount
results_fa[inst,end] = inst; 
updf = DataFrame(results_fa, :auto);
CSV.write(fname,updf)
