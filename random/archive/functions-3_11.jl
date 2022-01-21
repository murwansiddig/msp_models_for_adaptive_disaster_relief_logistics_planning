##############################################################################################################################
##############################################################################################################################

#clairvoyance model functions 

function deterministic_model()
 
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni,t=1:T] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:T] <= f_cap[i,ii]
               0 <= y[i=1:Ni,j=1:Nj,t=1:T]
               0 <= z[j=1:Nj,t=1:T]
               0 <= v[i=1:Ni,t=1:T]
            end
          );

    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:T)
              +sum(sum(ch[i,t]*x[i,t] for i=1:Ni) for t=1:T)
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t] for t=1:T)
              +sum(sum(sum(ca[i,j,t]*y[i,j,t] for j=1:Nj) for i=1:Ni) for t=1:T)
              +sum(sum(z[j,t] for j=1:Nj)*p for t=1:T)
              +sum(sum(v[i,t] for i=1:Ni)*q for t=1:T)       
               );


    #######################
    #Define the constraints.
    dCons = Dict(); #a dictonary to store all the demand constraint
    for t=1:T
        for i=1:Ni
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           +sum(y[i,j,t] for j=1:Nj)
                           +v[i,t]
                           == x_0[i]
                            );
               @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x_0[i]);
            else
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           +sum(y[i,j,t] for j=1:Nj)   
                           +v[i,t]
                           == x[i,t-1]
                            );
               @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x[i,t-1]);                
            end
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        
    end

    return m, x, f, y, z, v, dCons
end

function clairvoyant_eval()
    start=time();
    
    OS_paths = Matrix(CSV.read("OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    objs_cv = zeros(nbOS);
    procurmnt_cost = zeros(T); transport_cost = zeros(T); inventory_cost = zeros(T); dshortage_cost = zeros(T);
    desalvege_cost = zeros(T);
    
    xval_cv = Array{Any,1}(undef,nbOS); fval_cv = Array{Any,1}(undef,nbOS);
    yval_cv = Array{Any,1}(undef,nbOS); zval_cv = Array{Any,1}(undef,nbOS); vval_cv = Array{Any,1}(undef,nbOS);

    for s=1:nbOS
        #find the period when the hurricane made landfall && intensity != 1
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
        if τ === nothing
            continue 
        else
            k_t = OS_paths[s,τ] # state corresponding to the landfall
            m = OS_M[s]; # realization from OS path corresponding to layer 2
            for j=1:Nj, t=1:T
                if t == τ
                    set_normalized_rhs(dCons_cv[τ,j], SCEN[k_t][j,m]); #set the RHS of demand constraint
                else
                    set_normalized_rhs(dCons_cv[t,j], 0); #set the RHS of demand constraint
                end
            end
  
            #solve the model
            optimize!(m_cv)

            #check the status 
            status = termination_status(m_cv);
            if status != MOI.OPTIMAL
                println(" Model in clairvoyant is ", status)
                exit(0);
            else
                xval_cv[s] = value.(x_cv);
                fval_cv[s] = value.(f_cv);                 
                yval_cv[s] = value.(y_cv);
                zval_cv[s] = value.(z_cv);
                vval_cv[s] = value.(v_cv);
                objs_cv[s] = objective_value(m_cv);
                for t=1:T
                    procurmnt_cost[t] += (sum(fval_cv[s][N0,i,t] for i=1:Ni)*h[t])/nbOS
                    transport_cost[t] += sum(sum(cb[i,ii,t]*fval_cv[s][i,ii,t] for ii=1:Ni) for i=1:N0)/nbOS
                    inventory_cost[t] += sum(ch[i,t]*xval_cv[s][i,t] for i=1:Ni)/nbOS
                    dshortage_cost[t] += (sum(zval_cv[s][j,t] for j=1:Nj)*p)/nbOS;
                    desalvege_cost[t] += (sum(vval_cv[s][i,t] for i=1:Ni)*q)/nbOS;            
                end
            end
        end
    end
    
    #writeCost("procurmnt_cost",procurmnt_cost,"cv")
    #writeCost("transport_cost",transport_cost,"cv")
    #writeCost("inventory_cost",inventory_cost,"cv")
    #writeCost("dshortage_cost",dshortage_cost,"cv")
    #writeCost("desalvege_cost",desalvege_cost,"cv")
    
    cv_bar = mean(objs_cv);
    cv_std = std(objs_cv);
    cv_low = cv_bar-1.96*cv_std/sqrt(nbOS);
    cv_high = cv_bar+1.96*cv_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", cv_bar, "±", [cv_low,cv_high]);
    elapsed = time() - start;
    vals = [xval_cv, fval_cv, yval_cv, zval_cv, vval_cv];

    return objs_cv, cv_bar, cv_low, cv_high, elapsed#, vals
end



function writeCost(cost_type,costs,alg)
    fname = string("./costFactors/",alg,"_cost_factor/",cost_type,"/0.1.csv")
    df = CSV.read(fname,DataFrame);
    reSults = Matrix(df);
    reSults[inst,:] = costs
    updf = DataFrame(reSults, :auto);
    CSV.write(fname,updf)
end


##############################################################################################################################
##############################################################################################################################

#MSP fully adaptive model functions 


function MC_sample(current_state)
    states = collect(1:1:K);
    k = sample(states, Weights(P_joint[current_state,:]));
    return k
end

function MSP_fa_update_RHS(k_t,t,xval,m)
    for i=1:Ni        
        set_normalized_rhs(FB1Cons_fa[t,k_t][i], xval[i,t-1]);
        set_normalized_rhs(FB2Cons_fa[t,k_t][i], xval[i,t-1]);
    end 
    for j=1:Nj
        if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
            set_normalized_rhs(dCons_fa[t,k_t][j], SCEN[k_t][j,m]);
        else
            set_normalized_rhs(dCons_fa[t,k_t][j], 0);
        end
    end
end


function termination_check(iter,relative_gap,LB,start,cutviol_iter)
    flag = 0;
    Elapsed = time() - start;
    if iter > max_iter
        flag = 1
        println("max iteration is reached")
    elseif Elapsed > time_limit
        flag = 2
        println("time limit is reached")
    elseif cutviol_iter > cutviol_maxiter
        flag = 3
        println("cut violation is reached")
    else
        if iter > stall
            relative_gap = (LB[iter]-LB[iter-stall])/max(1e-10,abs(LB[iter-stall]));
            if relative_gap < ϵ
                flag = 4
                println("the LB is not making significant progress");
            end
        end
    end
    return flag, Elapsed
end


function stage_t_state_k_problem(t)
 
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni] <= f_cap[i,ii]
               0 <= y[i=1:Ni,j=1:Nj]
               0 <= z[j=1:Nj]
               0 <= v[i=1:Ni]
               0 <= ϴ
            end
          );

    #######################
    #Define the objective.
    @objective(m,
               Min,
                     sum(sum(cb[i,ii,t]*f[i,ii] for ii=1:Ni) for i=1:N0) 
                    +sum(ch[i,t]*x[i] for i=1:Ni)
                    +sum(f[N0,i] for i=1:Ni)*h[t]
                    +sum(sum(ca[i,j,t]*y[i,j] for j=1:Nj) for i=1:Ni)
                    +sum(z[j] for j=1:Nj)*p
                    +sum(v[i] for i=1:Ni)*q
                    +ϴ
               );
    
    #######################
    #Define the constraints.
    FB1Cons = Dict(); #a dictonary to store flow-balance constraints 1
    FB2Cons = Dict(); #a dictonary to store  flow-balance constraints 2
    dCons = Dict(); #a dictonary to store all the demand constraint

    for i=1:Ni
        if t == 1
            FB1Cons[i] = @constraint(m, 
                                        x[i]
                                        +sum(f[i,j] for j=1:Ni if j != i)
                                        -sum(f[j,i] for j=1:N0 if j != i)
                                        +sum(y[i,j] for j=1:Nj)
                                        +v[i]             
                                        == x_0[i]
                                    );
            FB2Cons[i] = @constraint(m, sum(f[i,j] for j=1:Ni if j != i) <= x_0[i]);
        else
            FB1Cons[i] = @constraint(m, 
                                        x[i]
                                        +sum(f[i,j] for j=1:Ni if j != i)
                                        -sum(f[j,i] for j=1:N0 if j != i)
                                        +sum(y[i,j] for j=1:Nj) 
                                        +v[i]                
                                        == 0
                                    );
            FB2Cons[i] = @constraint(m, sum(f[i,j] for j=1:Ni if j != i) <= 0);  
        end
    end
    for j=1:Nj
       dCons[j] = @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) >= 0);
    end
    
    return m, x, f, y, z, v, ϴ, dCons, FB1Cons, FB2Cons
end

function define_models()
    m = Dict(); x = Dict(); f = Dict(); y = Dict(); z = Dict(); v = Dict(); 
    ϴ = Dict(); dCons = Dict(); FB1Cons = Dict(); FB2Cons = Dict(); 

    for  t=1:T, k=1:K
        if t == 1
            m[t,k_init], x[t,k_init], f[t,k_init], y[t,k_init], z[t,k_init], v[t,k_init],
            ϴ[t,k_init], dCons[t,k_init], FB1Cons[t,k_init], FB2Cons[t,k_init] = stage_t_state_k_problem(t);
        else
            if k in absorbing_states
                continue 
            else
                m[t,k], x[t,k], f[t,k], y[t,k], z[t,k], v[t,k],
                ϴ[t,k], dCons[t,k], FB1Cons[t,k], FB2Cons[t,k] = stage_t_state_k_problem(t);
            end            
        end
    end    
    return m, x, f, y, z, v, ϴ, dCons, FB1Cons, FB2Cons
end

function FOSDDP_forward_pass_oneSP_iteration(lb,xval,thetaval)
    
    k_t = copy(k_init);
    in_sample = [k_t]; #what is the state in the first stage
    for t=1:T
        #the state is known in the first stage; if not sample a new state k 
        if t>1
            #sample a new state k
            k_t = MC_sample(in_sample[t-1]);
            push!(in_sample,k_t)
            # if k_t is absorbing no need to do any computation
            if k_t in absorbing_states
                continue 
            end
            #update the RHS
            MSP_fa_update_RHS(k_t,t,xval,rand(1:M));
        end
               
            
        #solve the model
        optimize!(m_fa[t,k_t])

        #check the status 
        status = termination_status(m_fa[t,k_t]);
        if status != MOI.OPTIMAL
            println(" in Forward Pass")
            println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
            return 
        else
            #collect values
            xval[:,t] = value.(x_fa[t,k_t]);
            thetaval[t] = value(ϴ_fa[t,k_t]);
            if t == 1
                lb = objective_value(m_fa[t,k_t]);
            end
        end
    end
    return xval, thetaval, lb, in_sample
end




function FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample)

    cutviolFlag = 0;
    for t=T:-1:2
        #initialize
        Q = zeros(K); #list for all the optimal values
        pi1 = zeros(K,Ni); #list for all the dual multiplies of the first set of constraints
        pi2 = zeros(K,Ni); #list for all the dual multiplies of the second set of constraints
        
        sample_n = in_sample[t-1]; #the states observed at time t-1

        for k=1:K
            if k in absorbing_states
                Q[k] = 0;
                continue
            elseif S[k][3] == Nc-1
                Qtemp = zeros(M);
                
                pi1temp = zeros(M,Ni); 
                pi2temp = zeros(M,Ni); 
                
                for m=1:M
                    MSP_fa_update_RHS(k,t,xval,m);
                    #solve the model
                    optimize!(m_fa[t,k])

                    #check the status 
                    status = termination_status(m_fa[t,k]);
                    
                    if status != MOI.OPTIMAL
                        println(" in Backward Pass")
                        println("Model in stage =", t, " and state = ", k, ", in forward pass is ", status)
                        return 
                    else
                        #collect values
                        Qtemp[m] = objective_value(m_fa[t,k]);
                        for i=1:Ni
                            pi1temp[m,i] = dual(FB1Cons_fa[t,k][i]);
                            pi2temp[m,i] = dual(FB2Cons_fa[t,k][i]);
                        end
                    end                    
                end
                Q[k] = sum(Qtemp)/M;
                for i=1:Ni
                    pi1[k,i] = sum(pi1temp[:,i])/M;
                    pi2[k,i] = sum(pi2temp[:,i])/M;
                end
            else
				# Here we just update xval, not rhs, so doesn't matter which m to use...
                MSP_fa_update_RHS(k,t,xval,rand(1:M));
                #solve the model
                optimize!(m_fa[t,k])

                #check the status 
                status = termination_status(m_fa[t,k]);
                if status != MOI.OPTIMAL
                    println(" in Backward Pass")
                    println("Model in stage =", t, " and state = ", k, ", in forward pass is ", status)
                    return 
                else
                    #collect values
                    Q[k] = objective_value(m_fa[t,k]);
                    for i=1:Ni
                        pi1[k,i] = dual(FB1Cons_fa[t,k][i]);
                        pi2[k,i] = dual(FB2Cons_fa[t,k][i]);
                    end
                end                 
            end
        end

        
        for n = 1:K
            if  n ∉ absorbing_states
                if t-1 == 1 && n != k_init
                    continue
                end
                #what is the expected cost value 
                Qvalue = sum(Q[k]*P_joint[n,k]  for k=1:K);

                # check if cut is violated at the sample path encountered in the forward pass
                if n == sample_n && (Qvalue-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
                    cutviolFlag = 1;
                end

                # we are doing cut sharing so we will add the cut regardless
                @constraint(m_fa[t-1,n],
                ϴ_fa[t-1,n]
                -sum(sum((pi1[k,i]+pi2[k,i])*x_fa[t-1,n][i] for i=1:Ni)*P_joint[n,k] for k=1:K)
                >=
                Qvalue-sum(sum((pi1[k,i]+pi2[k,i])*xval[i,t-1] for i=1:Ni)*P_joint[n,k] for k=1:K)
                );

            end
        end
        
    end
    return cutviolFlag;
end


function train_models_offline()
    #set the RHS of the first_stage problem
    for i=1:Ni
        set_normalized_rhs(FB1Cons_fa[1,k_init][i], x_0[i]);
        set_normalized_rhs(FB2Cons_fa[1,k_init][i], x_0[i]);
    end

    #intialize stuff
    train_time = 0;
    relative_gap = 1e10;
    lb=0;
    LB = [];
    xval = zeros(Ni,T);
    thetaval = zeros(T);
    iter = 0;
    cutviol_iter = 0;
    start=time();
    while true
        iter+=1
        #forward pass
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration(lb,xval,thetaval);
        push!(LB,lb)
        println("LB = ", lb)
        #termination check
        flag, Elapsed = termination_check(iter,relative_gap,LB,start,cutviol_iter);
        if flag != 0
            train_time = Elapsed
            break
        end

        #backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample)
        if cutviolFlag == 1
            cutviol_iter = 0;
        else
            cutviol_iter = cutviol_iter + 1;
        end
    end
    return LB, train_time, iter
end

function FOSDDP_eval_offline()
    
    start=time();
    OS_paths = Matrix(CSV.read("OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    objs_fa = zeros(nbOS,T);
    
    procurmnt_cost = zeros(T); transport_cost = zeros(T); inventory_cost = zeros(T); dshortage_cost = zeros(T);
    desalvege_cost = zeros(T);
    
    xval_fa = Array{Any,2}(undef,nbOS,T); fval_fa = Array{Any,2}(undef,nbOS,T);
    yval_fa = Array{Any,2}(undef,nbOS,T); zval_fa = Array{Any,2}(undef,nbOS,T); vval_fa = Array{Any,2}(undef,nbOS,T);

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
                    
                    procurmnt_cost[t] += (sum(fval_fa[s,t][N0,i] for i=1:Ni)*h[t])/nbOS
                    transport_cost[t] += sum(sum(cb[i,ii,t]*fval_fa[s,t][i,ii] for ii=1:Ni) for i=1:N0)/nbOS
                    inventory_cost[t] += sum(ch[i,t]*xval_fa[s,t][i] for i=1:Ni)/nbOS
                    dshortage_cost[t] += (sum(zval_fa[s,t][j] for j=1:Nj)*p)/nbOS;
                    desalvege_cost[t] += (sum(vval_fa[s,t][i] for i=1:Ni)*q)/nbOS;                               
                end
            end
        end        
    end

    #writeCost("procurmnt_cost",procurmnt_cost,"fa")
    #writeCost("transport_cost",transport_cost,"fa")
    #writeCost("inventory_cost",inventory_cost,"fa")
    #writeCost("dshortage_cost",dshortage_cost,"fa")
    #writeCost("desalvege_cost",desalvege_cost,"fa")
    
    fa_bar = mean(sum(objs_fa[:,t] for t=1:T));
    fa_std = std(sum(objs_fa[:,t] for t=1:T));
    fa_low = fa_bar-1.96*fa_std/sqrt(nbOS);
    fa_high = fa_bar+1.96*fa_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", fa_bar, " ± ", [fa_low,fa_high]);
    elapsed = time() - start;
    vals = [xval_fa, fval_fa, yval_fa, zval_fa, vval_fa];
    
    return objs_fa, fa_bar, fa_low, fa_high, elapsed#, vals
end

##############################################################################################################################
##############################################################################################################################

function RH_2SSP_define_models(t_roll,nbstages1,nbstages2,x_init,ξ)
    #define first stage (master problem) model  
    master, x, f, y1, z1, v1, θ = RH_2SSP_first_stage(t_roll,nbstages1,x_init,ξ)
    
    #define second stage (subproblem) optimality model
    subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons = RH_2SSP_second_stage(t_roll,nbstages1,nbstages2)
    
    #define second stage (subproblem) feasability model
    feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons = feasability_problem(t_roll,nbstages1,nbstages2)
    
    
    return master, x, f, y1, z1, v1, θ, subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons, feas_m, feas_y, feas_z, feas_v, feas_ConsFB, feas_dCons, feas_rCons
end


function RH_2SSP_first_stage(t_roll,nbstages1,x_init,ξ)
    
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
                begin
                   0 <= x[i=1:Ni,t=1:nbstages1] <= x_cap[i]
                   0 <= f[i=1:N0,ii=1:Ni,t=1:nbstages1] <= f_cap[i,ii]
                   0 <= y[i=1:Ni,j=1:Nj]
                   0 <= z[j=1:Nj]
                   0 <= v[i=1:Ni]
                   0 <= θ 
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t_roll-1+t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:nbstages1)
              +sum(sum(ch[i,t_roll-1+t]*x[i,t] for i=1:Ni) for t=1:nbstages1)
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t_roll-1+t] for t=1:nbstages1)
              +sum(sum(ca[i,j,t_roll]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(v[i] for i=1:Ni)*q
              +θ 
               );
    
    #######################
    #Define the constraints.
    for t=1:nbstages1
        for i=1:Ni
            if t == 1
                @constraint(m, x[i,t]
                               +sum(f[i,j,t] for j=1:Ni if j != i)
                               -sum(f[j,i,t] for j=1:N0 if j != i)
                               +sum(y[i,j] for j=1:Nj) 
                               +v[i]
                               == x_init[i]
                            );
                @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x_init[i]);                
            else
               @constraint(m, sum(f[i,j,t] for j=1:Ni if j != i) <= x[i,t-1]);                
               #@constraint(m, x[i,t-1]
               #              -x[i,t]
               #              -sum(f[i,j,t] for j=1:Ni if j != i)
               #              +sum(f[j,i,t] for j=1:N0 if j != i)>=0);
            end
        end
    end
    
    #in case the first stage corresponds to the period where the hurricane make landfall
    for j=1:Nj
        @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) >= ξ[j]);
    end
    
    return m, x, f, y, z, v, θ 
    
end

function RH_2SSP_second_stage(t_roll,nbstages1,nbstages2)
    
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
                begin
                   0 <= y[i=1:Ni,j=1:Nj,t=1:nbstages2]
                   0 <= z[j=1:Nj,t=1:nbstages2]
                   0 <= v[i=1:Ni,t=1:nbstages2]
                   0 >= reimbursement[t=1:nbstages2]
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
              +sum(sum(sum(ca[i,j,t_roll+t]*y[i,j,t] for j=1:Nj) for i=1:Ni) for t=1:nbstages2)
              +sum(sum(z[j,t] for j=1:Nj)*p for t=1:nbstages2)
              +sum(sum(v[i,t] for i=1:Ni)*q for t=1:nbstages2)
              +sum(reimbursement[t] for t=1:nbstages2)
               );

    #######################
    #Define the constraints.
    ConsFB = Dict(); #a dictonary to store first set of constraints
    dCons = Dict(); #a dictonary to store all the demand constraint
    rCons = Dict(); #a dictonary to store all the reimbursement constraint
    
    for t=1:nbstages2
        for i=1:Ni
            ConsFB[t,i] = @constraint(m, sum(y[i,j,t] for j=1:Nj)+v[i,t] == 0); #x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)            
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        rCons[t] = @constraint(m, reimbursement[t] == 0)
    end
    
    return m, y, z, v, reimbursement, ConsFB, dCons, rCons
end


function initialize(nbstages1,nbstages2,s,t_roll)
    LB = -1e10;
    UB = 1e10;
    iter = 0;
    xval = zeros(Ni,nbstages1);
    fval = zeros(N0,Ni,nbstages1);
    θval = 0;
    start=time();
    Elapsed = time()-start;
    nbscen = 1000
    allscen = Matrix(CSV.read("OOS.csv",DataFrame));
    allscenM = Matrix(CSV.read("inOOS.csv",DataFrame))
    scen = allscen[collect(1:convert(Int,10000/nbscen):10000),1:T]
    scenM = allscenM[collect(1:convert(Int,10000/nbscen):10000),1]

    if t_roll > 1
        for n=1:nbscen
            scenM[n] = rand(1:M)
            
            for t=2:t_roll
                scen[n,t] = allscen[s,t]
            end
            
            for t=t_roll+1:T
                scen[n,t] = MC_sample(scen[s,t-1]);
            end
        end
    end
    qprob = fill(1/nbscen,nbscen)
    return LB, UB, iter, xval, fval, θval, start, Elapsed, nbscen, scen, scenM, qprob
end




function solve_first_stage(LB,xval,fval,θval,master,subproblem,x,f,θ)
    optimize!(master) #solve the model
    status_master = termination_status(master); #check the status 
    if status_master != MOI.OPTIMAL
        println("Master problem status is: ", status_master, " === Oops! :/")
        return
        #exit(0);
    else
        #update the values
        LB = objective_value(master);
        xval = value.(x);
        fval = value.(f);
        θval = value(θ);
    end

    return LB, xval, fval, θval
end

function solve_scen_subproblem(Q,pi1,pi2,n,subproblem,ConsFB,dCons,rCons,nbstages2)
    flag = 0
    optimize!(subproblem) #solve the model
    status_subproblem = termination_status(subproblem); #check the status 
    if status_subproblem != MOI.OPTIMAL
        #println("Subproblem problem status is: ", status_subproblem, " === Oops! :/")
        #println("subproblem = ", subproblem)
        if status_subproblem == MOI.INFEASIBLE
            flag = -1
        end
        #exit(0);
    else
        #update the values
        Q[n] = objective_value(subproblem);
        pi1temp = zeros(nbstages2,Ni)
        pi2temp = zeros(nbstages2)
        
        for t=1:nbstages2
            for i=1:Ni
                pi1temp[t,i] = dual(ConsFB[t,i]);
            end
            pi2temp[t] = dual(rCons[t]);            
        end
        pi1[n] = pi1temp
        pi2[n] = pi2temp
    end
    return Q, pi1, pi2, flag
end

function solve_second_stage(t_roll,nbstages1,nbstages2,xval,fval,θval,scen,scenM,nbscen,qprob,master,subproblem,x,f,θ,y,z,v,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons)
    flag = 0;
    Q = zeros(nbscen); #list for all the optimal values
    pi1 = Array{Any,1}(undef,nbscen); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,1}(undef,nbscen); #list for all the dual multiplies of the third set of constraints
    Qbar = 0;
    
    τ = nothing
    for n=1:nbscen
        #identify the period where the hurricane makes landfall 
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);
        
        #update the RHS
        RH_2SSP_update_RHS(τ,scen[n,:],scenM[n],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        
        #solve the subproblem and store the dual information
        Q, pi1, pi2, flag = solve_scen_subproblem(Q,pi1,pi2,n,subproblem,ConsFB,dCons,rCons,nbstages2)
        if flag == -1
            #generate feasability cut
generate_feasability_cut(feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons,n,scen,scenM,nbstages1,xval,fval,t_roll,nbstages2,master,x,f)
            break
        end
    end
    Qbar = 0;
    if flag != -1
        Qbar = sum(Q[n]*qprob[n] for n=1:nbscen)
        pi1bar = sum(pi1[n]*qprob[n] for n=1:nbscen)
        pi2bar = sum(pi2[n]*qprob[n] for n=1:nbscen)

        #add a cut (if any)
        if (Qbar-θval)/max(1e-10,abs(Qbar)) > ϵ
            if τ === nothing
                @constraint(master,
                θ-sum(pi1bar[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
                >= 
             Qbar-sum(pi1bar[t-1,i]*(xval[i,t-1]-xval[i,t]-sum(fval[i,j,t] for j=1:Ni if j != i)+sum(fval[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
                );
            else
                @constraint(master,
                            θ-sum(pi1bar[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
                             -sum(pi2bar[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t>τ)
                            >= 
                         Qbar-sum(pi1bar[t-1,i]*(xval[i,t-1]-xval[i,t]-sum(fval[i,j,t] for j=1:Ni if j != i)+sum(fval[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
                             -sum(pi2bar[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t>τ)
                            );
            end

            flag = 1;
        end
    end
    return flag, Qbar
end


function RH_2SSP_update_RHS(τ,k_t,m,nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
    for t=2:nbstages1
        for i=1:Ni
            set_normalized_rhs(ConsFB[t-1,i],max(0,xval[i,t-1]
                         -xval[i,t]
                         -sum(fval[i,j,t] for j=1:Ni if j != i)
                         +sum(fval[j,i,t] for j=1:N0 if j != i))
                         );
        end
        
        for j=1:Nj
            if t == τ && k_t ∉ absorbing_states 
                set_normalized_rhs(dCons[t-1,j], SCEN[k_t[τ]][j,m]);
            else
                set_normalized_rhs(dCons[t-1,j], 0);
            end
        end
        
        if τ === nothing|| t <= τ
            set_normalized_rhs(rCons[t-1], 0)
        else
            set_normalized_rhs(rCons[t-1],
                -(sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)
                +sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)  
                +sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1])
                                );
        end
    end
end

function RH_2SSP_solve_roll(s,t_roll,nbstages1,nbstages2,master,subproblem,x,f,θ,y,z,v,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons)
    
    LB, UB, iter, xval, fval, θval, start, Elapsed, nbscen, scen, scenM, qprob = initialize(nbstages1,nbstages2,s,t_roll)
    
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ && Elapsed < time_limit 
        iter+=1;
        #println("iter = ", iter)        
        # solve first stage
        LB, xval, fval, θval = solve_first_stage(LB,xval,fval,θval,master,subproblem,x,f,θ)
        #println("LB = ", LB)
        #println("UB = ", UB)
        #println("===================")

        if nbstages2 > 0
            # solve second stage 
            flag, Qbar = solve_second_stage(t_roll,nbstages1,nbstages2,xval,fval,θval,scen,scenM,nbscen,qprob,master,subproblem,x,f,θ,y,z,v,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons)
            if flag != -1
                UB = min(LB-θval+Qbar,UB);
            end
            Elapsed = time()-start;
        else
            break
        end

    end


    return LB, UB, iter, xval, fval, θval, Elapsed
end

function feasability_problem(t_roll,nbstages1,nbstages2)
    
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
                begin
                   0 <= y[i=1:Ni,j=1:Nj,t=1:nbstages2]
                   0 <= z[j=1:Nj,t=1:nbstages2]
                   0 <= v[i=1:Ni,t=1:nbstages2]
                   0 >= reimbursement[t=1:nbstages2]
                   0 <= a_plus[i=1:Ni,t=1:nbstages2]
                   0 <= a_minus[i=1:Ni,t=1:nbstages2]
                end
              );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(a_plus[i,t] for i=1:Ni) for t=1:nbstages2)
              +sum(sum(a_minus[i,t] for i=1:Ni) for t=1:nbstages2)
               );

    #######################
    #Define the constraints.
    ConsFB = Dict(); #a dictonary to store first set of constraints
    dCons = Dict(); #a dictonary to store all the demand constraint
    rCons = Dict(); #a dictonary to store all the reimbursement constraint
    
    for t=1:nbstages2
        for i=1:Ni
            ConsFB[t,i] = @constraint(m, sum(y[i,j,t] for j=1:Nj)+v[i,t]+a_plus[i,t]-a_minus[i,t]== 0); #x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)            
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        rCons[t] = @constraint(m, reimbursement[t] == 0)
    end

    return m, y, z, v, ConsFB, dCons, rCons
end


function generate_feasability_cut(feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons,n,scen,scenM,nbstages1,xval,fval,t_roll,nbstages2,master,x,f)
    #identify the period where the hurricane makes landfall 
    τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);

    #update the RHS
    RH_2SSP_update_RHS(τ,scen[n,:],scenM[n],nbstages1,feas_ConsFB,feas_dCons,feas_rCons,xval,fval,t_roll)
    
    
    optimize!(feas_m) #solve the model
    status_feasability = termination_status(feas_m); #check the status 
    if status_feasability != MOI.OPTIMAL
        println("Feasability problem status is: ", status_feasability, " === Oops! :/")
        println("subproblem = ", feas_m)
        exit(0);
    else
        #update the values
        Q = objective_value(feas_m);
        pi1 = zeros(nbstages2,Ni)
        pi2 = zeros(nbstages2)
        
        for t=1:nbstages2
            for i=1:Ni
                pi1[t,i] = dual(feas_ConsFB[t,i]);
            end
            pi2[t] = dual(feas_rCons[t]);            
        end
        if τ === nothing
            @constraint(master,
            0>=sum(pi1[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
            );
        else
            @constraint(master,        
            0>=sum(pi1[t-1,i]*(x[i,t-1]-x[i,t]-sum(f[i,j,t] for j=1:Ni if j != i)+sum(f[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)+sum(pi2[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t>τ));
        end
    end
end

function static_2SSP_evaluation()
    start=time();
    OS_paths = Matrix(CSV.read("OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    optimize!(master) #solve the model
    #update the values
    objs_st2SSP = fill(objective_value(master)-value(θ),nbOS,T);
    xval = value.(x);
    fval = value.(f);


    for s=1:nbOS
        #identify the period where the hurricane makes landfall 
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
        
        #update the RHS
        RH_2SSP_update_RHS(τ,OS_paths[s,1:T],OS_M[s],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        
        
        Q = zeros(nbOS); #list for all the optimal values
        pi1 = Array{Any,1}(undef,nbOS); #list for all the dual multiplies of the first set of constraints
        pi2 = Array{Any,1}(undef,nbOS); #list for all the dual multiplies of the third set of constraints

        #solve the subproblem and store the dual information
        Q, pi1, pi2, flag = solve_scen_subproblem(Q,pi1,pi2,s,subproblem,ConsFB,dCons,rCons,nbstages2)
        yval = value.(y2);
        zval = value.(z2);
        vval = value.(v2);
        rval = value.(r2);
        for t=1:nbstages2
            objs_st2SSP[s,t+t_roll] = (sum(sum(ca[i,j,t_roll+t]*yval[i,j,t] for j=1:Nj) for i=1:Ni)
                               +sum(zval[j,t] for j=1:Nj)*p
                               +sum(vval[i,t] for i=1:Ni)*q
                               +rval[t])
        end
    end
        
    st2SSP_bar = mean(sum(objs_st2SSP[:,t] for t=1:T));
    st2SSP_std = std(sum(objs_st2SSP[:,t] for t=1:T));
    st2SSP_low = st2SSP_bar-1.96*st2SSP_std/sqrt(nbOS);
    st2SSP_high = st2SSP_bar+1.96*st2SSP_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", st2SSP_bar, " ± ", [st2SSP_low,st2SSP_high]);
    elapsed = time() - start;
    #vals = [xval_st2SSP, fval_st2SSP, yval_st2SSP, zval_st2SSP, vval_st2SSP];
    
    return objs_st2SSP, st2SSP_bar, st2SSP_low, st2SSP_high, elapsed#, vals
    
end

function rolling_horizon()
    start=time();
    s = 1
    t_roll = 1
    nbstages1 = T
    nbstages2 = nbstages1-1 
    x_init = x_0
    ξ = zeros(Nj)
    #define the the model.
    master, x, f, y1, z1, v1, θ, subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons, feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons = RH_2SSP_define_models(t_roll,nbstages1,nbstages2,x_init,ξ)

    #solve the model.
    RH_2SSP_solve_roll(s,t_roll,nbstages1,nbstages2,master,subproblem,x,f,θ,y1,z1,v1,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons);


    OS_paths = Matrix(CSV.read("OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    optimize!(master) #solve the model
    #update the values
    xval = value.(x);
    fval = value.(f);
    x_init = xval[:,1];
    y1val = value.(y1);
    z1val = value.(z1);
    v1val = value.(v1);

    temp_obj = (sum(sum(cb[i,ii,t_roll]*fval[i,ii,1] for ii=1:Ni) for i=1:N0) 
                                     +sum(ch[i,t_roll]*xval[i,1] for i=1:Ni)
                                     +sum(fval[N0,i,1] for i=1:Ni)*h[t_roll]
                                     +sum(sum(ca[i,j,t_roll]*y1val[i,j] for j=1:Nj) for i=1:Ni)
                                     +sum(z1val[j] for j=1:Nj)*p
                                     +sum(v1val[i] for i=1:Ni)*q
                                    );
    objs_RH2SSP = fill(temp_obj,nbOS,T);


    for s=1:nbOS
        println("s = ", s)
        x_init = deepcopy(xval[:,1]);

        for t_roll=2:T-1
            nbstages1 = T-t_roll+1
            nbstages2 = T-t_roll

            ξ = zeros(Nj)
            k_t = OS_paths[s,t_roll]
            if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
                ξ = SCEN[k_t][:,OS_M[s]]
            end
            #define the the model.
            master, x, f, y1, z1, v1, θ, subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons, feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons = RH_2SSP_define_models(t_roll,nbstages1,nbstages2,x_init,ξ)

            #solve the model.
            RH_2SSP_solve_roll(s,t_roll,nbstages1,nbstages2,master,subproblem,x,f,θ,y1,z1,v1,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons);

            #implement xₜ, pay cost, and pass xₜ to new t+1.
            optimize!(master) #solve the model
            
            xval = value.(x);
            fval = value.(f);
            x_init = xval[:,1];
            y1val = value.(y1);
            z1val = value.(z1);
            v1val = value.(v1);

#            objs_RH2SSP[s,t_roll] = objective_value(master)-value(θ);
            objs_RH2SSP[s,t_roll] = (sum(sum(cb[i,ii,t_roll]*fval[i,ii,1] for ii=1:Ni) for i=1:N0) 
                                     +sum(ch[i,t_roll]*xval[i,1] for i=1:Ni)
                                     +sum(fval[N0,i,1] for i=1:Ni)*h[t_roll]
                                     +sum(sum(ca[i,j,t_roll]*y1val[i,j] for j=1:Nj) for i=1:Ni)
                                     +sum(z1val[j] for j=1:Nj)*p
                                     +sum(v1val[i] for i=1:Ni)*q
                                    );
        end

        t_roll = T
        nbstages1 = 1
        nbstages2 = 0
        ξ = zeros(Nj)
        k_t = OS_paths[s,T]
        if S[k_t][3] == Nc-1 && k_t ∉ absorbing_states
            ξ = SCEN[k_t][:,OS_M[s]]
        end
        master, x, f, y1, z1, v1, θ = RH_2SSP_first_stage(t_roll,nbstages1,x_init,ξ)
        #implement xₜ, pay cost, and pass xₜ to new t+1.
        optimize!(master) #solve the model

        objs_RH2SSP[s,t_roll] = objective_value(master)-value(θ);
    end
    
    RH2SSP_bar = mean(sum(objs_RH2SSP[:,t] for t=1:T));
    RH2SSP_std = std(sum(objs_RH2SSP[:,t] for t=1:T));
    RH2SSP_low = RH2SSP_bar-1.96*RH2SSP_std/sqrt(nbOS);
    RH2SSP_high = RH2SSP_bar+1.96*RH2SSP_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", RH2SSP_bar, " ± ", [RH2SSP_low,RH2SSP_high]);
    elapsed = time() - start;
    
    return objs_RH2SSP, RH2SSP_bar, RH2SSP_low, RH2SSP_high, elapsed
    
end