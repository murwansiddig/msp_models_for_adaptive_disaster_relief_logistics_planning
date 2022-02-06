#MSP fully adaptive model functions 

#Define stage-t problem
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

###############################################################
###############################################################

#Define the model
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

###############################################################
###############################################################

#Train model: forward pass
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

###############################################################
###############################################################

#Train model: backward pass
function FOSDDP_backward_pass_oneSP_iteration(lb,xval,thetaval,in_sample)
    cutviolFlag = 0;
    for t=T:-1:2
        #initialize
        Q = zeros(K); #list for all the optimal values
         #list for all the dual multiplies of the first and second set of constraints
        pi1 = zeros(K,Ni); pi2 = zeros(K,Ni); 
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

###############################################################
###############################################################

#Train model
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

###############################################################
###############################################################

#evaluate model
function FOSDDP_eval_offline()
    
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
    fa_bar = mean(sum(objs_fa[:,t] for t=1:T));
    fa_std = std(sum(objs_fa[:,t] for t=1:T));
    fa_low = fa_bar-1.96*fa_std/sqrt(nbOS);
    fa_high = fa_bar+1.96*fa_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", fa_bar, " ± ", [fa_low,fa_high]);
    elapsed = time() - start;
    vals = [xval_fa, fval_fa, yval_fa, zval_fa, vval_fa];
    
    return objs_fa, fa_bar, fa_low, fa_high, elapsed#, vals
end


###############################################################
###############################################################

#update RHS of flow-balance and demand constraint
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

