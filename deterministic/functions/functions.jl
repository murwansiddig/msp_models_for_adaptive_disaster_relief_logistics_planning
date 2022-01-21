# This function defines a non-terminal stage for the MSP model with deterministic landfall
function non_terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,t)
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni] <= f_cap[i,ii] 
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
              +ϴ
               );

    #######################
    #Define the constraints.
    
    #initialize two arrays to store the constraints which containt x_{t-1}
    FB1 = Array{Any,1}(undef,Ni);
    FB2 = Array{Any,1}(undef,Ni);
    
    #create the following constraints for every SP i
    #initialize the RHS to be zero for now, and will change it later 
    for i=1:Ni
        FB1[i] = @constraint(m, 
                             x[i]
                            +sum(f[i,j] for j=1:Ni if j != i)
                            -sum(f[j,i] for j=1:N0 if j != i)
                            == 0
                            );
        FB2[i] = @constraint(m, 
                            sum(f[i,j] for j=1:Ni if j != i)
                            <= 0
                            );
    end

    return m, x, f, ϴ, FB1, FB2
end


# This function defines the terminal stage for the MSP model with deterministic landfall
function terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,t)
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
               0 <= d[j=1:Nj]
               0 <= z[j=1:Nj]
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
              +sum(x[i]-sum(y[i,j] for j=1:Nj) for i=1:Ni)*q
               );
    

    #######################
    #Define the constraints.
    
    #initialize two arrays to store the constraints which containt x_{t-1}
    FB1 = Array{Any,1}(undef,Ni);
    FB2 = Array{Any,1}(undef,Ni);
    
    #create the following constraints for every SP i
    #initialize the RHS to be zero for now, and will change it later 
    for i=1:Ni
        FB1[i] = @constraint(m, 
                             x[i]
                            +sum(f[i,j] for j=1:Ni if j != i)
                            -sum(f[j,i] for j=1:N0 if j != i)
                            == 0
                            );
        FB2[i] = @constraint(m, 
                            sum(f[i,j] for j=1:Ni if j != i)
                            <= 0
                            );
        @constraint(m,
                    sum(y[i,j] for j=1:Nj)
                    <=x[i]
                   );
    end
    for j=1:Nj
        # this constraint ensures that the flow sent to an SP is not more than the realized demand
        @constraint(m,
                    sum(y[i,j] for i=1:Ni)
                    <=d[j]
                   );
        # this constraint computes the unsatisfied demand
        @constraint(m,
                    d[j]-sum(y[i,j] for i=1:Ni)
                    <=z[j]
                   );
    end
    
    # this constraint ensures that items cannot be shipped directly from the MDC to SP
    @constraint(m,
            sum(f[N0,i] for i=1:Ni) == 0
            );
    
    return m, x, f, y, d, z, FB1, FB2
end

# This function defines all the models for the MSP with deterministic landfall
function define_models(K,T,Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q)
    # first we initialize the list where we store all the models
    model = Array{Any,2}(undef,T,K); # list to store all the models for every stage and Markovian state
    x = Array{Any,2}(undef,T,K); # list to store all the x variables for every stage and Markovian state
    f = Array{Any,2}(undef,T,K); # ..............
    theta = Array{Any,2}(undef,T-1,K); # ..............
    y = Array{Any,1}(undef,K);
    d = Array{Any,1}(undef,K);
    z = Array{Any,1}(undef,K);
    FB1 = Array{Any,2}(undef,T,K);
    FB2 = Array{Any,2}(undef,T,K);
    
    #Then we define the a model for every stage and Markovian state state 
    for k=1:K, t=1:T
        if t < T
            # define use the function non_terminal_stage_single_period_problem
            model[t,k], x[t,k], f[t,k], theta[t,k], FB1[t,k], FB2[t,k] = non_terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,t);
        else
            #use the function terminal_stage_single_period_problem
            model[t,k], x[t,k], f[t,k], y[k], d[k], z[k], FB1[t,k], FB2[t,k] = terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,t);
        end
    end
    return model, x, f, theta, y, d, z, FB1, FB2
end



# this is a function which takes the current Markovian state as input and returns a randomly generaterd MC
function MC_sample(current_state)
    states = collect(1:1:K);
    k = sample(states, Weights(P_joint[current_state,:]));
    return k
end


# this is a function which updates the RHS of the flow balance constraints in the MSP model with deterministic landfall
function update_RHS(k_t,t,FB1,FB2,xval)
    for i=1:Ni
        set_normalized_rhs(FB1[t,k_t][i], xval[i,t-1]);
        set_normalized_rhs(FB2[t,k_t][i], xval[i,t-1]);
    end
end


# this function runs one iteration of the forward pass of the MC-SDDP in the MSP model with deterministic landfall
function FOSDDP_forward_pass_oneSP_iteration(model,x,f,theta,y,d,z,FB1,FB2,k,T,xval,thetaval,lb)
    
    #initial list where we will store the sample path encountered during the forward pass
    in_sample = [];
    #what is the state in the first stage
    k_t = k; # k is the current stage and is taken as an input
    for t=1:T
        #the state is known in the first stage; if not sample a new state k 
        if t>1
            #sample a new state k
            k_temp = in_sample[t-1];#the state at the previous period
            k_t = MC_sample(k_temp); #sample the new state k_t
            
            #update the RHS
            update_RHS(k_t,t,FB1,FB2,xval);
        end
        
        #solve the model
        optimize!(model[t,k_t])

        #check the status 
        status = termination_status(model[t,k_t]);
        if status != MOI.OPTIMAL
            println(" in Forward Pass")
            println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
            exit(0);
        else
            #collect values
            xval[:,t] = value.(x[t,k_t]);
            if t == 1
                #update the lower bound value
                lb = objective_value(model[t,k_t]);
            end
            if t < T
                # if it's a non-terminal stage update the value of theta
                thetaval[t] = value(theta[t,k_t])
            end
        end
        # add the current state to the sample path
        push!(in_sample,k_t)
    end
    
    return xval, thetaval, lb, in_sample
end

function FOSDDP_backward_pass_oneSP_iteration(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample)
    cutviolFlag = 0;
    for t=T:-1:2
        if t == T
            flag1 = terminal_stage_cut(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample,t);
            if flag1 == 1
                cutviolFlag = 1;
            end
        else
            flag2 = non_terminal_stage_cut(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample,t);
            if flag2 == 1
                cutviolFlag = 1;
            end
        end        
    end
    return cutviolFlag;
end



# this a function which runs check for termination based on 4 criterion: time, number of iterations, LB progress and cut violation 
function termination_check(iter,relative_gap,LB,start,cutviol_iter)
    flag = 0;
    if iter > stall
        relative_gap = (LB[iter]-LB[iter-stall])/max(1e-10,abs(LB[iter-stall]));
    end
    Elapsed = time() - start;
    #println("cutviol_iter = ", cutviol_iter);
    if relative_gap < ϵ ||  iter > max_iter || Elapsed > time_limit || cutviol_iter > cutviol_maxiter
        if relative_gap > ϵ
            if iter > max_iter
                println("terminated because iteration limit was hit")
            else
                if Elapsed > time_limit
                    println("terminated because time limit was hit")
                else
                    println("terminated because cut violation stopping criterion was hit");
                end
            end
        end
        flag = 1;
    end
    return flag, Elapsed, relative_gap
end


#This is a function which generates a cut from T to T-1 in the Backward pass
function terminal_stage_cut(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample,t)
    #initialize
    returnval = 0;
    Q = Array{Any,2}(undef,K,M); #list for all the optimal values
    Q_prob = Array{Any,2}(undef,K,M); #list for all the probabilities
    pi1 = Array{Any,2}(undef,K,M); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,2}(undef,K,M); #list for all the dual multiplies of the second set of constraints

    sample_n = in_sample[t-1]; #the states observed at time t-1
    for k=1:K
        update_RHS(k,t,FB1,FB2,xval)
        for m=1:M
            #calculate the conditional probability
            for j=1:Nj
                fix(d[k][j], SCEN[k][j,m]; force=true);
            end

            #solve the model
            optimize!(model[t,k])

            #check the status 
            status = termination_status(model[t,k]);
            if status != MOI.OPTIMAL
                println(" in Backward Pass")
                println("Model in stage =", t, " and state = ", k, ", in forward pass is ", status)
                exit(0);
            else
                #collect values
                Q[k,m] = objective_value(model[t,k]);
                pi1[k,m] = dual.(FB1[t,k]);
                pi2[k,m] = dual.(FB2[t,k]);
            end                    
        end
    end
    for n = 1:K
        #what is the expected cost value 
        Qvalue = sum(sum(Q[k,m]*P_joint[n,k]*1.0/M for m=1:M) for k=1:K);    
        
        # check if cut is violated at the sample path encountered in the forward pass
        if n == sample_n && (Qvalue-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
            returnval = 1;
        end
        
        # we are doing cut sharing so we will add the cut regardless
        @constraint(model[t-1,n],
            theta[t-1,n]
            -sum(sum(sum((pi1[k,m][i]+pi2[k,m][i])*x[t-1,n][i] for i=1:Ni)*P_joint[n,k]*1.0/M for m=1:M) for k=1:K)
            >=
            Qvalue-sum(sum(sum((pi1[k,m][i]+pi2[k,m][i])*xval[i,t-1] for i=1:Ni)*P_joint[n,k]*1.0/M for m=1:M) for k=1:K)
        );
        
    end 
    return returnval;
end

#This is a function which generates a cut from any t<T to t-1 in the Backward pass
function non_terminal_stage_cut(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample,t)
    #initialize
    returnval = 0;
    Q = Array{Any,1}(undef,K); #list for all the optimal values
    Q_prob = Array{Any,1}(undef,K); #list for all the probabilities
    pi1 = Array{Any,1}(undef,K); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,1}(undef,K); #list for all the dual multiplies of the second set of constraints
    
    sample_n = in_sample[t-1]; #the states observed at time t-1
    for k=1:K
        update_RHS(k,t,FB1,FB2,xval)
        #solve the model
        optimize!(model[t,k])

        #check the status 
        status = termination_status(model[t,k]);
        if status != MOI.OPTIMAL
            println(" in Backward Pass")
            println("Model in stage =", t, " and state = ", k, ", in forward pass is ", status)
            exit(0);
        else
            #collect values
            Q[k] = objective_value(model[t,k]);
            pi1[k] = dual.(FB1[t,k]);
            pi2[k] = dual.(FB2[t,k]);
        end                    
    end
    
    for n = 1:K
        #what is the expected cost value 
        Qvalue = sum(Q[k]*P_joint[n,k]  for k=1:K);

        
        # check if cut is violated at the sample path encountered in the forward pass
        if n == sample_n && (Qvalue-thetaval[t-1])/max(1e-10,abs(thetaval[t-1])) > ϵ
            returnval = 1;
        end
        
        # we are doing cut sharing so we will add the cut regardless
        @constraint(model[t-1,n],
        theta[t-1,n]
        -sum(sum((pi1[k][i]+pi2[k][i])*x[t-1,n][i] for i=1:Ni)*P_joint[n,k] for k=1:K)
        >=
        Qvalue-sum(sum((pi1[k][i]+pi2[k][i])*xval[i,t-1] for i=1:Ni)*P_joint[n,k] for k=1:K)
        );
        
    end 
    return returnval; 
end



function train_models_offline(model,x,f,theta,y,d,z,FB1,FB2,k,T)
    #set the RHS of the first_stage problem
    for i=1:Ni
        set_normalized_rhs(FB1[1,k][i], x_0[i]);
        set_normalized_rhs(FB2[1,k][i], x_0[i]);
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
        xval, thetaval, lb, in_sample = FOSDDP_forward_pass_oneSP_iteration(model,x,f,theta,y,d,z,FB1,FB2,k,T,xval,thetaval,lb);
        push!(LB,lb)
        println("LB = ", lb)
        #termination check
        flag, Elapsed, relative_gap = termination_check(iter,relative_gap,LB,start,cutviol_iter);
        if flag == 1
            train_time = Elapsed
            break
        end
        #backward pass (if not terminated)
        cutviolFlag = FOSDDP_backward_pass_oneSP_iteration(model,x,f,theta,y,d,z,FB1,FB2,T,xval,thetaval,in_sample);
        if cutviolFlag == 1
            cutviol_iter = 0;
        else
            cutviol_iter = cutviol_iter + 1;
        end
    end
    return LB, train_time, iter
end



# this function runs SDDP for the MSP model with deterministic landfall
function FOSDDP_eval_offline(model,x,f,theta,y,d,z,FB1,FB2,k,T,nbOS)
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    actions = [];
    start=time();
    costs = zeros(nbOS,T);
    Z = [];
        
    procurmnt_amount = zeros(T);
    transport_amount = zeros(T);
    inventory_amount = zeros(T);
    dshortage_amount = 0;
    desalvege_amount = 0;
    
    for s=1:nbOS
        temp_action = [];
        cx = 0;
        xval = zeros(Ni,T);
    
        for t=1:T
            k_t = OS_paths[s,t];
            #the state is known in the first stage; if not sample a new state k 
            if t>1
                #update the RHS
                update_RHS(k_t,t,FB1,FB2,xval);
            end
            if t == T
                m = OS_M[s];
                for j=1:Nj
                    fix(d[k_t][j], SCEN[k_t][j,m]; force=true);
                end
            end 
            #solve the model
            optimize!(model[t,k_t])

            #check the status 
            status = termination_status(model[t,k_t]);
            if status != MOI.OPTIMAL
                println(" in Forward Pass")
                println("Model in stage =", t, " and state = ", k_t, ", in forward pass is ", status)
                exit(0);
            else
                #collect values
                xval[:,t] = value.(x[t,k_t]);
                
                procurmnt_amount[t] += (sum(value.(f[t,k_t])[N0,i] for i=1:Ni))/nbOS
                transport_amount[t] += sum(sum(value.(f[t,k_t])[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
                inventory_amount[t] += sum(xval[i,t] for i=1:Ni)/nbOS
                
                if t < T
                    costs[s,t] = objective_value(model[t,k_t]) - value(theta[t,k_t]);
                    push!(temp_action,value.(f[t,k_t])[N0,:]);
                else
                    costs[s,t] = objective_value(model[t,k_t]);
                    push!(temp_action,-value.(z[k_t]));
                    dshortage_amount += (sum(value.(z[k_t])[j] for j=1:Nj))/nbOS;
                    desalvege_amount += (sum(xval[i,t]-sum(value.(y[k_t])[i,j] for j=1:Nj) for i=1:Ni))/nbOS;
                end
                cx+=costs[s,t];
            end
        end
        push!(Z,cx)
        push!(actions,temp_action)
    end
    cfa_fname = string("./output/cost_factor_analysis.csv")
    df = CSV.read(cfa_fname,DataFrame);
    reSults = convert.(Float64,Matrix(df)); 
    reSults[inst,1:T] = procurmnt_amount
    updf = DataFrame(reSults, :auto);
    CSV.write(cfa_fname,updf)
    
    
    COSTs = [procurmnt_amount, transport_amount, inventory_amount, dshortage_amount, desalvege_amount];
    elapsed = time() - start;
    
    UB_bar = mean(Z);
    UB_std = std(Z);
  
    UB_low = UB_bar-1.96*UB_std/sqrt(nbOS);
    UB_high = UB_bar+1.96*UB_std/sqrt(nbOS);
    
    println("μ = ", UB_bar);
    println("μ ± 1.96*σ/√NS = ", [UB_low,UB_high]);

    return costs, UB_bar, UB_low, UB_high, elapsed, actions, COSTs
end

function deterministic_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,x_0)
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables. Note that we disable a "direct shipment" from MDC to DPs at the terminal stage, hence y variables are not defined for N0
    @variables(m,
            begin
               0 <= x[i=1:Ni,t=1:T] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:T] <= f_cap[i,ii]
               0 <= y[i=1:Ni,j=1:Nj]
               0 <= d[j=1:Nj]
               0 <= z[j=1:Nj]
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:T)
              +sum(sum(ch[i,t]*x[i,t] for i=1:Ni) for t=1:T)
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t] for t=1:T)
              +sum(sum(ca[i,j,T]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(x[i,T]-sum(y[i,j] for j=1:Nj) for i=1:Ni)*q
               );
    
    #######################
    #Define the constraints.
    for i=1:Ni
        for t=1:T
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x_0[i]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x_0[i]
                            );
            else
                
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x[i,t-1]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x[i,t-1]
                            );                
            end
        end
       @constraint(m,
            sum(y[i,j] for j=1:Nj)
            <=x[i,T]
           );
    end
    for j=1:Nj
        @constraint(m,
                    sum(y[i,j] for i=1:Ni)
                    <=d[j]
                   );
        @constraint(m,
                    d[j]-sum(y[i,j] for i=1:Ni)
                    <=z[j]
                   );
    end
	# Assuming no flow from MDC at the terminal stage
    @constraint(m,
        sum(f[N0,i,T] for i=1:Ni) == 0
        );
    
    return m, x, f, y, d, z
end


function clairvoyant_eval(model_cv,x_cv,f_cv,y_cv,d_cv,z_cv,T,nbOS)
    
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    k = OS_paths[1,1];
    
    
    start=time();
    costs = zeros(nbOS);
    
    for s=1:nbOS   
        k_t = OS_paths[s,T]
        m = OS_M[s];
        for j=1:Nj
            fix(d_cv[j], SCEN[k_t][j,m]; force=true);
        end
        
        #solve the model
        optimize!(model_cv)
        
        #check the status 
        status = termination_status(model_cv);
        if status != MOI.OPTIMAL
            println(" Model in clairvoyant is ", status)
            exit(0);
        else
            #collect values
            costs[s] = objective_value(model_cv)
        end
    end
    
    
    cv_bar = mean(costs);
    cv_std = std(costs);
  
    cv_low = cv_bar-1.96*cv_std/sqrt(nbOS);
    cv_high = cv_bar+1.96*cv_std/sqrt(nbOS);
    
    println("μ = ", cv_bar);
    println("μ ± 1.96*σ/√NS = ", [cv_low,cv_high]);
    
    elapsed = time() - start;
    return costs, cv_bar, cv_low, cv_high, elapsed 
end

function MVP_solve(model_mvp,x_mvp,f_mvp,y_mvp,d_mvp,z_mvp,k,T)
    start=time();
    xval_mvp = zeros(Ni,T);
    
    P_temp = P_joint^(T-1);
    P_temp_c = deepcopy(P_temp);
    #normalize the probabilities
    for k=1:K, kk=1:K
        P_temp[k,kk] = P_temp_c[k,kk]/sum(P_temp_c[k,:])
    end
    d_bar = sum(sum(SCEN[kk][:,m]*(1/M) for m=1:M)*P_temp[k,kk] for kk=1:K);
    
    for j=1:Nj
        fix(d_mvp[j], d_bar[j]; force=true);
    end
    
    #solve the model
    optimize!(model_mvp)
    
    #check the status 
    status = termination_status(model_mvp);
    if status != MOI.OPTIMAL
        println(" Model in mean value problem is ", status)
        exit(0);
    else
        #collect values
        xval_mvp = value.(x_mvp)
    end
    
    elapsed_solve=time()-start;
    
    return xval_mvp, elapsed_solve
    
end


function MVP_eval(model_mvp,x_mvp,f_mvp,y_mvp,d_mvp,z_mvp,k,T,xval_mvp,nbOS,OS_paths,OS_M)
    start=time();
    costs = zeros(nbOS);
    
    #fix all the x values using the MVP solution 
    for i=1:Ni, t=1:T
        fix(x_mvp[i,t], xval_mvp[i,t]; force=true);
    end
    
    for s=1:nbOS   
        k_t = OS_paths[s,T]
        m = OS_M[s];
        for j=1:Nj
            fix(d_mvp[j], SCEN[k_t][j,m]; force=true);
        end
        
        #solve the model
        optimize!(model_mvp)
        
        #check the status 
        status = termination_status(model_mvp);
        if status != MOI.OPTIMAL
            println(" Model in clairvoyant is ", status)
            exit(0);
        else
            #collect values
            costs[s] = objective_value(model_mvp)
        end
    end
    
    mvp_bar = mean(costs);
    mvp_std = std(costs);
  
    mvp_low = mvp_bar-1.96*mvp_std/sqrt(nbOS);
    mvp_high = mvp_bar+1.96*mvp_std/sqrt(nbOS);
    
    println("μ = ", mvp_bar);
    println("μ ± 1.96*σ/√NS = ", [mvp_low,mvp_high]);
    
    elapsed_eval = time() - start;
    
    return costs, mvp_bar, mvp_low, mvp_high, elapsed_eval 
end


function static_twostage_first_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,x_0)
    #######################
    #Define the first-stage (multi-period) model for the static two-stage SP model
     m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x[i=1:Ni,t=1:(T-1)] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:(T-1)] <= f_cap[i,ii]
               0 <= theta
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:(T-1))
              +sum(sum(ch[i,t]*x[i,t] for i=1:Ni) for t=1:(T-1))
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t] for t=1:(T-1))
              +theta
               );
    
    #######################
    #Define the constraints.
    for i=1:Ni
        for t=1:(T-1)
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x_0[i]
                            );
                @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x_0[i]
                            );
            else
                
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x[i,t-1]
                            );
                @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x[i,t-1]
                            );                
            end
        end
    end
    return m, x, f, theta
end

function rolling_twostage_first_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,t0,T)
    #######################
    #Define the first-stage (multi-period) model for each roll of the rolling two-stage SP model, assuming that the current stage is t0
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));
    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x_prev[i=1:Ni]
               0 <= x[i=1:Ni,t=1:(T-t0)] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:(T-t0)] <= f_cap[i,ii]
               0 <= ϴ
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t+t0-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:(T-t0))
              +sum(sum(ch[i,t+t0-1]*x[i,t] for i=1:Ni) for t=1:(T-t0))
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t+t0-1] for t=1:(T-t0))
              +ϴ
               );
    
    #######################
    #Define the constraints.
    for i=1:Ni
        for t=1:(T-t0)
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x_prev[i]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x_prev[i]
                            );
            else
                
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x[i,t-1]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x[i,t-1]
                            );                
            end
        end
    end
    return m, x, ϴ, x_prev, f
end

function rolling_twostage_eval_online(T,nbOS,x0)

    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    k = OS_paths[1,1];    
    
    start=time();
    costs = zeros(nbOS,T);
    Z = [];
    actions = [];
    
    procurmnt_amount = zeros(T);
    transport_amount = zeros(T);
    inventory_amount = zeros(T);
    dshortage_amount = 0;
    desalvege_amount = 0;
    
    # define the terminal-stage problem
    model_T, x_T, f_T, y_T, d_T, z_T, FB1_T, FB2_T = terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T);

    # solve the first roll model (same for all sample paths!)
    m1, x1, theta1, x_prev1, f1 = rolling_twostage_first_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,1,T);
    for i=1:Ni
        fix(x_prev1[i], x0[i]; force=true);
    end

    init_k = OS_paths[1,1];
    flag1 = 1;
    xval1 = zeros(Ni,T-1);
    thetaval1 = 0;
    LB1 = 0;
    imm_amount1 = 0;
    
    while flag1 == 1
        optimize!(m1);
        status = termination_status(m1);
        if status != MOI.OPTIMAL
            println(" in master problem")
            exit(0);
        else
            #collect values
            LB1 = objective_value(m1);
            xval1 = value.(x1);
            fval1 = value.(f1);
            thetaval1 = value(theta1);
            imm_amount1 = (sum(sum(cb[i,ii,1]*fval1[i,ii,1] for ii=1:Ni) for i=1:N0)
             +sum(ch[i,1]*xval1[i,1] for i=1:Ni)
             +sum(fval1[N0,i,1] for i=1:Ni)*h[1]);
        end
        flag1 = terminal_stage_cut_two_stage(m1,x1,theta1,model_T,d_T,FB1_T,FB2_T,1,T,xval1,thetaval1,init_k);
    end
    

    f_val1 = value.(f1)[N0,:,1];
    fval_all1 = value.(f1)[:,:,1];
    for s=1:nbOS
        # starting the rolling horizon procedure one for each sample path
        println("Start sample path #", s); 
        
        cx = imm_amount1;
        costs[s,1] = cx;
        xvals = zeros(Ni,T); # xval records the state variables in each stage on the sample path
        xvals[:,1] = xval1[:,1]; # the first one has been computed (same for all sample paths)
        
        
        temp_action = [];
        push!(temp_action,f_val1)
        
        procurmnt_amount[1] += (sum(f_val1[i] for i=1:Ni)*h[1])/nbOS
        transport_amount[1] += sum(sum(cb[i,ii,1]*fval_all1[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
        inventory_amount[1] += sum(ch[i,1]*xval1[i] for i=1:Ni)/nbOS
        
        for t=2:(T-1)
            k_t = OS_paths[s,t];
            # create a two-stage SP model for each stage t
            m_roll, x_roll, theta_roll, x_prev, f_roll = rolling_twostage_first_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,t,T);
            #update the RHS
            for i=1:Ni
                fix(x_prev[i], xvals[i,t-1]; force=true);
            end
            #solve the two-stage model
            init_k = k_t;
            flag = 1;
            xval_running = zeros(Ni,T-t);
            thetaval = 0;
            LB = 0;
            imm_amount_roll=0;
            while flag == 1
                optimize!(m_roll);
                status = termination_status(m_roll);
                if status != MOI.OPTIMAL
                    println(" in the final out-of-sample evaluation at roll = ", t);
                    exit(0);
                else
                    #collect values
                    LB = objective_value(m_roll);
                    xval_running = value.(x_roll);
                    fval_running = value.(f_roll);
                    thetaval = value(theta_roll);
                    imm_amount_roll = (sum(sum(cb[i,ii,t]*fval_running[i,ii,1] for ii=1:Ni) for i=1:N0)
                                     +sum(ch[i,t]*xval_running[i,1] for i=1:Ni)
                                     +sum(fval_running[N0,i,1] for i=1:Ni)*h[t]);
                end
                flag = terminal_stage_cut_two_stage(m_roll,x_roll,theta_roll,model_T,d_T,FB1_T,FB2_T,t,T,xval_running,thetaval,init_k);
            end

            xvals[:,t] = xval_running[:,1];
            costs[s,t] = imm_amount_roll;
            cx+=costs[s,t];
            f_val_roll = value.(f_roll)[:,:,1];
            
            push!(temp_action,value.(f_roll)[N0,:,1])
    
            procurmnt_amount[t] += (sum(f_val_roll[N0,i] for i=1:Ni)*h[t])/nbOS
            transport_amount[t] += sum(sum(cb[i,ii,t]*f_val_roll[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
            inventory_amount[t] += sum(ch[i,t]*xval_running[i] for i=1:Ni)/nbOS
        end

        # last stage?
        for i=1:Ni
            set_normalized_rhs(FB1_T[i], xvals[i,T-1]);
            set_normalized_rhs(FB2_T[i], xvals[i,T-1]);
        end    
        for j=1:Nj
            fix(d_T[j], SCEN[OS_paths[s,T]][j,OS_M[s]]; force=true);
        end
        #solve the model
        optimize!(model_T);
        status = termination_status(model_T);
        if status != MOI.OPTIMAL
            println(" in the final out-of-sample evaluation at roll = ", T);
            exit(0);
        else
            costs[s,T] = objective_value(model_T);
            cx += costs[s,T];
        end
        push!(Z,cx)
        push!(temp_action,-value.(z_T))
        push!(actions,temp_action)
        fval_T = value.(f_T);
        xval_T = value.(x_T);
        zval_T = value.(z_T);
        yval_T = value.(y_T);
      
        procurmnt_amount[T] += (sum(fval_T[N0,i] for i=1:Ni)*h[T])/nbOS
        transport_amount[T] += sum(sum(cb[i,ii,T]*fval_T[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
        inventory_amount[T] += sum(ch[i,T]*xval_T[i] for i=1:Ni)/nbOS

        dshortage_amount += (sum(zval_T[j] for j=1:Nj)*p)/nbOS;
        desalvege_amount += (sum(xval_T[i]-sum(yval_T[i,j] for j=1:Nj) for i=1:Ni)*q)/nbOS;
    end
    
    COSTs = [procurmnt_amount, transport_amount, inventory_amount, dshortage_amount, desalvege_amount];

    elapsed = time() - start;
    
    UB_bar = mean(Z);
    UB_std = std(Z);
  
    UB_low = UB_bar-1.96*UB_std/sqrt(nbOS);
    UB_high = UB_bar+1.96*UB_std/sqrt(nbOS);
    
    println("μ = ", UB_bar);
    println("μ ± 1.96*σ/√NS = ", [UB_low,UB_high]);
    
    return costs, UB_bar, UB_low, UB_high, elapsed, actions, COSTs
end



function terminal_stage_cut_two_stage(model,x,theta,model_sub,d,FB1,FB2,t,T,xval,thetaval,init_k)
    # Just a single model: model_sub, FB1, FB2, etc. passed along here
    # t: when this two-stage SP model is solved
    returnval = 0;
    Q = Array{Any,2}(undef,K,M); #list for all the optimal values
    pi1 = Array{Any,2}(undef,K,M); #list for all the dual multiplies of the first set of constraints
    pi2 = Array{Any,2}(undef,K,M); #list for all the dual multiplies of the second set of constraints

    for k=1:K
        for m=1:M
            for i=1:Ni
                set_normalized_rhs(FB1[i], xval[i,T-t]);
                set_normalized_rhs(FB2[i], xval[i,T-t]);
            end
            #calculate the conditional probability
            for j=1:Nj
                fix(d[j], SCEN[k][j,m]; force=true);
            end

            #solve the model
            optimize!(model_sub)

            #check the status 
            status = termination_status(model_sub);
            if status != MOI.OPTIMAL
                println(" in static two stage cut generation")
                exit(0);
            else
                #collect values
                Q[k,m] = objective_value(model_sub);
                pi1[k,m] = dual.(FB1);
                pi2[k,m] = dual.(FB2);
            end                    
        end
    end
    Qvalue = sum(sum(Q[k,m]*P_terminals[t][init_k,k]*1.0/M for m=1:M) for k=1:K);    

    # check if cut is violated at the sample path encountered in the forward pass
    if (Qvalue-thetaval)/max(1e-10,abs(thetaval)) > ϵ
        #add a cut if any needed
        @constraint(model,
        theta
        -sum(sum(sum((pi1[k,m][i]+pi2[k,m][i])*x[i,T-t] for i=1:Ni)*P_terminals[t][init_k,k]*1.0/M for m=1:M) for k=1:K)
        >=
        Qvalue-sum(sum(sum((pi1[k,m][i]+pi2[k,m][i])*xval[i,T-t] for i=1:Ni)*P_terminals[t][init_k,k]*1.0/M for m=1:M) for k=1:K)
        );
        returnval = 1;
    end
    return returnval;
end



function two_stage_static_train()
    #define the first-stage problem 
    model_2s, x_2s, f_2s, theta_2s = static_twostage_first_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,x_0);

    # define the terminal-stage problem
    model_T, x_T, f_T, y_T, d_T, z_T, FB1_T, FB2_T = terminal_stage_single_period_problem(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T);
    
    init_k = copy(k_init);
    flag = 1;
    xval = zeros(Ni,(T-1));
    thetaval = 0;
    LB = 0;
    start = time();
    iter = 0;
        
    while flag == 1
        iter = iter + 1;
        optimize!(model_2s);
        status = termination_status(model_2s);
        if status != MOI.OPTIMAL
            println(" in master problem")
            exit(0);
        else
            #collect values
            LB = objective_value(model_2s);
            xval = value.(x_2s);
            thetaval = value(theta_2s);
            #println("LB = ", LB);
            #println("thetaval = ", thetaval);
        end
        flag = terminal_stage_cut_two_stage(model_2s,x_2s,theta_2s,model_T,d_T,FB1_T,FB2_T,1,T,xval,thetaval,init_k);
    end
    
    Elapsed = time()-start;
    
    return LB, thetaval, iter, xval, Elapsed;
    
end


function two_stage_static_eval(xval_twoSP)
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    
    model_twoSP, x_twoSP, f_twoSP, y_twoSP, d_twoSP, z_twoSP = deterministic_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,x_0);
  
    start=time();
    costs = zeros(nbOS);

    actions = [];
    procurmnt_amount = zeros(T);
    transport_amount = zeros(T);
    inventory_amount = zeros(T);
    dshortage_amount = 0;
    desalvege_amount = 0;

    #fix all the x values using the twoSP solution 
    for i=1:Ni, t=1:T-1
        fix(x_twoSP[i,t], xval_twoSP[i,t]; force=true);
    end
    
    for s=1:nbOS  
        k_t = OS_paths[s,T]
        m = OS_M[s];
        for j=1:Nj
            fix(d_twoSP[j], SCEN[k_t][j,m]; force=true);
        end
        
        #solve the model
        optimize!(model_twoSP)
        
        #check the status 
        status = termination_status(model_twoSP);
        if status != MOI.OPTIMAL
            println(" Model in clairvoyant is ", status)
            exit(0);
        else
            #collect values
            costs[s] = objective_value(model_twoSP)
        end
        
        
        temp_action = [];
        xval_all = value.(x_twoSP)
       
        for t=1:T
            fval_all = value.(f_twoSP)[:,:,t]
            procurmnt_amount[t] += (sum(fval_all[N0,i] for i=1:Ni)*h[t])/nbOS
            transport_amount[t] += sum(sum(cb[i,ii,t]*fval_all[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
            inventory_amount[t] += sum(ch[i,t]*xval_all[i,t] for i=1:Ni)/nbOS
            
            if t < T
                f_valt = value.(f_twoSP)[N0,:,t]
                push!(temp_action,f_valt)
            else
                z_valt = value.(z_twoSP);
                push!(temp_action,-z_valt)
                
                dshortage_amount += (sum(z_valt[j] for j=1:Nj)*p)/nbOS;
                desalvege_amount += (sum(xval_all[i,t]-sum(value.(y_twoSP)[i,j] for j=1:Nj) for i=1:Ni)*q)/nbOS;
                
            end
        end

        
        push!(actions,temp_action)

    end
    COSTs = [procurmnt_amount, transport_amount, inventory_amount, dshortage_amount, desalvege_amount];

    twoSP_bar = mean(costs);
    twoSP_std = std(costs);
  
    twoSP_low = twoSP_bar-1.96*twoSP_std/sqrt(nbOS);
    twoSP_high = twoSP_bar+1.96*twoSP_std/sqrt(nbOS);
    
    println("μ = ", twoSP_bar);
    println("μ ± 1.96*σ/√NS = ", [twoSP_low,twoSP_high]);
    
    elapsed_eval = time() - start;
    
    return costs, twoSP_bar, twoSP_low, twoSP_high, elapsed_eval, actions, COSTs
end



function deterministic_RH_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,t0)
    #######################
    #Define the model, assuming that the current stage is t0, Note that we disable a "direct shipment" from MDC to DPs, hence y variables are not defined for N0
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0));

    #######################
    #Define the variables.
    @variables(m,
            begin
               0 <= x_prev[i=1:Ni]
               0 <= x[i=1:Ni,t=1:(T-t0+1)] <= x_cap[i]
               0 <= f[i=1:N0,ii=1:Ni,t=1:(T-t0+1)] <= f_cap[i,ii]
               0 <= y[i=1:Ni,j=1:Nj]
               0 <= d[j=1:Nj]
               0 <= z[j=1:Nj]
            end
          );
    
    #######################
    #Define the objective.
    @objective(m,
               Min,
               sum(sum(sum(cb[i,ii,t+t0-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0) for t=1:(T-t0+1))
              +sum(sum(ch[i,t+t0-1]*x[i,t] for i=1:Ni) for t=1:(T-t0+1))
              +sum(sum(f[N0,i,t] for i=1:Ni)*h[t+t0-1] for t=1:(T-t0+1))
              +sum(sum(ca[i,j,T]*y[i,j] for j=1:Nj) for i=1:Ni)
              +sum(z[j] for j=1:Nj)*p
              +sum(x[i,T-t0+1]-sum(y[i,j] for j=1:Nj) for i=1:Ni)*q
               );
    
    #######################
    #Define the constraints.
    for i=1:Ni
        for t=1:(T-t0+1)
            if t == 1
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x_prev[i]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x_prev[i]
                            );
            else
                
                @constraint(m, 
                            x[i,t]
                           +sum(f[i,j,t] for j=1:Ni if j != i)
                           -sum(f[j,i,t] for j=1:N0 if j != i)
                           == x[i,t-1]
                            );
               @constraint(m, 
                            sum(f[i,j,t] for j=1:Ni if j != i)
                            <= x[i,t-1]
                            );                
            end
        end
       @constraint(m,
            sum(y[i,j] for j=1:Nj)
            <=x[i,T-t0+1]
           );
    end
    for j=1:Nj
        @constraint(m,
                    sum(y[i,j] for i=1:Ni)
                    <=d[j]
                   );
        @constraint(m,
                    d[j]-sum(y[i,j] for i=1:Ni)
                    <=z[j]
                   );
    end
    # Assuming no flow from MDC at the terminal stage
    @constraint(m,
        sum(f[N0,i,T-t0+1] for i=1:Ni) == 0
        );
    
    return m, x_prev, x, f, y, d, z
end


function rolling_MVP_eval_online(T,nbOS,OS_paths,OS_M,x0)
    start=time();
    costs = zeros(nbOS,T);
    Z = [];

    actions = [];
    procurmnt_amount = zeros(T);
    transport_amount = zeros(T);
    inventory_amount = zeros(T);
    dshortage_amount = 0;
    desalvege_amount = 0;
    
    #solve the first stage problem outside the loop because it's the same for all sample paths
    x1_initial = x0;
    s = 1;
    t0 = 1;
    imm_amount1 = 0;
    #define the determistic MVP problem at time t0 
    model_mvp, x_prev_mvp, x_mvp, f_mvp, y_mvp, d_mvp, z_mvp = deterministic_RH_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,t0);
    #fix the state variable value at time t0
    for i=1:Ni
        fix(x_prev_mvp[i], x1_initial[i]; force=true);
    end
    #which state are we currently in?
    current_k = OS_paths[s,t0];
    #what is the expected demand at time T
    d_bar = sum(sum(SCEN[kk][:,m]*(1/M) for m=1:M)*P_terminals[t0][current_k,kk] for kk=1:K);
    for j=1:Nj
        fix(d_mvp[j], d_bar[j]; force=true);
    end
   
    #solve the problem
    optimize!(model_mvp);
    status = termination_status(model_mvp);
    
    if status != MOI.OPTIMAL
        println(" in the out-of-sample evaluation at roll = ", t0);
        exit(0);
    else
        xval_running = value.(x_mvp);
        fval_running = value.(f_mvp);

        imm_amount1 = (sum(sum(cb[i,ii,t0]*fval_running[i,ii,1] for ii=1:Ni) for i=1:N0)
                         +sum(ch[i,t0]*xval_running[i,1] for i=1:Ni)
                         +sum(fval_running[N0,i,1] for i=1:Ni)*h[t0]); 
        x1_initial = xval_running[:,1];
        
        costs[s,t0] = imm_amount1
    end
    
    f_val1 = fval_running[N0,:,1];
    fval_all1 = fval_running[:,:,1];    
    
    
    for s=1:nbOS
        # starting the rolling horizon procedure one for each sample path
        println("Start sample path #", s); 
        x_initial = deepcopy(x1_initial);
        costs[s,1] = imm_amount1;
        
        #initialize the sample path cost
        path_amount = imm_amount1;
        
        
        temp_action = [];
        push!(temp_action,f_val1)
        
        procurmnt_amount[1] += (sum(f_val1[i] for i=1:Ni)*h[1])/nbOS
        transport_amount[1] += sum(sum(cb[i,ii,1]*fval_all1[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
        inventory_amount[1] += sum(ch[i,1]*x_initial[i] for i=1:Ni)/nbOS
        
        for t0=2:T
            #initialize the stage cost
            imm_amount = 0;
            
            #define the determistic MVP problem at time t0 
            model_mvp, x_prev_mvp, x_mvp, f_mvp, y_mvp, d_mvp, z_mvp = deterministic_RH_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,t0);
            
            #fix the state variable value at time t0
            for i=1:Ni
                fix(x_prev_mvp[i], x_initial[i]; force=true);
            end
            #which state are we currently in?
            current_k = OS_paths[s,t0];
            #what is the expected demand at time T
            d_bar = sum(sum(SCEN[kk][:,m]*(1/M) for m=1:M)*P_terminals[t0][current_k,kk] for kk=1:K);
            for j=1:Nj
                fix(d_mvp[j], d_bar[j]; force=true);
            end
            
            #solve the problem
            optimize!(model_mvp);
            status = termination_status(model_mvp);
            
            if status != MOI.OPTIMAL
                println(" in the out-of-sample evaluation at roll = ", t0);
                exit(0);
            else
                #collect values
                if t0<T
                    xval_running = value.(x_mvp);
                    fval_running = value.(f_mvp);
                    
                    imm_amount = (sum(sum(cb[i,ii,t0]*fval_running[i,ii,1] for ii=1:Ni) for i=1:N0)
                                     +sum(ch[i,t0]*xval_running[i,1] for i=1:Ni)
                                     +sum(fval_running[N0,i,1] for i=1:Ni)*h[t0]); 
                    x_initial = xval_running[:,1];
                    push!(temp_action,value.(f_mvp)[N0,:,1]);
                    
                    procurmnt_amount[t0] += (sum(fval_running[N0,i,1] for i=1:Ni)*h[t0])/nbOS
                    transport_amount[t0] += sum(sum(cb[i,ii,t0]*fval_running[i,ii,1] for ii=1:Ni) for i=1:N0)/nbOS
                    inventory_amount[t0] += sum(ch[i,t0]*xval_running[i,1] for i=1:Ni)/nbOS
                else
                    fval_T = value.(f_mvp)[:,:,(T-t0+1)]
                    xval_T = value.(x_mvp)[:,(T-t0+1)];                    
                    zval_T = value.(z_mvp);
                    yval_T = value.(y_mvp);
                    imm_amount = objective_value(model_mvp);
                    push!(temp_action,-zval_T);
                    
                    procurmnt_amount[t0] += (sum(fval_T[N0,i] for i=1:Ni)*h[t0])/nbOS
                    transport_amount[t0] += sum(sum(cb[i,ii,t0]*fval_T[i,ii] for ii=1:Ni) for i=1:N0)/nbOS
                    inventory_amount[t0] += sum(ch[i,t0]*xval_T[i] for i=1:Ni)/nbOS
                    
                    dshortage_amount += (sum(zval_T[j] for j=1:Nj)*p)/nbOS;
                    desalvege_amount += (sum(xval_T[i]-sum(yval_T[i,j] for j=1:Nj) for i=1:Ni)*q)/nbOS;
                end
                
                path_amount += imm_amount
                costs[s,t0] = imm_amount
            end
        end
        push!(Z,path_amount)
        push!(actions,temp_action)
    end
    
    COSTs = [procurmnt_amount, transport_amount, inventory_amount, dshortage_amount, desalvege_amount];

    elapsed = time() - start;
    
    UB_bar = mean(Z);
    UB_std = std(Z);
  
    UB_low = UB_bar-1.96*UB_std/sqrt(nbOS);
    UB_high = UB_bar+1.96*UB_std/sqrt(nbOS);

    println("μ = ", UB_bar);
    println("μ ± 1.96*σ/√NS = ", [UB_low,UB_high]);
    
    return costs, UB_bar, UB_low, UB_high, elapsed, actions, COSTs
end

function expected_number_of_stages(P_landfall)
    Q = P_landfall[1:end-1,1:end-1];
    M = inv(I-Q);
    exp_Ts = M*ones(size(M)[1]);
    exp_Ts = round.(Int,exp_Ts)
    return exp_Ts
end