#static 2SSP models

#Define first-stage master problem
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
               -1e10 <= θ 
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
           #@constraint(m, v[i]<=x[i,t])
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
            #=   
                @constraint(m, x[i,t-1]
                             -x[i,t]
                             -sum(f[i,j,t] for j=1:Ni if j != i)
                             +sum(f[j,i,t] for j=1:N0 if j != i)==0);
            =#   
            end
        end
    end
    
    #in case the first stage corresponds to the period where the hurricane make landfall
    for j=1:Nj
        @constraint(m, z[j]+sum(y[i,j] for i=1:Ni) >= ξ[j]);
    end
    
    return m, x, f, y, z, v, θ 
    
end


###############################################################
###############################################################

#Define second-stage scenario supbproblem
function RH_2SSP_second_stage(t_roll,nbstages1,nbstages2)
    
    #######################
    #Define the model.
    m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GRB_ENV), "OutputFlag" => 0, "Presolve" => 0));

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


###############################################################
###############################################################

#define a feasibility problem
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
            ConsFB[t,i] = @constraint(m, sum(y[i,j,t] for j=1:Nj)+v[i,t]+a_plus[i,t]-a_minus[i,t]== 0);         
        end
        for j=1:Nj
           dCons[t,j] = @constraint(m, z[j,t]+sum(y[i,j,t] for i=1:Ni) >= 0);
        end
        rCons[t] = @constraint(m, reimbursement[t] == 0)
    end

    return m, y, z, v, ConsFB, dCons, rCons
end


###############################################################
###############################################################

#defines the two-stage SP models
function RH_2SSP_define_models(t_roll,nbstages1,nbstages2,x_init,ξ)
    #define first stage (master problem) model  
    master, x, f, y1, z1, v1, θ = RH_2SSP_first_stage(t_roll,nbstages1,x_init,ξ)
    
    #define second stage (subproblem) optimality model
    subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons = RH_2SSP_second_stage(t_roll,nbstages1,nbstages2)
    
    #define second stage (subproblem) feasability model
    feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons = feasability_problem(t_roll,nbstages1,nbstages2)
    
    
    return master, x, f, y1, z1, v1, θ, subproblem, y2, z2, v2, r2, ConsFB, dCons, rCons, feas_m, feas_y, feas_z, feas_v, feas_ConsFB, feas_dCons, feas_rCons
end

###############################################################
###############################################################

#initialize parameter for two-stage model
function initialize(nbstages1,nbstages2,s,t_roll)
    LB = -1e10; UB = 1e10; iter = 0; θval = 0;
    xval = zeros(Ni,nbstages1);fval = zeros(N0,Ni,nbstages1);
    start=time();
    Elapsed = time()-start;
    allscen = Matrix(CSV.read("./data/OOS.csv",DataFrame));
    allscenM = Matrix(CSV.read("./data/inOOS.csv",DataFrame))
    scen = allscen[collect(1:convert(Int,10000/nbscen):10000),1:T]
    scenM = allscenM[collect(1:convert(Int,10000/nbscen):10000),1]

    if t_roll > 1
        for n=1:nbscen
            scenM[n] = rand(1:M)
            
            for t=2:t_roll
                scen[n,t] = allscen[s,t]
            end
            
            for t=t_roll+1:T
                scen[n,t] = MC_sample(scen[n,t-1]);
            end
        end
    end
    qprob = fill(1/nbscen,nbscen)
    return LB, UB, iter, xval, fval, θval, start, Elapsed, nbscen, scen, scenM, qprob
end

###############################################################
###############################################################

#solves the two-stage SP model
function RH_2SSP_solve_roll(s,t_roll,nbstages1,nbstages2,master,subproblem,x,f,θ,y,z,v,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons)
    
    LB, UB, iter, xval, fval, θval, start, Elapsed, nbscen, scen, scenM, qprob = initialize(nbstages1,nbstages2,s,t_roll)
    
    while (UB-LB)*1.0/max(1e-10,abs(LB)) > ϵ && Elapsed < time_limit 
        iter+=1;
        # solve first stage
        LB, xval, fval, θval = solve_first_stage(LB,xval,fval,θval,master,subproblem,x,f,θ);
        if nbstages2 > 0
            # solve second stage 
            flag, Qbar = solve_second_stage(t_roll,nbstages1,nbstages2,xval,fval,θval,scen,scenM,nbscen,qprob,master,subproblem,x,f,θ,y,z,v,ConsFB,dCons,rCons,feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons)
            if flag != -1
                UB = min(LB-θval+Qbar,UB);
            end
            Elapsed = time()-start;
        else
            println("exiting here")
            break
        end
    end

    return LB, UB, iter, xval, fval, θval, Elapsed
end

###############################################################
###############################################################

#solves the first-stage problem
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

###############################################################
###############################################################

#solves the second-stage problem
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
        if τ === nothing     
            RH_2SSP_update_RHS(τ,scen[n,end],scenM[n],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        else
            RH_2SSP_update_RHS(τ,scen[n,τ],scenM[n],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        end
        
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
                             -sum(pi2bar[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*f[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*x[i,t] for i=1:Ni)+sum(f[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t_roll+t-1>τ)
                            >= 
                         Qbar-sum(pi1bar[t-1,i]*(xval[i,t-1]-xval[i,t]-sum(fval[i,j,t] for j=1:Ni if j != i)+sum(fval[j,i,t] for j=1:N0 if j != i)) for i=1:Ni, t=2:nbstages1)
                             -sum(pi2bar[t-1]*(sum(sum(cb[i,ii,t_roll+t-1]*fval[i,ii,t] for ii=1:Ni) for i=1:N0)+sum(ch[i,t_roll+t-1]*xval[i,t] for i=1:Ni)+sum(fval[N0,i,t] for i=1:Ni)*h[t_roll+t-1]) for t=2:nbstages1 if t_roll+t-1>τ)
                            );
            end

            flag = 1;
        end
    end
    return flag, Qbar
end


###############################################################
###############################################################
#solves scenario subproblem of the second stage
function solve_scen_subproblem(Q,pi1,pi2,n,subproblem,ConsFB,dCons,rCons,nbstages2)
    flag = 0
    optimize!(subproblem) #solve the model
    status_subproblem = termination_status(subproblem); #check the status 
    if status_subproblem != MOI.OPTIMAL
        #println("Subproblem problem status is: ", status_subproblem, " === Oops! :/")
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

###############################################################
###############################################################

#updates the RHS of the flow-balance and demand constraints 
function RH_2SSP_update_RHS(τ,k_t,m,nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
    for t=2:nbstages1
        for i=1:Ni
            set_normalized_rhs(ConsFB[t-1,i],xval[i,t-1]
                         -xval[i,t]
                         -sum(fval[i,j,t] for j=1:Ni if j != i)
                         +sum(fval[j,i,t] for j=1:N0 if j != i)
                         );         
        end
        
        for j=1:Nj
            if t_roll+t-1 == τ && k_t ∉ absorbing_states           
                set_normalized_rhs(dCons[t-1,j], SCEN[k_t][j,m]);
            else
                set_normalized_rhs(dCons[t-1,j], 0);
            end
        end
        
        if τ === nothing || t_roll+t-1 <= τ
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

###############################################################
###############################################################

#generate a feasability cut
function generate_feasability_cut(feas_m,feas_y,feas_z,feas_v,feas_ConsFB,feas_dCons,feas_rCons,n,scen,scenM,nbstages1,xval,fval,t_roll,nbstages2,master,x,f)
    #identify the period where the hurricane makes landfall 
    τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, scen[n,:]);

    #update the RHS
    if τ === nothing
        RH_2SSP_update_RHS(τ,scen[n,end],scenM[n],nbstages1,feas_ConsFB,feas_dCons,feas_rCons,xval,fval,t_roll)
    else
        RH_2SSP_update_RHS(τ,scen[n,τ],scenM[n],nbstages1,feas_ConsFB,feas_dCons,feas_rCons,xval,fval,t_roll)
    end

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