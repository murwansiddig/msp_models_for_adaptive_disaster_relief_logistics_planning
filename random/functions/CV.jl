#clairvoyance model functions 

#Define the model
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

###############################################################
###############################################################

#evaluate the model
function clairvoyant_eval()
    start=time();
    
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    objs_cv = zeros(nbOS);
    
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
             
            optimize!(m_cv); #solve the model        
            status = termination_status(m_cv); #check the status 
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
            end
        end
    end

    cv_bar = mean(objs_cv);
    cv_std = std(objs_cv);
    cv_low = cv_bar-1.96*cv_std/sqrt(nbOS);
    cv_high = cv_bar+1.96*cv_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", cv_bar, "±", [cv_low,cv_high]);
    elapsed = time() - start;
    vals = [xval_cv, fval_cv, yval_cv, zval_cv, vval_cv];

    return objs_cv, cv_bar, cv_low, cv_high, elapsed#, vals
end