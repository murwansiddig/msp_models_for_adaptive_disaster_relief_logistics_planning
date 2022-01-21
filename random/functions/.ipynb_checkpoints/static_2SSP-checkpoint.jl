#evaluate the two-stage SP in a static fashion

function static_2SSP_evaluation()
    start=time();
    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    optimize!(master) #solve the model
    #update the values
    xval = value.(x);
    fval = value.(f);
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
    objs_st2SSP = fill(temp_obj,nbOS,T);
    for s=1:nbOS
        #identify the period where the hurricane makes landfall 
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
        
        #update the RHS
        if τ === nothing
            RH_2SSP_update_RHS(τ,OS_paths[s,end],OS_M[s],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        else
            RH_2SSP_update_RHS(τ,OS_paths[s,τ],OS_M[s],nbstages1,ConsFB,dCons,rCons,xval,fval,t_roll)
        end
        
        Q = zeros(nbOS); #list for all the optimal values
        #list for all the dual multiplies of the first and second set of constraints
        pi1 = Array{Any,1}(undef,nbOS); 
        pi2 = Array{Any,1}(undef,nbOS); 
        #solve the subproblem and store the dual information
        Q, pi1, pi2, flag = solve_scen_subproblem(Q,pi1,pi2,s,subproblem,ConsFB,dCons,rCons,nbstages2)
        
        yval = value.(y2);
        zval = value.(z2);
        vval = value.(v2);
        rval = value.(r2);
       
        

        
        for t=1:nbstages2
            objs_st2SSP[s,t_roll+t] = (sum(sum(cb[i,ii,t_roll+t]*fval[i,ii,t_roll+t] for ii=1:Ni) for i=1:N0) 
                                     +sum(ch[i,t_roll+t]*xval[i,t_roll+t] for i=1:Ni)
                                     +sum(fval[N0,i,t_roll+t] for i=1:Ni)*h[t_roll+t]
                                     +sum(sum(ca[i,j,t_roll+t]*yval[i,j,t] for j=1:Nj) for i=1:Ni)
                                     +sum(zval[j,t] for j=1:Nj)*p
                                     +sum(vval[i,t] for i=1:Ni)*q
                                     +rval[t]
                                    );
        end
    end
      
    st2SSP_bar = mean(sum(objs_st2SSP[:,t] for t=1:T));
    st2SSP_std = std(sum(objs_st2SSP[:,t] for t=1:T));
    st2SSP_low = st2SSP_bar-1.96*st2SSP_std/sqrt(nbOS);
    st2SSP_high = st2SSP_bar+1.96*st2SSP_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", st2SSP_bar, " ± ", [st2SSP_low,st2SSP_high]);
    elapsed = time() - start;
    
    return objs_st2SSP, st2SSP_bar, st2SSP_low, st2SSP_high, elapsed#, vals
    
end