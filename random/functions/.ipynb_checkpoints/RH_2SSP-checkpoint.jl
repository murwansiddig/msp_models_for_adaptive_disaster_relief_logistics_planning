# #evaluate the two-stage SP in a rolling-horizon fashion

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

    OS_paths = Matrix(CSV.read("./data/OOS.csv",DataFrame)); #read the out-of-sample file
    OS_M = Matrix(CSV.read("./data/inOOS.csv",DataFrame))[:,1] #read the second layer OOS
    
    optimize!(master) #solve the model
    #update the values
    xval = value.(x);
    fval = value.(f);
    xval1 = deepcopy(xval[:,1])
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
    objs_RH2SSP = zeros(nbOS,T);
    objs_RH2SSP[:,1] .= temp_obj
    for s=1:nbOS
        println("s = ", s)
        τ = findfirst(x -> S[x][3] == Nc-1 && x ∉ absorbing_states, OS_paths[s,1:T]);
        x_init = deepcopy(xval1[:,1]);
        for t_roll=2:T-1
            if τ !== nothing && t_roll <= τ
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
                x_init = deepcopy(xval[:,1]);
                y1val = value.(y1);
                z1val = value.(z1);
                v1val = value.(v1);
                objs_RH2SSP[s,t_roll] = (sum(sum(cb[i,ii,t_roll]*fval[i,ii,1] for ii=1:Ni) for i=1:N0) 
                                         +sum(ch[i,t_roll]*xval[i,1] for i=1:Ni)
                                         +sum(fval[N0,i,1] for i=1:Ni)*h[t_roll]
                                         +sum(sum(ca[i,j,t_roll]*y1val[i,j] for j=1:Nj) for i=1:Ni)
                                         +sum(z1val[j] for j=1:Nj)*p
                                         +sum(v1val[i] for i=1:Ni)*q
                                        );

            else
                continue 
            end
        end
    end
    ########


    ##########
    RH2SSP_bar = mean(sum(objs_RH2SSP[:,t] for t=1:T));
    RH2SSP_std = std(sum(objs_RH2SSP[:,t] for t=1:T));
    RH2SSP_low = RH2SSP_bar-1.96*RH2SSP_std/sqrt(nbOS);
    RH2SSP_high = RH2SSP_bar+1.96*RH2SSP_std/sqrt(nbOS);
    println("μ ± 1.96*σ/√NS = ", RH2SSP_bar, " ± ", [RH2SSP_low,RH2SSP_high]);
    elapsed = time() - start;
  
    return objs_RH2SSP, RH2SSP_bar, RH2SSP_low, RH2SSP_high, elapsed
    
end