#miscellaneous functions 

#sample a markovian state
function MC_sample(current_state)
    states = collect(1:1:K);
    k = sample(states, Weights(P_joint[current_state,:]));
    return k
end


#training termination check
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


# function to create an out-of-sample for given initial state k_init
function create_OSpaths(k_init)
    OS_paths = Matrix(CSV.read("OOS.csv",DataFrame)); #read the out-of-sample file
    rr = size(OS_paths)[1];
    cc = size(OS_paths)[2];
    OS_paths = fill(k_init,rr,cc);
    for s=1:rr
        for t=2:cc
            OS_paths[s,t] = MC_sample(OS_paths[s,t-1]) 
        end
    end
    df = DataFrame(OS_paths, :auto);
    CSV.write("OOS.csv",df)
end

#function to save lp files
function save_lp(lp,name)
    open("./output/lp/$name.lp", "w") do f
        print(f, lp)
    end
end

#function save an n×m matrix as an csv file
function save_csv(m_data, m_colname, f_dir, m_fname)
    fname = f_dir*m_fname*".csv"
    df = DataFrame(m_data, m_colname);
    CSV.write(fname,df)
end