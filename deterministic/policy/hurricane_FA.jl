function hurricane_FA()
    #Fully adaptive
    #Defne the model 
    model, x, f, theta, y, d, z, FB1, FB2 = define_models(K,T,Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q);

    #train the model
    k = copy(k_init)
    LB, train_time, iter = train_models_offline(model,x,f,theta,y,d,z,FB1,FB2,k,T);
    println("***********************************************")
    println("***********************************************")
    println("finished training!")
    println("LB = ", LB[end])
    println("training time = ", train_time)
    println("number of iterations = ", iter)


    #evaluate the model 
    costs, UB_bar, UB_low, UB_high, elapsed, actions, COSTs = FOSDDP_eval_offline(model,x,f,theta,y,d,z,FB1,FB2,k,T,nbOS);

    fname = "./output/benchmark/FAresults.csv"
    df = CSV.read(fname,DataFrame);
    reSults = Matrix(df);
    reSults[inst,1] = LB[end]
    reSults[inst,2] = UB_bar
    reSults[inst,3] = UB_bar-UB_low
    reSults[inst,4] = train_time
    reSults[inst,5] = elapsed
    reSults[inst,6] = iter

    updf = DataFrame(reSults, :auto);
    CSV.write(fname,updf)
    
end

@time hurricane_FA();
println("############################################################")
println("############################################################")