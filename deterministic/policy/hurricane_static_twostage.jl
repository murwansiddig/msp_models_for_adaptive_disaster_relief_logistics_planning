function hurricane_static_twostage()
    #static two-stage model
    
    #train the model
    LB_twoSP, thetaval_twoSP, iter_twoSP, xval_twoSP, train_Elapsed_twoSP = two_stage_static_train();
    println("***********************************************")
    println("finished training!")
    println("LB = ", LB_twoSP)
    println("training time = ", train_Elapsed_twoSP)
    println("number of iterations = ", iter_twoSP)
    
    #evaluate the model
    costs_twoSP, UB_twoSP_bar, UB_twoSP_low, UB_twoSP_high, eval_Elapsed_twoSP, actions, COSTs = two_stage_static_eval(xval_twoSP);

    fname = "./output/benchmark/static2SPresults.csv"
    df = CSV.read(fname,DataFrame);
    
    reSults = Matrix(df);
    reSults[inst,1] = LB_twoSP[end]
    reSults[inst,2] = UB_twoSP_bar
    reSults[inst,3] = UB_twoSP_bar-UB_twoSP_low
    reSults[inst,4] = train_Elapsed_twoSP
    reSults[inst,5] = eval_Elapsed_twoSP
    reSults[inst,6] = iter_twoSP
    updf = DataFrame(reSults, :auto);
    CSV.write(fname,updf)
end

@time hurricane_static_twostage();
println("############################################################")
println("############################################################")