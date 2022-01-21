function hurricane_CV()
    #clairvoyant solution 

    #define the clairvoyant model 
    model_cv, x_cv, f_cv, y_cv, d_cv, z_cv = deterministic_model(Ni,Nj,N0,x_cap,cb,ch,h,ca,p,q,T,x_0);

    #evaluate the clairvoyant solution
    costs_cv, cv_bar, cv_low, cv_high, elapsed_cv =  clairvoyant_eval(model_cv,x_cv,f_cv,y_cv,d_cv,z_cv,T,nbOS);

    if inst == 1
        sleep(50)
    elseif inst == 2
        sleep(5)
    elseif inst == 3
        sleep(10)
    elseif inst == 4
        sleep(20)
    elseif inst == 5
        sleep(25)
    elseif inst == 6
        sleep(30)
    elseif inst == 7
        sleep(35)
    elseif inst == 8
        sleep(40)
    elseif inst == 9
        sleep(45)
    end
        
    fname = "./output/benchmark/CVresults.csv"
    df = CSV.read(fname,DataFrame);
    reSults = Matrix(df);
    reSults[inst,1] = 0
    reSults[inst,2] = cv_bar
    reSults[inst,3] = cv_bar-cv_low
    reSults[inst,4] = 0
    reSults[inst,5] = elapsed_cv
    #reSults[inst,6] = 0

    updf = DataFrame(reSults, :auto);
    CSV.write(fname,updf)
    
end

@time temp_vec = hurricane_CV();
println("############################################################")
println("############################################################")