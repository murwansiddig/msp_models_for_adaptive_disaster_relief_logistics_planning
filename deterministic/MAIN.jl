#change the currnet working directory to where the file is executed from
cwd = "/"*relpath((@__FILE__)*"/..","/");
cd(cwd) 
#comment this line if don't have bots setup
#include(homedir()*"/misc/tgtext.jl");

#1: instance number 
#2: number of supply points
#3: number of demand points
#4: number of stages
#5: number of scenarios in the out-of-sample
#6: initial state
#7: number of stages

PARAMS = ARGS;
inst = parse(Int, PARAMS[1][1:end]); #instance number 
Ni = parse(Int, PARAMS[2][1:end]); #number of supply points 
Nj = parse(Int, PARAMS[3][1:end]); #number of demand points
factor = parse(Float64, PARAMS[4][1:end]);
nbOS = parse(Int, PARAMS[5][1:end]); #number of sample paths used in out-of-sample
k_init = parse(Int, PARAMS[6][1:end]); #initial state
Tmax = parse(Int, PARAMS[7][1:end]); #T_max

include("packages.jl");
include("./functions/functions.jl");
include("./data/data.jl");

#create gurobi environment
const GRB_ENV = Gurobi.Env();

#Clairvoyance solution 
include("./policy/hurricane_CV.jl");

#Fully adaptive model
include("./policy/hurricane_FA.jl");
 
#Rolling-horizon 2SSP
include("./policy/hurricane_rolling_twostage.jl")

#Static 2SSP
include("./policy/hurricane_static_twostage.jl");


#tg_sendtext("Julia: $instname is DONE!"); #comment this line if don't have bots setup

