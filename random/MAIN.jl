#change the currnet working directory to where the file is executed from
cwd = "/"*relpath((@__FILE__)*"/..","/");
cd(cwd) 
#comment this line if don't have bots setup
#include(homedir()*"/misc/tgtext.jl");

#1: instance number 
#2: number of supply points
#3: number of demand points
#4: cost factor
#5: number of scenarios in the out-of-sample 
#6: initial state
#7: number of sample paths when solving two-stage
PARAMS = ARGS;    
if length(ARGS) < 2
    PARAMS = ["3", "2", "4", "5", "10", "65", "100"];
end;

instname = PARAMS[1]*"-"*PARAMS[2]*"SPs-"*PARAMS[3]*"DPs-"*PARAMS[4]*"ρ";
#tg_sendtext("Julia: Now starting instance $instname"); #comment this line if don't have bots setup
inst = parse(Int, PARAMS[1]); #instance number 
Ni = parse(Int, PARAMS[2]); N0=Ni+1; #number of supply points 
Nj = parse(Int, PARAMS[3]); #number of demand points
factor = parse(Float64, PARAMS[4]); #cost scaling factor
nbOS = parse(Int, PARAMS[5]); #cost scaling factor
k_init = parse(Int, PARAMS[6]); #cost scaling factor
nbscen = parse(Int, PARAMS[7]); #number of sample paths when solving two-stage

include("packages.jl");
include("./data/data.jl");
include("./functions/functions.jl");
#create gurobi environment
const GRB_ENV = Gurobi.Env();

#Clairvoyance solution 
#include("./policy/hurricane_CV.jl");

#Fully adaptive model
include("./policy/hurricane_FA.jl");
 
#Rolling-horizon 2SSP
#include("./policy/hurricane_rolling_twostage.jl")

#Static 2SSP
#include("./policy/hurricane_static_twostage.jl");

#sensitivity analysis
#include("SENS.jl");


#tg_sendtext("Julia: $instname is DONE!"); #comment this line if don't have bots setup
