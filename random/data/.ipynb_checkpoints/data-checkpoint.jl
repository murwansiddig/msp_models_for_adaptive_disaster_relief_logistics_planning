Random.seed!(2021);
nodes = CSV.read("./data/nodes.csv",DataFrame); # read the nodes locations
x_cap = Matrix(nodes)[1:Ni,5]*(Nj/Ni); #capacity of each SP
x_0 = zeros(Ni); #initial items at different SPs
f_cap = fill(Inf,N0,Ni);

#locations of the different facilities:
MDC = [350,450];
x_low = 0; x_up = 700; x_int = 100;
y_low = 0; y_up = 100;
L = [(r,r+x_int) for r in collect(x_low:x_int:x_up)];

# probability distributions:
P_intesity = Matrix(CSV.read("./data/intensity.csv",DataFrame)) #intensity MC
P_location = Matrix(CSV.read("./data/location.csv",DataFrame)) #location MC
P_landfall = Matrix(CSV.read("./data/landfall_7.csv",DataFrame)) #landfall MC

Na = size(P_intesity)[1]; #intensity MC number of states
Nb = size(P_location)[1]; #location MC number of states
Nc = size(P_landfall)[1]; #landfall MC number of states
T = copy(Nc); #define T_max

K = Na*Nb*Nc; #number of state in the joint MC
P_joint = zeros(K,K); #initialize the joint probability distribution MC
S = Array{Any,1}(undef,K); #list with elements [intensity,location]
absorbing_states = []; # list of absorbing states

k1 = 0; #counter for the number of states
for k=1:Na, l=1:Nb, f=1:Nc
    global k1 +=1
    k2 = 0;
    for n=1:Na, m=1:Nb, j=1:Nc
        k2 +=1
        P_joint[k1,k2] = P_intesity[k,n]*P_location[l,m]*P_landfall[f,j];
    end
    S[k1] = [k,l,f]; 
    if k == 1 || f == Nc
        push!(absorbing_states,k1)
    end
end

#normalize the probabilities
P_temp = deepcopy(P_joint);
for k=1:K, kk=1:K
    P_joint[k,kk] = P_temp[k,kk]/sum(P_temp[k,:])
end
    
# Create the transition probability from stage 1 to stage T (applying the C-K equation)
P_terminals = Matrix{Float64}[]
for t = 1:T  
    P_terminal = P_joint^(T-t);
    P_terminal_c = deepcopy(P_terminal);
    #normalize the probabilities
    for k=1:K, kk=1:K
       P_terminal[k,kk] = P_terminal_c[k,kk]/sum(P_terminal_c[k,:])
    end
    push!(P_terminals, P_terminal);
end    

#list for the coordinates of the different supply points
SP = []; 
for i=1:Ni
    push!(SP,[nodes[i,1],nodes[i,2]+100]);
end

#list for the coordinates of the different demand points
DP = []; 
for j=1:Nj
    push!(DP,[nodes[j,3],nodes[j,4]]);
end

#fuel cost
fuel = 0.0038;

#unit cost of transporting/rerouting items from MDC/SP i to/between SP i' 
cb = Array{Float64,3}(undef,N0,Ni,T);
for i=1:N0, ii=1:Ni, t=1:T
    if i < N0
        cb[i,ii,t] = fuel*norm(SP[i]-SP[ii],2)*(1+factor*(t-1))
    else
        cb[i,ii,t] = fuel*norm(MDC-SP[ii],2)*(1+factor*(t-1))
    end 
end

#unit cost of transporting items from MDC/SP i to/between a demand point j
ca = Array{Float64,3}(undef,N0,Nj,T);
for i=1:N0, j=1:Nj, t=1:T
    if i < N0
        ca[i,j,t] = fuel*norm(SP[i]-DP[j],2)*(1+factor*(t-1))
    else
        ca[i,j,t] = fuel*norm(MDC-DP[j],2)*(1+factor*(t-1))
    end 
end

base = 5; # base unit cost for logistic costs
h = Array{Float64,1}(undef,T); #unit cost for purchasing a relief item
ch = Array{Float64,2}(undef,Ni,T); #unit cost for holding an item at SP i
for t=1:T
    h[t] = base*(1+factor*(t-1));
    ch[:,t] = fill(0.2*base,Ni);
end
p = 80*base; #penalty cost for failing to serve the demand
q = -0.05*base; #salvage cost for each unit of overstock.
D_max = 400; #the maximum demand that can happen


#split the sample evenly towards the left and the right of the range
M = 10;
layer2 = []; 
for l in L
    cover = l;
    list = [];
    right = [cover[1],cover[2]]
    parts_right = (right[2]-right[1])/(2*M);
    count_right = right[1];
    count_2 = 0
    while true    
        count_2+=1
        count_right +=parts_right
        if count_right >= right[2]
            break
        elseif isodd(count_2)
            push!(list,count_right)
        end
    end
    push!(layer2,list);
end

SCEN = Array{Any,1}(undef,K); #SCEN is the list for each state k 
c_max = 300; #the largest possible radius of a hurricane

# mapping state to demand
for k=1:K
    #initialize a scenario matrix (scen) for each DP j and each possible point m in layer2;
    scen = zeros(Nj,M); #we will create a scenario list for each state k
    a = S[k][1]; #what is the observed intensity state
    l = S[k][2]; #what is the observed location state
    
    #what are the coordinates of where the hurricane made landfall [xx,yy]
    #note that yy=0 since the huricane makes landfall in the cost by assumption
    #since we observed state location l, here are the M possible locations
    for m=1:M
        predicted = layer2[l][m]; #this is the predicted x_coordinates for the landfall location
        xx_coord = max(0,min(predicted,x_up)); #we already now x in can't be smaller than 0 and bigger than 325; 
        landfall = [xx_coord,0]
        #now lets calculate the demand from each DP to the m location 
        for j=1:Nj
            #how far did destination to DPj j
            c_j = norm(landfall-DP[j],2);
            if c_j <= c_max
                scen[j,m] = D_max*(1-(c_j/c_max))*(a-1)^2/((Na-1)^2)
            else
                scen[j,m] = 0;
            end
        end
    end
    SCEN[k] = scen;
end

##############################
max_iter = 100000;
stall = 500;
cutviol_maxiter = 100000;
nbhrs = 3;
time_limit = nbhrs*60^2;
#time_limit = nbhrs*60^5;
ϵ = 1e-5;
