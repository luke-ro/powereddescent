using OrdinaryDiffEq
using LinearAlgebra 
using Plots
using Convex
using ECOS


const e_sig = [0.0; 0; 0; 1];
const E = [Matrix(1.0I,6,6) zeros(Float64,6,1)];
const F = [zeros(Float64,1,6) 1];
const E_u = [Matrix(1.0I,3,3) zeros(Float64,3,1)]

struct ProblemParameters
    g
    alpha
    m_dr0
    m_wet
    I_sp
    T_bar
    T_1
    T_2
    n
    phi
    y0
    S #######################################
    v # General constraint funciton parameters
    c #
    a #######################################
end

"""
calcAB(dt,alpha)

calculate the A and B matrices for a given timestep Dt and 
mass consumption rate constant alpha
"""
function calcAB(Dt,alpha)
    Ac = zeros(Float64,7,7);
    Ac[1:3,4:6] = Matrix(1.0I, 3, 3);
    
    Bc = zeros(Float64, 7, 4);
    Bc[4:6,1:3] = Matrix(1.0I, 3, 3);
    Bc[7,4] = alpha;

    function intB!(dB, B, p, t)
        temp = exp(Ac*(Dt-t))*Bc;
        dB[:] = reshape(temp,:,1);
        # println(dB)
 
    end

    A = exp(Ac*Dt);
    prob = ODEProblem(intB!,zeros(Float64,28,1),(0,Dt));
    B = solve(prob, Tsit5()); 
    # print(size(B))
    B = reshape(B[:,1,end],7,4);

    return A,B;
end

function makePHI(A,k)
    return A^k;
end

function makePSI(A,B,k,N;m=7,n=4)
    # PSI = zeros(Float64,m*N,n*N);
    # for i in 2:(k+1)
    #     for j in 1:i-1
    #         PSI[(i-1)*m+1:i*m, (j-1)*n+1:j*n] = A^(i-2-(j-1))*B; 
    #     end
    # end
    println("In makePsi:")
    print("N: "); println(N);
    print("n: "); println(n);
    PSI = zeros(Float64,m,n*(N+1)) # number of columns must match number of rows of eta
    for i in 1:k
        PSI[:,(i-1)*n+1:i*n] = A^(k-1-(i-1))*B; 
    end
    print("size(psi): "); println(size(PSI))
    return PSI;
end

function makeLambda(A,B,k;m=7,n=4)

    LAMBDA = zeros(Float64,m,n);
    if k==0
        return LAMBDA
    end

    for i in 1:k
        LAMBDA += A^(i-1)*B;
    end
    return LAMBDA;
end

function makeUpsilon(N,k;p=4)
    UPSILON = zeros(Float64,p,p*(N+1));
    UPSILON[:,(k)*p+1:(k+1)*p] = Matrix(1.0I, p, p); #TODO what is this when k==0?
    # UPSILON[:,(k-1)*m+1:k*m] = I;
    return UPSILON
end

function makeomega(w,N)
    omega = w*repeat(e_sig,outer=N+1)
    return omega
end

function makeConstraintEq50!(constraints,N,eta,E_u,e_sig;m=7,n=4)
    # UPSILON = makeUpsilon(0)
    # constraints = Constraint[sumsquares(E_u*UPSILON*eta)<=transpose(e_sig)*UPSILON*eta]
    for k in 0:N
        UPSILON = makeUpsilon(N,k)
        push!(constraints,sumsquares(E_u*UPSILON*eta)<=transpose(e_sig)*UPSILON*eta)
    end
    return
end

function makez0(t,ps::ProblemParameters)
    rho2 = ps.n*ps.T_2*cos(ps.phi) #eq3
    return log(ps.m_wet-ps.alpha*rho2*t)
end

function makeMu(num,k,Dt,ps::ProblemParameters)
    t = Dt*k;
    if num==1
        rho = ps.n*ps.T_1*cos(ps.phi) #eq3
    elseif num==2
        rho = ps.n*ps.T_2*cos(ps.phi) #eq3
    else
        throw("Unrecognized number for Mu in MakeMu")
    end
    z0 = makez0(t,ps) #eq35
    mu = rho*exp(-z0) #eq33
    return mu
end

"""

"""
function addConstraints!(constraints,A,B,N,Dt,eta,omega,params::ProblemParameters;m=7,n=4)
    n_s = size(a,1) #number of constraints for ||S*x-v||+c^T*x + a <= 0

    #constraint in eq50
    for k in 0:N
        UPSILON = makeUpsilon(N,k)
        push!(constraints,sumsquares(E_u*UPSILON*eta)<=transpose(e_sig)*UPSILON*eta)
    end

    for k in 1:N
        t = k*Dt
        mu1_k = makeMu(1,k,Dt,params)
        mu2_k = makeMu(2,k,Dt,params)
        PSI_k = makePSI(A,B,k,N)
        z0_k = makez0(t,params)
        UPSILON_k = makeUpsilon(N,k)
        PHI_k =  makePHI(A,k)
        LAMBDA_k = makeLambda(A,B,k)
        XI_k = PHI_k*params.y0 + LAMBDA_k*vcat(params.g, [0])

        println("here:")
        # print("size(XI_k+PHI_k*eta): "); println(size(XI_k+PHI_k*eta))
        println(N)
        print(size(XI_k+(PSI_k*eta)))
        
        # eq51 part 1
        temp = mu1_k*(1-((F*(XI_k+PSI_k*eta)-z0_k))+((F*(XI_k+PSI_k*eta)-z0_k)*(F*(XI_k+PSI_k*eta)-z0_k))/2) <= transpose(e_sig)*UPSILON_k*eta
        push!(constraints,temp)
        # eq51 part 2
        temp = transpose(e_sig)*UPSILON_k*eta <= mu2_k*(1-(F*(XI_k+PSI_k*eta)-z0_k))
        push!(constraints,temp)
        
        # eq52 part 1
        rho1 = params.n*params.T_1*cos(params.phi)
        rho2 = params.n*params.T_2*cos(params.phi)
        temp = log(params.m_wet-params.alpha*rho2*t) <= F*(XI_k+PSI_k*eta)
        push!(constraints,temp)

        #eq 52 part 1
        temp = F*(XI_k+PSI_k*eta) <= log(params.m_wet-params.alpha*rho1*t)
        push!(constraints,temp)

        #||S*x-v||+c^T*x + a <= 0 constraints
        for j = 1:n_s
            temp = sumsquares(params.S[j]*E*(XI_k+PSI_k*eta)-params.v[j]) + transpose(params.c[j])*E*(XI_k+PSI_k*eta) + params.a[j] <= 0
            push!(constraints,temp)
        end

        return 
    end
end

function solveProblem(N,Dt,params)
    A,B = calcAB(Dt,params.alpha);
    m = size(B,2);
    n = size(B,1);

    # N = Variable();
    eta = Variable(m*(N+1),1);

    w = Dt;
    omega = makeomega(w,N)#probably need to turn this into a 1d Vector

    constraints = []
    addConstraints!(constraints,A,B,N,Dt,eta,omega,params)

    # for k = 0:N
    #     UPSILON = makeUpsilon(k);
    #     E_u*UPSILON<=transpose(e_sig)*UPSILON*
    # end



    # for i in 0:N
    #     UPSILON = makeUpsilon(k);
    #     PHI = makePHI(A,k);
    #     PSI = makePSI(A,B,k);
    #     LAMBDA = makeLambda(A,B,k);
    # end

    cost = transpose(omega)*eta;
    problem = minimize(cost,constraints)
    solver = ECOS.Optimizer

    Convex.solve!(problem, solver, silent=true)

    eta_opt = evaluate(eta)
    # N_opt = evaluate(N)
    
    return eta_opt, A, B
end

function eom(t,x,u,A,B,ps::ProblemParameters)
    dx = A*x + B*(vcat(ps.g,[0]) + u)
    return dx
end

function calcTrajectory(Dt,N,A,B,eta,ps::ProblemParameters)
    x = []
    u = []
    for k = 1:N
        t = k*Dt
        mu1_k = makeMu(1,k,Dt,ps)
        mu2_k = makeMu(2,k,Dt,ps)
        PSI_k = makePSI(A,B,k,N)
        z0_k = makez0(t,ps)
        UPSILON_k = makeUpsilon(N,k)
        PHI_k =  makePHI(A,k)
        LAMBDA_k = makeLambda(A,B,k)
        XI_k = PHI_k*ps.y0 + LAMBDA_k*vcat(ps.g, [0])

        state = XI_k+PSI_k*eta
        control = UPSILON_k*eta
            

        push!(x,state[1:6])
        push!(u,control[1:3])
    end
    return reduce(hcat,x),u
end



Dt = 0.5;
N = 10

g = [-3.7114; 0; 0];
alpha = 0.25;
m_dry = 1505; #kg
m_wet = 1905; #kg
I_sp = 225; #s
T_bar = 3.1; #kN
T_1 = 0.3*T_bar;
T_2 = 0.8*T_bar;
n = 6;
phi = 27*3.1415/180;
r0 = [1.5,0,2]*1000 #km -> m
dr0 = [-75,0,100] #m/s
y0 = vcat(r0,dr0,log(m_wet))
S = []
v = []
c = []
a = []
params = ProblemParameters(g,alpha,m_dry,m_wet,I_sp,T_bar,T_1,T_2,n,phi,y0,S,v,c,a)
print(params)




eta_opt, A, B = solveProblem(N,Dt,params)
# [t,x] = simulateProblem(A,B,)

eta = reshape(eta_opt,:,4)
x,u = calcTrajectory(Dt,N,A,B,eta_opt,params::ProblemParameters)

print(x)

p = plot(transpose(x))
display(p)
# plot(eta[:,1:3])
