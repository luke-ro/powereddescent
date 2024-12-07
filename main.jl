using OrdinaryDiffEq
using LinearAlgebra 
using Plots
using Convex
using ECOS


"""
calcAB(dt,alpha)

calculate the A and B matrices for a given timestep Dt and 
mass consumption rate constant alpha
"""
function calcAB(Dt,alpha)
    Ac = zeros(Float64,7,7);
    Ac[1:3,2:4] = Matrix(1.0I, 3, 3);
    
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

function makePSI(A,B,k;m=7,n=4)
    PSI = zeros(Float64,m*(k+2),n*(k+2));
    for i in 2:(k+1)
        for j in 1:i-1
            PSI[(i-1)*m+1:i*m, (j-1)*n+1:j*n] = A^(i-2-(j-1))*B; 
        end
    end
    return PSI;
end

function makeLambda(A,B,k;m=7,n=4)

    LAMBDA = zeros(Float64,m,1);
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

function makeConstraintEq50!(constraints,N,eta,E_u,e_sig;m=7,n=4)
    # UPSILON = makeUpsilon(0)
    # constraints = Constraint[sumsquares(E_u*UPSILON*eta)<=transpose(e_sig)*UPSILON*eta]
    for k in 0:N
        UPSILON = makeUpsilon(N,k)
        push!(constraints,sumsquares(E_u*UPSILON*eta)<=transpose(e_sig)*UPSILON*eta)
    end
    return
end

function solveProblem(N,Dt,alpha)
    A,B = calcAB(Dt,alpha);
    m = size(B,2);
    n = size(B,1);

    # N = Variable();
    eta = Variable(m*(N+1),1);

    e_sig = [0.0; 0; 0; 1];
    E = [Matrix(1.0I,6,6) zeros(Float64,6,1)];
    F = [zeros(Float64,1,6) 1];
    w = Dt;
    omega = w*repeat(e_sig,outer=N+1)#probably need to turn this into a 1d Vector
    E_u = [Matrix(1.0I,3,3) zeros(Float64,3,1)]

    constraints = []

    makeConstraintEq50!(constraints,N,eta,E_u,e_sig)

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
    N_opt = evaluate(N)

    return eta_opt,N_opt
end


Dt = 0.5;
N = 10
alpha = 0.25;
g = [-3.7114; 0; 0];
m_dry = 1505; #kg
m_wet = 1905; #kg
I_sp = 225; #s
T_bar = 3.1; #kN
T_1 = 0.3*T_bar;
T_2 = 0.8*T_bar;
n = 6;
phi = 27*3.1415/180;


solveProblem(N,Dt,alpha)


