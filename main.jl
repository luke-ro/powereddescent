using OrdinaryDiffEq
using LinearAlgebra 
using Plots
using Convex
using ECOS
using Clarabel
using LaTeXStrings
using MathOptInterface


const e_sig = [0.0; 0; 0; 1];
const E = [Matrix(1.0I,6,6) zeros(Float64,6,1)];
const F = [zeros(Float64,1,6) 1];
const E_u = [Matrix(1.0I,3,3) zeros(Float64,3,1)]

mutable struct ProblemParameters
    g
    alpha
    m_dry
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
    Bc[7,4] = -alpha;

    # function intB!(dB, B, p, t)
    #     temp = exp(Ac*(Dt-t))*Bc;
    #     dB[:] = reshape(temp,:,1);
    #     # println(dB)
 
    # end

    # A = exp(Ac*Dt);
    # prob = ODEProblem(intB!,zeros(Float64,28,1),(0,Dt));
    # B = solve(prob, Tsit5()); 
    # # print(size(B))
    # B = reshape(B[:,1,end],7,4);

    A_hat = [Ac Bc; 
            zeros(4,11)]
    
    temp = exp(A_hat*Dt)
    A = temp[1:7,1:7]
    B = temp[1:7,8:11]

    return A,B;
end

function makePHI(A,k)
    return A^k;
end

function makePSI(A,B,k,N;m=7,n=4)
    PSI = zeros(Float64,m,n*(N+1)) # number of columns must match number of rows of eta
    for i in 1:k
        PSI[:,(i-1)*n+1:i*n] = A^(k-1-(i-1))*B; 
    end
    # print("size(psi): "); println(size(PSI))
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
    omega = w*repeat(e_sig, outer=N+1)
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

function makeRho(num,ps::ProblemParameters)
    if num==1
        return rho = ps.n*ps.T_1*cos(ps.phi) #eq3
    else
        return rho = ps.n*ps.T_2*cos(ps.phi) #eq3
    end
end

function makeMu(num,k,Dt,ps::ProblemParameters)
    t = Dt*k;
    if num==1
        rho = makeRho(num,ps)
    elseif num==2
        rho = makeRho(num,ps)
    else
        throw("Unrecognized number for Mu in MakeMu")
    end
    z0 = makez0(t,ps) #eq35
    mu = rho*exp(-z0) #eq33
    return mu
end

function simulateTrajectory(simulation_Dt,guidance_update_Dt,ps::ProblemParameters)
    F,G = calcAB(simulation_Dt, ps.alpha)

    guidance_Dk = Int(round(guidance_update_Dt/simulation_Dt));

    x_hist = []
    t_hist = []
    control_hist = []

    x = ps.y0
    eta = []
    t = 0
    k = 0
    completed = true
    while x[1] > 0.0
        #guidnace loop
        if k % guidance_Dk == 0
            ps.m_wet = exp(x[7])
            eta,N,feasible = runGuidance(x,simulation_Dt,ps)
            if !feasible
                completed = false
                break
            end
        end

        i = (k % guidance_Dk)+1
        u = eta[(i-1)*4+1: i*4]

        x = F*x + G*(u .+ vcat(ps.g,[0]))

        push!(x_hist,x)
        push!(t_hist,t)
        push!(control_hist, u)

        println(x)
        t += simulation_Dt
        k += 1
    end
    return transpose(reduce(hcat,t_hist)), reduce(hcat,x_hist), reduce(hcat,control_hist), completed
end

function runGuidance(x,Dt,ps::ProblemParameters)
    # run low fidelity optimization
    eta,N,A,B,solved = optimizeProblem(x,Dt,ps)
    # print("Solved: ")
    # println(solved)

    # run higher fidelity optimization
    # eta_prime,cost,A,B,solved = optimizeProblem(x,N,Dt,ps)

    #return control
    return eta,N,solved
end

function addConstraints!(constraints,x0,A,B,N,Dt,eta,omega,params::ProblemParameters;m=7,n=4)
    n_s = size(a,1) #number of constraints for ||S*x-v||+c^T*x + a <= 0

    #constraint in eq50
    for k in 0:N
        UPSILON = makeUpsilon(N,k)
        temp = norm(E_u*UPSILON*eta) <= transpose(e_sig)*UPSILON*eta
        push!(constraints,temp)
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
        XI_k = PHI_k*x0 + LAMBDA_k*vcat(params.g, [0])

        # println("here:")
        # print("size(XI_k+PHI_k*eta): "); println(size(XI_k+PHI_k*eta))
        # println(N)
        # print(size(XI_k+(PSI_k*eta)))
        
        # eq51 part 1
        temp = mu1_k*(1-(F*(XI_k+PSI_k*eta)-z0_k)+((F*(XI_k+PSI_k*eta)-z0_k)*(F*(XI_k+PSI_k*eta)-z0_k))/2) <= transpose(e_sig)*UPSILON_k*eta
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
        for j in 1:n_s
            temp = norm(params.S[j]*E*(XI_k+PSI_k*eta)-params.v[j]) + transpose(params.c[j])*E*(XI_k+PSI_k*eta) + params.a[j] <= 0
            push!(constraints,temp)
        end

        if k == N
            temp = [Matrix(1.0I,6,6) zeros(6,1)]*(XI_k+PSI_k*eta) == 0*ones(6,1)
            push!(constraints,temp)

            temp = [0 1 0 0;
                    0 0 1 0]*UPSILON_k*eta == [0;0]
            push!(constraints,temp)

            temp = [1 0 0 0]*UPSILON_k*eta >= 0
            push!(constraints,temp)
        end
    end

    return 
end

function optimizeProblem(x0,Dt,ps::ProblemParameters; N_min=[], N_max=[])
    rho1 = makeRho(1,ps)
    rho2 = makeRho(2,ps)
    t_min = (ps.m_wet-ps.m_dry)*norm(x0[4:6])/rho2
    t_max = (ps.m_wet-ps.m_dry)/ps.alpha/rho1
    if N_min == []
        N_min = Int(floor(t_min/Dt)+1)
    end

    if N_max == []
        N_max = Int(floor(t_max/Dt)+1)
    end

    cost_min = typemax(Float64)
    eta_opt = []
    A = []
    B = []
    N_opt = N_min
    for N in N_min:N_max
        # print("Checking N=$(N): ")
        eta, cost, A, B, feasible = solveSubproblem(x0,N,Dt,ps)
        if feasible && cost < cost_min 
            # println("More optimal. Cost of $(cost)")
            cost_min = cost
            eta_opt = eta
            N_opt = N
        else
            # println()
        end
    end 

    solved = true
    if eta_opt == []
        solved = false
    end
    return eta_opt,N_opt,A,B,solved
end

function solveSubproblem(x0,N,Dt,params)
    A,B = calcAB(Dt,params.alpha);
    m = size(B,2);
    n = size(B,1);

    # N = Variable();
    eta = Variable(m*(N+1),1);

    w = Dt;
    omega = makeomega(w,N)#probably need to turn this into a 1d Vector

    constraints = []
    addConstraints!(constraints,x0,A,B,N,Dt,eta,omega,params)

    cost = transpose(omega)*eta;
    problem = minimize(cost,constraints)
    solver = ECOS.Optimizer

    Convex.solve!(problem, solver, silent=true)

    eta_opt = evaluate(eta)
    # N_opt = evaluate(N)

    feasible = true
    if problem.status == MathOptInterface.INFEASIBLE
        feasible = false    
    end 
    
    cost = transpose(omega)*eta_opt
    return eta_opt, cost, A, B, feasible
end

function eom(t,x,u,A,B,ps::ProblemParameters)
    dx = A*x + B*(vcat(ps.g,[0]) + u)
    return dx
end

function calcTrajectory(Dt,N,A,B,eta,ps::ProblemParameters)
    x = []
    u = []
    m = []
    throttle = []
    theta_hist = []
    for k = 0:N
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
        mass = exp(state[7])
        theta = acos(dot(control[1:3]/norm(control[1:3]),[1;0;0]))

        push!(x,state[1:6])
        push!(u,control[1:3])
        push!(m,mass)
        push!(throttle,control[4]*mass/ps.T_bar/cos(ps.phi)/6)
        push!(theta_hist,theta)
    end
    return reduce(hcat,x),reduce(hcat,u),reduce(hcat,m),reduce(hcat,throttle),reduce(hcat,theta_hist)
end

function makePlots(t,pos,vel,acc,force,throttle,theta,params::ProblemParameters)
    data = []
    fs = 8

    # println(t)
    # println(throttle)
    # p1 = plot([0; 1],[0; 1])
    p1 = plot(t,transpose(pos)/1000,tickfontsize=fs,linewidth=2,legend=false)#,labels=["x" "y" "z"])
    ylabel!("Position [km]",labelfontsize=fs)
    p2 = plot(t,transpose(vel),tickfontsize=fs,linewidth=2,legend=false)
    ylabel!("Velocity [m/s]",labelfontsize=fs)
    p3 = plot(t,transpose(acc),tickfontsize=fs,linewidth=2,legend=false)
    ylabel!("Acceleration [m/s/s]",labelfontsize=fs)
    p4 = plot(t,transpose(force),tickfontsize=fs,linewidth=2,legend=false)
    ylabel!("Control Force",labelfontsize=fs)
    p5 = plot(t,transpose(throttle),tickfontsize=fs,linewidth=2,legend=false)
    ylabel!("Throttle %",labelfontsize=fs)
    p6 = plot(t,transpose(theta).*180/3.1415,tickfontsize=fs,linewidth=2,legend=false)
    ylabel!(L"\theta ",labelfontsize=fs)

    # Add each plot to the subplot grid
    p = plot(p1,p2,p3,p4,p5,p6,layout=(3,2))
  
    # Display the plot
    display(p)
end

function calcTrajParams(t_hist,x_hist,u_hist,ps::ProblemParameters)
    pos_hist = []
    vel_hist = []
    acc_hist = []
    force_hist = []
    throttle_hist = []
    theta_hist = []

    for i in 1:size(t_hist,1)
        t = t_hist[i]
        x = x_hist[:,i]
        u = u_hist[:,i]
        m = exp(x[7])
        throttle = u[4]*m/ps.T_bar/cos(ps.phi)/ps.n
        theta = acos(dot(u[1:3]/norm(u[1:3]),[1;0;0]))
        
        
        push!(pos_hist,x[1:3])
        # println(x)
        push!(vel_hist,x[4:6])
        push!(acc_hist,u[1:3])
        push!(force_hist,u[1:3].*m)
        push!(throttle_hist,throttle)
        push!(theta_hist,theta)

        # print("Throttle hist:")
        # println(reduce(hcat,throttle_hist))
    end
    return reduce(hcat,pos_hist),reduce(hcat,vel_hist),reduce(hcat,acc_hist),reduce(hcat,force_hist),reduce(hcat,throttle_hist),reduce(hcat,theta_hist)
end



g = [-3.7114; 0; 0];
m_dry = 1505; #kg
m_wet = 1905; #kg
I_sp = 225; #s
T_bar = 3.1*1000; #kN -> N
T_1 = 0.3*T_bar;
T_2 = 0.8*T_bar;
n = 6;
phi = 27.0*3.1415/180;
r0 = [1.5,0,2]*1000 #km -> m
dr0 = [-75,0,0] #m/s
# r0 = [2,0,0]
# dr0 = [-50,0,0]
y0 = vcat(r0,dr0,log(m_wet))
alpha = 1/(I_sp*9.807*cos(phi));
theta_tilde = 86.0*3.1415/180;
S = [[0 1 0 0 0 0;
     0 0 1 0 0 0]]
v = [[0;0]]
c = [[-tan(theta_tilde); 0; 0; 0; 0; 0]]
a = [[0]]
params = ProblemParameters(g,alpha,m_dry,m_wet,I_sp,T_bar,T_1,T_2,n,phi,y0,S,v,c,a)
# print(params)

Dt = 5;
# eta_opt, N, A, B, solved = optimizeProblem(Dt,params)
# [t,x] = simulateProblem(A,B,)
# print(reshape(eta_opt,4,:))

t_hist,x_hist,u_hist,completed = simulateTrajectory(4,4,params)
pos_hist,vel_hist,acc_hist,force_hist,throttle_hist,theta_hist = calcTrajParams(t_hist,x_hist,u_hist,params)
makePlots(t_hist,pos_hist,vel_hist,acc_hist,force_hist,throttle_hist,theta_hist,params)
# print(t_hist)

optimizeProblem(y0,Dt,params)

# eta = transpose(reshape(eta_opt,4,:))
# x,u,mass_hist,throttle_hist,theta_hist = calcTrajectory(Dt,N,A,B,eta_opt,params::ProblemParameters)
# println(size(x))
# t=0:Dt:(N)*Dt
# makePlots(t,eta_opt,x,u,mass_hist,throttle_hist,theta_hist,params)

