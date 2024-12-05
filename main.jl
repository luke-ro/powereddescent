using OrdinaryDiffEq
using LinearAlgebra 
using Plots

function calcAB(dt,alpha)
    Ac = zeros(Float64,7,7);
    Ac[1:3,2:4] = Matrix(1.0I, 3, 3);
    
    Bc = zeros(Float64, 7, 4);
    Bc[4:6,1:3] = Matrix(1.0I, 3, 3);
    Bc[7,4] = alpha;

    function intB!(dB, B, p, t)
        temp = exp(Ac*(dt-t))*Bc;
        dB[:] = reshape(temp,:,1);
        # println(dB)
 
    end

    A = exp(Ac*dt);
    prob = ODEProblem(intB!,zeros(Float64,28,1),(0,dt));
    B = solve(prob, Tsit5()); 
    # print(size(B))
    B = reshape(B[:,1,end],7,4);

    return A,B;
end

Dt = 0.5;
alpha = 0.25;
A,B = calcAB(Dt,.alpha);

