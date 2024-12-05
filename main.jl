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

function makeUpsilon(k;m=7,n=4)
    UPSILON = zeros(Float64,m,m*k);
    UPSILON[:,(k-1)*m+1:k*m] = Matrix(0.0I, m, m);
    # UPSILON[:,(k-1)*m+1:k*m] = I;
end

k = 3;
Dt = 0.5;
alpha = 0.25;
# A,B = calcAB(Dt,alpha);
# A = [1 2 3; 1 2 4; 2 2 2]     # 3×3 Matrix{Float64}
A = [1 2; 3 4]
B = [5;6];
UPSILON = makeUpsilon(k;m=2,n=1);
PHI = makePHI(A,k);
PSI = makePSI(A,B,k;m=2,n=1);
LAMBDA = makeLambda(A,B,k;m=2,n=1);


