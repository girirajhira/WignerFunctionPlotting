
using Base.Threads


function meshgrid(x,y)
    return ones(size(y))*x', y*ones(size(x))' 
end


function GenLaguerre(n::Int, alpha, x)
    # Generalised Laguerre polynomial 
    if n>1 # Gets called the most, so the default case.
        lagpast = 1;
        lagnow = 1 + alpha .- x;
        lagtemp =  0*lagnow
        for k = 1:(n-1)
            lagtemp = @. ((2*k+1+alpha-x)*lagnow - (k+alpha)*lagpast)/(k+1);
            lagpast, lagnow = lagnow, lagtemp
        end
        return lagnow;
    elseif  n==0 
        return 1 .+ 0*x
    else  # n == 1
        return 1 + alpha .- x
    end
end


function fastfact(n::Int)
    # returns factorial of n,
    # If n>11 then the Stirling Approximation is used
    e = exp(1)
    if n < 11
        return factorial(n)
    else
        return @fastmath sqrt(2*3.14159*n)*(n/2.7182818)^n * (1 + 1/(12*n) + 1/(288*n^2))
        # return @fastmath sqrt(2*pi*n)*(n/e)^n * (1 + 1/(12*n))
    end
end



function Wnm(m::Int, n::Int, x, p)
    # Weyl Symbols of |n><m|, only real if n=m
    hbar = 1
    alpha = (x.^2 .+ p.^2)
    if n == m
        pf = (-1)^m/(pi*hbar)
        t2 = exp.(-alpha/hbar)
        t3 = GenLaguerre(m, n-m, alpha/(hbar/2))
        return @. pf*t2*t3
    else
        pf = sqrt(fastfact(m)/fastfact(n))*(-1)^m/(pi*hbar)
        t1 = @. exp(1im*(m-n)*atan.(p, x))*(alpha/(hbar/2))^((n-m)/2)
        t2 = exp.(-alpha/hbar)
        t3 = GenLaguerre(m, n-m, alpha/(hbar/2))
        return @. pf*t1*t2*t3
    end
    
end


function WignerFunction(rho, x, p)
    # Given a density matrix rho, this function returns 
    # given Wigner function evaluated on the grid x,p
    xmesh, pmesh = meshgrid(x,p)
    nmax, mmax = size(rho)
    wf = zeros(size(xmesh))
    tolerance = 1e-12
    # Diagonal Terms
    wfarr = zeros((nmax,length(x), length(p)))
    @simd for n in 1:nmax
        if abs2(rho[n,n]) > tolerance
             wf = @. wf + real(rho[n,n]*Wnm(n-1, n-1, xmesh, pmesh))
        end
    end
    # Off Diagonal
    @simd for n in 1:(nmax)
        for m in 1:(n-1)
            if abs2(rho[m,n]) > tolerance
                wf = @. wf + 2*real(rho[m, n]*Wnm(m-1, n-1, xmesh, pmesh))
            end 
        end
    end

    # Normalisation
    return real(wf) #/(sum(real(wf))*abs(x[2]-x[1])*abs(p[2]-p[1]))
end



# Optimised
# +++++++++++++++++++++++++++++++++++++++++++++++++++
function GenLaguerreOptimised(n, alpha, x)
    # Generalised Laguerre polynomial 
    if n==0
        return 1 .+ 0*x
    elseif n==1
        return 1 + alpha .- x
    else
        lagpast = 1;
        lagnow = 1 + alpha .- x;
        lagtemp =  fill(0.0, size(x))
        for kk = 2:n
            k = kk-1
            lagtemp = ((2*k+1+alpha.-x).*lagnow .- (k+alpha)*lagpast)/(k+1);
            lagpast = lagnow;
            lagnow = lagtemp;
        end
        return lagnow;
    end
end



function WnmOptimised(lastlag, m, n, x, p)
    hbar = 1
    alpha = (x.^2 .+ p.^2)
    if n == m
        pf = (-1)^m/(pi*hbar)
        t2 = exp.(-alpha/hbar)
        t3 = GenLaguerreOptimised(lastlag, m, n-m, alpha/(hbar/2))
        return pf.*t2.*t3, t3
    else
        pf = sqrt(fastfact(m)/fastfact(n))*(-1)^m/(pi*hbar)
        t1 = exp.(1im*(m-n)*atan.(p, x)).*(alpha/(hbar/2)).^((n-m)/2)
        t2 = exp.(-alpha/hbar)
        t3 = GenLaguerreOptimised(lastlag, m, n-m, alpha/(hbar/2))
        return pf.*t1.*t2.*t3, t3
    end
   
end


function WignerFunctionFast(rho, x, p)
    xmesh, pmesh = meshgrid(x,p)
    nmax, mmax = size(rho)
    wf = zeros(size(xmesh))
    tolerance = 1e-12
    # Useful definaitons
    hbar = 1
    alpha = @. (xmesh^2 + pmesh^2)
    expAlpha2 = exp.(-alpha/hbar)
    theta = atan.(pmesh,xmesh)
    # Diagonal Terms
    lag1 = @. 1 + 0*xmesh
    lag2 = @. 1 - alpha/(hbar/2)
    for n in 0:(nmax-1)
        if n > 1
            lag2, lag1 = @. ((2*n+1 - alpha/(hbar/2))*lag2 - n*lag1)/(n+1), lag2
            if abs2(rho[n+1,n+1]) > tolerance
                wnn = @. (-1)^(n)/(pi*hbar)*expAlpha2*lag2
                wf = @. wf + real(rho[n+1,n+1]*wnn)
            end
        else
            if abs2(rho[n+1,n+1]) > tolerance
                wnn = @. (-1)^(n)/(pi*hbar)*expAlpha2*GenLaguerre(n, 0, alpha/(hbar/2))
                wf = @. wf + real(rho[n+1,n+1]*wnn)
            end
        end
    end
    # Off Diagonal
    for n in 0:(nmax-1)
        for m in 0:(n-1)
            if abs2(rho[m+1,n+1]) > tolerance
                pf = sqrt(fastfact(m)/fastfact(n))*(-1)^m/(pi*hbar)
                t1 = @. exp(1im*(m-n)*theta)*(alpha/(hbar/2))^((n-m)/2)
                wmn = @. pf*t1*expAlpha2*GenLaguerre(m,n-m,alpha/(hbar/2))
                wf = @. wf + 2*real(rho[m+1, n+1]*wmn)
            end 
        end
    end
    return -wf
end




    # function WignerFunctionFast(rho, x, p)
    #     xmesh, pmesh = meshgrid(x,p)
    #     nmax, mmax = size(rho)
    #     wf = zeros(size(xmesh))
    #     tolerance = 1e-12
    #     # Useful definaitons
    #     hbar = 1
    #     alpha = @. (xmesh^2 + pmesh^2)
    #     expAlpha2 = exp.(-alpha/hbar)
    #     theta = atan.(pmesh,xmesh)
    #     # Diagonal Terms
    #     lag1 = @. 1 + 0*xmesh
    #     lag2 = @. 1 - x
    #     for n in 0:(nmax-1)
    #         if n > 1
    #             lag2, lag1 = @. ((2*n+1 - x)*lag2 - n*lag1)/(n+1), lag2
    #             if abs2(rho[n+1,n+1]) > tolerance
    #                 wnn = @. (-1)^(n)/(pi*hbar)*expAlpha2*GenLaguerre(n, 0, alpha/(hbar/2))
    #                 wf = @. wf + real(rho[n+1,n+1]*wnn)
    #             else 
    #                 wnn = @. (-1)^(n)/(pi*hbar)*expAlpha2*GenLaguerre(n, 0, alpha/(hbar/2))
    #                 wf = @. wf + real(rho[n+1,n+1]*wnn)
    #             end
    #         end
    #     end
    #     # Off Diagonal
    #     for n in 0:(nmax-1)
    #         for m in 0:(n-1)
    #             if abs2(rho[m+1,n+1]) > tolerance
    #                 pf = sqrt(fastfact(m)/fastfact(n))*(-1)^m/(pi*hbar)
    #                 t1 = @. exp(1im*(m-n)*theta)*(alpha/(hbar/2))^((n-m)/2)
    #                 wmn = @. pf*t1*expAlpha2*GenLaguerre(m,n-m,alpha/(hbar/2))
    #                 wf = @. wf + 2*real(rho[m+1, n+1]*wmn)
    #             end 
    #         end
    #     end


    # Normalisation
#     return -real(wf) #/(sum(real(wf))*abs(x[2]-x[1])*abs(p[2]-p[1]))
# end




function WignerFunction2(rho, x, p)
    xmesh, pmesh = meshgrid(x,p)
    alpha = xmesh+ 1im*pmesh
    wf = fill(0.0, size(xmesh))
    tolerance = 1e-16
    N, N2 = size(rho)
    A, Ad, psi0 = SingleMode(N)
    parity = exp(1im*pi*Ad*A)
    # Diagonal Terms
    for i in 1:length(x)
        for j in 1:length(p)
            D1 = exp((alpha[i, j])*Ad .- conj(alpha[i, j])*A)
            wf[i, j] = real(tr(D1*parity*D1'*rho))
        end
    end
  

    # Normalisation
    return real(wf) #/(sum(real(wf))*abs(x[2]-x[1])*abs(p[2]-p[1]))
end

function PlotWignerFunction(wf, x, p)
  heatmap(x, p, wf, cmap=:bwr)
end

# x = -2:0.32:2
# p = -2:0.32:2

# WignerFunction([1. 0.;0. 1.], x, p)
# @profile WignerFunction([1. 0.;0. 1.], x, p)
# ProfileView.view()

