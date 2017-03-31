#Load the Distributions package. Use `Pkg.install("Distributions")` to install first time.
using Distributions: TDist, ccdf

type regress_results
    coefs
    yhat
    res
    vcv
    tstat
    pval
end

# Keyword arguments are placed after semicolon.
# Symbols start with colon, e.g. `:symbol`.
function ols(y, X; corr=:none, lags::Int=Int(floor(size(X,1)^(1/4))))

    # β̂ = X \ y is more stable than β̂ = inv(X'*X) * X' \ y 
    # see notes at bottom of case 1 notebook
    β̂ = X \ y
    ŷ = X * β̂
    μ̂ = y - ŷ

    T, K = size(X)
    σ̂² = dot(μ̂, μ̂) / (T - K)

    #use correction for variance covariance
    if corr == :none
        vcv = σ̂² * inv(X'*X)
    elseif corr == :white
        vcv = newey_west(X, μ̂, 0)
    elseif corr == :newey_west
        vcv = newey_west(X, μ̂, lags)
    else
        error("wrong argument for correction keyword")
    end
    
    # T statistics for H₀: βᵢ = 0
    tstat = β̂ ./ sqrt(diag(vcv))

    # absolute value and times two for double sided test
    pval  = 2 * ccdf(TDist(T-K), abs(tstat))

    regress_results(β̂, ŷ, μ̂, vcv, tstat, pval)
end


function newey_west(X, μ̂, lags::Integer)

    XtXInv = inv(X'*X)
    T, K = size(X)
    
    if lags==0 # White estimator
        return XtXInv * X' * diagm(μ̂.^2) * X * XtXInv
    end

    vcv = zeros(K, K)
    for t = 1:T
        vcv += μ̂[t]^2 * (X[t,:] * X[t,:]')
    end
    for lag in 1:lags
        w = 1 - lag / (lags + 1)
        for t in (lag + 1):T
            # Calculates the off-diagonal terms
            vcv += w * μ̂[t] * μ̂[t-lag] * (X[t-lag,:]*X[t,:]' + X[t,:]*X[t-lag,:]')            
        end
    end
    vcv = XtXInv * vcv * XtXInv
end

function gls(y, X, Ω)

  P = chol(inv(Ω))
  return ols(P*y, P*X)

end

function gmm(y, X, Z; corr=:none, lags=nothing)

  T, Kx = size(X)
  T, Kz = size(Z)

  if corr==:none

    # Generalized 1-step IV estimator
    W = inv(Z'*Z)

  elseif corr==:white | corr==:newey_west

    if corr==:white
      lags=0
    end

    gmm1_res = gmm(y,X,Z;corr=:none)
    μ̂ = gmm1_res.res

    W = zeros(Kz, Kz)
    for lag in 0:lags
        w = 1 - lag / (lags + 1)
        for t in (lag + 1):T
            # Calculates the off-diagonal terms
            update = w * μ̂[t] * μ̂[t-lag] * (Z[t-lag,:]*Z[t,:]' + Z[t,:]*Z[t-lag,:]')
            W = W + update
        end
    end

  else
      error("wrong argument for correction keyword")
  end

  ZtX = Z'*X
  XtZ = X'*Z

  XtZ_W_ZtXInv = inv(XtZ*W*ZtX)
  β̂  = XtZ_W_ZtXInv*(XtZ*W*Z'*y)
  ŷ = X * β̂
  μ̂ = y - ŷ
  σ̂² = dot(μ̂, μ̂) / (T - Kz)
  vcv = σ̂² * XtZ_W_ZtXInv

  # T statistics for H₀: β₀ = 0
  tstat = β̂ ./ sqrt(diag(vcv))

  # absolute value and times two for double sided test
  pval  = 2 * ccdf(TDist(T-K), abs(tstat))

  return regress_results(β̂, ŷ, μ̂, vcv, tstat, pval)

end
