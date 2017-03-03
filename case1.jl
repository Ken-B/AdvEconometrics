using Distributions

type ols_results
  b::Array{Float64,2}
  yhat::Array{Float64,2}
  u::Array{Float64,2}
  vcv::Array{Float64,2}
  tstat::Array{Float64,2}
  pval::Array{Float64,2}
end

function gls(y,X; nwc=false,L=0)

  (T,k) = size(X)
  Ik = eye(k,k)
  XtX= X' * X
  XtXInv = XtX\Ik
  b = XtXInv * (X' * y)
  Px = X * XtXInv * X'
  yhat = Px * y
  u = (Ik - Px) * y
  sig2 = dot(u,u) / (T-k)

  if (nwc)
    vcv = nwywst(X,u,L)
  else
    vcv = sig2 * XtXInv
  end

  tstat = b ./ sqrt(diag(vcv))
  pval =

  return ols_results(b,yhat,u,vcv,tstat,pval)

end

# Autocovariance consistent variance correction (by G. Villa Cox & K. Fridirikson)
#   X:  Matrix of regressors
#   u:  OLS residual errors
#   L:  Maximum # of lags
function nwywst(X,u,L)

  (T,k) = size(X)
  vcv = zeros(k,k)

  for j in (0:L)

    w=1-(j/(L+1))
    for t in (j+1:T)
      if (j==0)
        # Calculates the S_0 portion
        upd=u(t)^2*(X(t,:)'*X(t,:))
      else
        # Calculates the off-diagonal terms
				upd=w*u(t)*u(t-j)*(X(t-j,:)'*X(t,:) + X(t,:)'*X(t-j,:))
      end
      vcv=vcv+upd
    end

  end

  # Newey-White Variance Covariance Matrix
	vcv=1/(T-k)*vcv
  xx1=(X'*X)\eye(k)
  vcv=T*xx1*vcv*xx1

  return vcv
end

function gls(y,X,Sigma)

  SigmaInv = Sigma\eye(size(Sigma))
  P=chol(SimgaInv)
  return ols(P*y,P*X)

end
