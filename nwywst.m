% Autocovariance consistent variance correction (by G. Villa Cox & K. Fridirikson)
%   X:  Matrix of regressors
%   u:  OLS residual errors
%   L:  Maximum # of lags
function [vcv]=nwywst(X,u,L)
    
    [T,k]=size(X);
        
    vcv=zeros(k);
	
	for j=0:L
		w=1-(j/(L+1));
		for t=j+1:T
			if (j==0)
				% Calculates the S_0 portion
				upd=u(t)^2*(X(t,:)'*X(t,:));
				vcv=vcv+upd;
			else
				% Calculates the off-diagonal terms
				upd=w*u(t)*u(t-j)*(X(t-j,:)'*X(t,:) + X(t,:)'*X(t-j,:));
				vcv=vcv+upd;
			end
		end
	end
	
	% Newey-White Variance Covariance Matrix
	vcv=1/(T-k) * vcv;
    xx1=(X'*X)\eye(k);
    vcv=T*xx1*vcv*xx1;
end

