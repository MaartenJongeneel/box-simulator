function [z, basis] = LCP(M,q)
% Solve the LCP using Lempke algorithm given as
% w >= 0
% z >= 0
% w = Mz + q
% z^T*w = 0
%
% INPUTS:    M       : Matrix in R^(n*n)
%            q       : Vector in R^(n)
%
% OUTPUTS:   w       : vector in R^(n)
%            z       : vector in R^(n)
% References
% [1]  The Linear Complementarity Problem, Richard W. Cottle, Jong-Shi 
%      Pang, and Richard E. Stone, 2009
%% Solve LCP

zer_tol = 1e-5;
piv_tol = 1e-8;
big = Inf;
maxiter = 100; %max iterations

n = size(q, 1);
e = ones(n,1);
z = zeros(2*n,1);

%If all elements are positive a trivial solution exists
if all(q >= 0.)
    basis = [];
    z = 0*q;
    return;
end

[tval,lvindex] = max(-q);
lvindex = lvindex(1);

t = 2*n+1;
basis = [n+1:n+lvindex-1 t n+lvindex+1:2*n];
x_B = q+tval*e;
x_B(lvindex) = tval;
B = -eye(n, n);
B(:,lvindex) = e;

leaving = n+lvindex;

for iter=1:maxiter
    if (leaving == t)
        z = zeros(2*n,1);
        z(basis) = x_B; z = z(1:n);
        basis = sort(basis(find(basis<=n)));
        bb = find(z<piv_tol);
        z(bb) = 0;
        basis = setdiff(basis, bb);
        return;
    elseif (leaving < n+1)
        entering = n+leaving;
        %Be = sparse(leaving,1,-1.0,n,1);
        Be = zeros(n, 1);
        Be(leaving) = -1;
    else
        entering = leaving - n;
        Be = M(:,entering);
    end
    d = B\Be;
    theta = big;
    for j=1:n
        adj = abs(d(j));
        if (d(j) > piv_tol)
        	tj = (x_B(j) + zer_tol)/adj;
        	theta = min(theta,tj);
        end
    end
    if (theta >= big)
        disp('unbounded ray');
        z(basis) = x_B; z = z(1:n);
        basis = sort(basis(find(basis<=n)));
        return;
    end

    best_piv = 0.;
    for j=1:n
        tj = big;
        adj = abs(d(j));
        if (d(j) > piv_tol)
        	tj = (x_B(j))/adj;
        	if (tj <= theta)
          	  if (basis(j) == t)
          	    ratio = tj;
        	    leaving = t;
        	    lvindex = j;
        	    break;
        	  elseif (adj > best_piv)
          	    best_piv = adj;
        	    ratio = tj;
        	    leaving = basis(j);
        	    lvindex = j;
              end
            end
        end
    end
    x_B = x_B - ratio*d;
    x_B(lvindex) = ratio;
    B(:,lvindex) = Be;
    basis(lvindex) = entering;
end

basis = sort(basis(find(basis<=n)));

end
