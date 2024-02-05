function [w,z] = LCP(M,q)
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
%% Solve the LCP using Lemke method. See [1] Section 4.4
maxits = 1e4;  %Maximum number of iterations
pivtol = 1e-8; %Maximum pivot tolerance
if min(q)>=0   %If all elements are positive a trivial solution exists
    w = q;
    z = zeros(size(q));
else
    n = length(q);
    %Create initial tableau
    tab = [eye(n), -M, -ones(n, 1), q];
    %Create the basis
    basis = 1:n;
    %First blocking variable is minimum element of q, namely w_r
    [~,row] = min(q);
    basis(row) = 2*n+1; % Replace that index of the basis
    
    %Perform the pivot <w_r,z_0>
    pivot = tab(row,:)/tab(row,2*n+1);
    tab = tab - tab(:,2*n+1)*pivot;
    tab(row,:) = pivot;
    
    %The driving variable is the compliment of w_r, namely z_r:
    col = row + n;
    %Perform complementary pivoting
    loopcount = 0;
    while max(basis) == 2*n+1 && loopcount < maxits%As long as z_0 ~= 0
        loopcount = loopcount +1;
        % Select the driving variable
        driving = tab(:,col);
        
        % Select the blocking variable
        nonzero = driving <= 0;  %Select nonzero entries
        alpha = tab(:,2*n+2)./driving;
        alpha(nonzero) = Inf; %Set the entries of alpha (for which eMs <=0) to Inf
        [~,row] = min(alpha);
        
        % If the driving variable is unblocked, then stop.
        if  sum(nonzero)~=n && abs(driving(row)) > pivtol
            
            %Perform the pivot <blocking,driving>
            pivot = tab(row,:)/tab(row,col);
            tab = tab - tab(:,col)*pivot;
            tab(row,:) = pivot;
            
            %Update the basis
            oldVar = basis(row);
            basis(row) = col;
            
            % The driving variable is the compliment of the preceding
            % blocking variable:
            if oldVar > n
                col = oldVar - n;
            else
                col = oldVar + n;
            end
        else
            break % Break out of the while loop
        end
    end
    % Return the solution
    vars = zeros(2*n+1,1);
    vars(basis) = tab(:,2*n+2).';
    w = vars(1:n,1);
    z = vars(n+1:2*n,1);
end
end

