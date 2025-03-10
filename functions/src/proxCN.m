function y=proxCN(x)
% Proximal point formulation for CN. See appendix A of my phd thesis for reference.
% prox_CN(x) = 0 for x=< 0
%            = x for x > 0
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Code
y=max(x,0);
end