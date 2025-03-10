function R = Rz(th)
%Rotation matrix around z with th degrees;
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Code
R = [cos(deg2rad(th)) -sin(deg2rad(th)) 0; sin(deg2rad(th)) cos(deg2rad(th)) 0; 0 0 1];
end