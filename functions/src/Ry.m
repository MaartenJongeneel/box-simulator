function R = Ry(th)
%Rotation matrix around y with th degrees;
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Code
R = [cos(deg2rad(th)) 0 sin(deg2rad(th)); 0 1 0;  -sin(deg2rad(th)) 0 cos(deg2rad(th))];
end