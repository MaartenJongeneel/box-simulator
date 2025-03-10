function R = Rx(th)
%Rotation matrix around x with th degrees;
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Code
R = [1 0 0; 0 cos(deg2rad(th)) -sin(deg2rad(th)); 0 sin(deg2rad(th)) cos(deg2rad(th))];
end