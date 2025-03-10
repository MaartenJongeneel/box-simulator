function [WN,WT] = CompW(AR_B,AR_C,vertices)
% This script compute the matrix containing the normal and tangential force directions.
%
% INPUTS:    AR_B        : 3x3 double, orientation of the box w.r.t world frame [deg]
%            AR_C        : 3x3 double, orientation of the contact surface w.r.t world frame [deg]
%            vertices    : 3xN double, positions of vertices of the box w.r.t box frame [m]
%
% OUTPUTS:   WN          : 1xN double, matrix with force directions in normal direction [-]
%            WT          : 2xN double, matrix with force directions in tangential direction [-]
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Code
tel = 1;
for ii = 1:length(vertices(1,:))
    w = (AR_C'*[AR_B -AR_B*hat(vertices(:,ii))])';
    WN(:,ii) = w(:,3);
    WT(:,tel:tel+1) = w(:,1:2);
    tel = tel+2;
end
end