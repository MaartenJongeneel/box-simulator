function [AH_B, BV_AB, PN, PT] = BoxSimulatorLCP(x,c,box,surface)
%% Box-simulator-FixedPoint:
%This script uses an LCP formulation
%for solving the nonlinear algebraic equations for contact impact and
%friction. The box has a body fixed frame B at its COM. Furthermore, the
%contact surface has its position and orientation defined by frame C, and
%we express all positions and velocities in terms of frame A, located at
%the camera coordinate frame of a camera, also plotted in the figure.
%
% INPUTS:    x.releasePosition   : 3x1 double, position of the box at release           [m]
%            x.releaseOrientation: 3x3 double, orientation of the box at release        [deg]
%            x.releaseLinVel     : 3x1 double, linear velocity of the box at release    [m/s]
%            x.releaseAngVel     : 3x1 double, angular velocity of the box at release   [rad/s]
%            c.eN                : 1x1 double, normal coefficient of restitution        [-]
%            c.eT                : 1x1 double, tangential coefficient of restitution    [-]
%            c.mu                : 1x1 double, coefficient of friction                  [-]
%            c.dt                : 1x1 double, Timestep at which the simulator runs     [1/s]
%            c.endtime           : 1x1 double, Simulation time you want to run          [s]
%            c.dimd              : 1x1 double, Discretization of the friction cone      [-]    
%            box                 : struct, with fields of box properties as
%                                  box.B_M_B   : 6x6 double intertia tensor of the box  []  
%                                  box.mass    : 1x1 double mass of the box             [kg]
%                                  box.vertices: 3xN double position of the contact     [m]
%                                                points w.r.t the body-fixed frame
%            surface             : Nx1 cell array, with each element a struct containing
%                                  dim      : 1x2 double, dimensions of the surface
%                                  speed    : 3x1 double, speed of the surface
%                                  transform: 4x4 transformation matrix, epxression the 
%                                             position of the surface w.r.t. the world frame
%
% OUTPUTS:   AH_B                : 4x4xN double, transformation matrices expressing the 
%                                  pose of the box w.r.t. world frame over time
%            BV_AB               : 6x1xN double, left trivialized velocity of the box over time
%            FN                  : Normal force acting on the box over time
%            FT                  : Tangential force acting on the box over time
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.

%% Constants and settings
g     = 9.81;                                 %Gravitational acceleration              [m/s^2]
PNfull = zeros(length(box.vertices),1);       %Initial guess for momenta PN            [N*s]
PTfull(1:length(box.vertices),1) = {[0;0]};   %initial guess for momenta PT            [N*s]
N = c.endtime/c.dt;                           %Run to this frame                       [-]

%% Retrieve info of box struct
M = box.B_M_B;
                           
%% Wrench acting on body B: expressed in B with orientation of A:
BA_fo  = [0; 0; -box.mass*g;];
BA_Tau = [0; 0; 0];
BA_f   = [BA_fo; BA_Tau];

%% Initial pose and velocity
%Initial pose of the box: Transformation matrix B w.r.t. frame A
AR_B    = x.releaseOrientation;%Initial orientation
Ao_B    = x.releasePosition;   %Initial position
AH_B(:,:,1) = [AR_B, Ao_B; zeros(1,3), 1]; %Homogeneous transformation matrix

%Initial left trivialized velocity of frame B w.r.t. frame A
Bv_AB = x.releaseLinVel;                      %Linear velocity at t0 of B wrt frame A  [m/s]
Bomg_AB = x.releaseAngVel;                    %Angular velocity at t0 of B wrt frame A [m/s]
BV_AB(:,1) = [Bv_AB; Bomg_AB];

%Discretizing the friction cone in dimd parts
dimd = c.dimd;
D = [];
angles = linspace(0,360,dimd+1);
for ss = 1:length(angles)-1
    theta = angles(ss);
    D = [D; [cos(deg2rad(theta)) sin(deg2rad(theta)) 0]];
end

%Construct the matrix E (diagonal with columns of ones, each of length dimd)
E = zeros(dimd*length(box.vertices),length(box.vertices));
ci = 1; %column index
ri = 1; %row index
for cp = 1:length(box.vertices)
    E(ri:ri+dimd-1,ci) = 1;
    ci = ci +1;
    ri = ri+dimd;
end

%% Dynamics
for t = 1:N %For each time step
    %Rewrite the body velocity as a 4x4 matrix in se(3)
    BV_ABs = hat(BV_AB(:,t));
    
    %Kinematics: Compute the configuration at the mid-time (tM)
    AH_Bm = AH_B(:,:,t)*expm(0.5*c.dt*BV_ABs);
    AR_Bm = AH_Bm(1:3,1:3);
    Ao_Bm = AH_Bm(1:3,4);
    AR_B  = AH_B(1:3,1:3,t);
    
    %Compute the wrench at tM
    B_fM = ([AR_Bm zeros(3); zeros(3) AR_Bm])'*BA_f;
    
    %And compute the gap-functions at tM column of normal contact distances
    gN=[];
    for jj=1:length(surface)
        %Tangential distance w.r.t surface frame
        gT = (surface{jj}.transform(1:3,1:2)'*(Ao_Bm + AR_Bm*box.vertices-surface{jj}.transform(1:3,4)));
        %Vertices outside the contact plane
        vert_o = find(sum((abs(gT)<(0.5*surface{jj}.dim'))==0));
        %Normal distance of the vertices w.r.t. the contact plane
        gN(:,jj) = (surface{jj}.transform(1:3,3)'*(Ao_Bm + AR_Bm*box.vertices-surface{jj}.transform(1:3,4)));

        %If the object is outside the plane area or underneath it, we do
        %not take it into account in solving the contact problem
        if (~isempty(vert_o)) 
            gN(vert_o,jj)=0;
        elseif (all(gN(:,jj)<0))
            gN(:,jj)=0;
        end
    end
       
    %Obtain the linear and angular velocity at tA
    vA = BV_AB(:,t);
    
    %If a gap function has been closed, contact has been made
    [IN,surf] = find(gN<0);
    if  IN > 0
        if length(surf)>1
            surf = unique(surf);
        end
        nn = length(IN);

        %Construct the diagonal friction matrix
        Mu = diag(repmat(c.mu,1,nn));

        indx = gN<0;
        WNAt =[]; WTAt=[]; WNMt=[]; WTMt=[];  WTM=[];
        for jj=1:length(surf)
            WTA{jj} = [];
            WTM{jj} = [];

            %Compute the matrix containing the normal and tangential force directions at tA and tM
            icp = indx(:,surf(jj)); %Indices of the box in contact with surface(jj)
            cp = find(icp>0); %Find the contact point that is in contact with that surface

            %So here we can compute them for all contact points at the time
            %for each surface
            [WNA{jj},~] = CompW(AR_B,surface{surf(jj)}.transform(1:3,1:3),box.vertices(:,icp));
            [WNM,~] = CompW(AR_Bm,surface{surf(jj)}.transform(1:3,1:3),box.vertices(:,icp));

            %Here, for each surface, we loop through the contact points
            for ic =1:length(cp)
               WTA{jj} = [WTA{jj} [D*[AR_B   -AR_B*hat(box.vertices(:,cp(ic)))]]'];
               WTM{jj} = [WTM{jj} [D*[AR_Bm -AR_Bm*hat(box.vertices(:,cp(ic)))]]'];               
            end
            
            WNAt = [WNAt WNA{jj}];
            WTAt = [WTAt WTA{jj}];
            WNMt = [WNMt WNM];
            WTMt = [WTMt WTM{jj}];
        end

        term1 = [hat(vA(4:6)), zeros(3); hat(vA(1:3)), hat(vA(4:6))]*M*vA*c.dt - B_fM*c.dt; 

        %Define the normal and tangential velocities at tA and tM. We
        %substranct from the relative velocity the velocity of the
        %contact surface.
        gammaNA=[]; gammaNE=[]; gammaTA=[]; gammaTE=[]; conv = [];
        for jj=1:length(surf)      
            % gammaNA = [gammaNA; (WNA{jj}'*vA-repmat(surface{surf(jj)}.speed(3,1),sum(indx(:,surf(jj))),1))];
            % gammaTA = [gammaTA; WTA{jj}'*vA];
            speed =surface{surf(jj)}.transform(1:3,1:3)*surface{surf(jj)}.speed;
            conv = [conv; repmat(D*speed,sum(indx(:,surf(jj))),1)];            
        end

        gammaNA = WNAt'*vA;
        gammaTA = WTAt'*vA;

        %Setup the matrices for the LCP (See Claude's thesis, page 128)
        C =[0*diag(ones(1,nn*dimd)), zeros(nn*dimd,nn),  E(1:nn*dimd,1:nn);
            zeros(nn,nn*dimd),      0.0*c.dt*diag(ones(1,nn)),       zeros(nn,nn);
            -E(1:nn*dimd,1:nn)',     Mu,                0*diag(ones(1,nn))];

        G =  [WTMt'; WNMt'; zeros(nn,6)];

        %In the vector below, the restitution parameters appear (not sure this is correct)
        v = [-c.eT*gammaTA+conv; -c.eN*gammaNA; zeros(nn,1)];
        u = M*vA-term1;

        A = C + G/M*G';
        b = G/M*u-v;
       
        %Solve the LCP formulation
        [x,~] = LCP(A,b);

        PT = x(1:dimd*nn);
        PN = x(dimd*nn+1:dimd*nn+nn);

        vE = vA + M\(B_fM*c.dt - [hat(vA(4:6)), zeros(3); hat(vA(1:3)), hat(vA(4:6))]*M*vA*c.dt + WNMt*PN + WTMt*PT);        
        BV_AB(:,t+1) = vE;        
    else
        %Update the velocity to the next time step
        vE = M\(B_fM*c.dt - [hat(vA(4:6)), zeros(3); hat(vA(1:3)), hat(vA(4:6))]*M*vA*c.dt) + vA;
        BV_AB(:,t+1) = vE;
    end
    
    %Complete the time step
    AH_B(:,:,t+1)  = AH_Bm*expm(0.5*c.dt*hat(BV_AB(:,t+1)));
end