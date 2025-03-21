% function [AH_B, BV_AB, FN, FT] = BoxSimulator(x,c,obj,surface)
%% Box-simulator-FixedPoint:
%This script uses an augmented Lagrangian approach (fixed-point iteration)
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
%            c.a                 : 1x1 double, Prox point auxilary parameter            [-]
%            c.tol               : 1x1 double, Error tol for fixed-point                [-]
%            obj                 : cell array where each object is a struct, with fields of box properties as
%                                  obj{ii}.B_M_B    : 6x6 double intertia tensor of the box  []  
%                                  obj{ii}.mass     : 1x1 double mass of the box             [kg]
%                                  obj{ii}.vertices : 3xN double position of the contact     [m]
%                                                    points w.r.t the body-fixed frame
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
% Copyright (c) 2025, Maarten Jongeneel
% All rights reserved.
%% Constants and settings
close all;clearvars;clc;
g     = 9.81;                                 %Gravitational acceleration              [m/s^2]
N = 1500;%1500;
c.dt = 1/1000;
c.eN = 0.2; %0.4;
c.eT = 0;
c.mu = 0.2;
c.a = 0.0001;
c.tol = 1e-7;
c.dimd = 8; %discretization of the friction cone
% Bv_AB = x.releaseLinVel;                      %Linear velocity at t0 of B wrt frame A  [m/s]
% Bomg_AB = x.releaseAngVel;                    %Angular velocity at t0 of B wrt frame A [m/s]
PNfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
LNfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
lnfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
PTfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
LTfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
ltfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
% N = round(c.endtime/c.dt);                           %Run to this frame                       [-]

%% Load the scene
% Multi3Box;
% MultiWall;
SCS700ChuteWidth;


figure('pos',[200,200,800,600]); 
hold on;
for ii = 1:length(obj)
    plotBox(obj{ii}.initAH_B,obj{ii},[0.7 0.7 0.7]);
    hold on
    % for jj = 1:6
    %     AH_C = obj{ii}.initAH_B*obj{ii}.surface{jj}.transform;
    %     tip = [AH_C(1:3,4)+0.1*AH_C(1:3,1) AH_C(1:3,4)+0.1*AH_C(1:3,2) AH_C(1:3,4)+0.1*AH_C(1:3,3)];
    %     plot3([AH_C(1,4) tip(1,1)],[AH_C(2,4) tip(2,1)],[AH_C(3,4) tip(3,1)],'r'); hold on
    %     plot3([AH_C(1,4) tip(1,2)],[AH_C(2,4) tip(2,2)],[AH_C(3,4) tip(3,2)],'g');
    %     plot3([AH_C(1,4) tip(1,3)],[AH_C(2,4) tip(2,3)],[AH_C(3,4) tip(3,3)],'b');
    % end
end
axis equal
% view(34,20);
view(-40,26);
ax = gca;
ax.Clipping = "off";
axis off;
% zoom(2)

%% Set up the equations of motion
% Gravity wrench acting on Body B
BA_f = [];
BV_AB =[];
B_M_B =[];
ii = 0;
for cnt = 1:length(obj) %For each object
    ii = ii+1;
    % Initial pose and velocity
    AH_B(:,:,ii,1) = obj{ii}.initAH_B;
    BA_fo  = [0; 0; -obj{ii}.mass*g;];
    BA_Tau = [0; 0; 0];
    BA_f   = [BA_f; BA_fo; BA_Tau];
    BV_AB = [BV_AB; obj{ii}.initBV_AB];
    % Mass matrix
    B_M_B = [B_M_B, zeros(length(B_M_B),6); zeros(6,length(B_M_B)), obj{ii}.B_M_B];
end

Nobj = length(obj); %Number of bodies in the simulation

%Setup index to track which surface of which body 
cnts=0; %count surface
for ii=1:Nobj
    for kk = 1:Nobj
        nsurf = length(obj{kk}.surface); %number of surfaces of object
        for jj = 1:nsurf
            cnts = cnts+1;
            cidx(1,cnts) = ii; %This is the body we consider
            cidx(2,cnts) = kk; %Which is in contact with this body
            cidx(3,cnts) = jj; %And the contact is on this surface
        end
    end
end

nn = 100; % Max number of contacts allowed
%Constraint stabilization matrices (See Claude's thesis, page 128 and 214)
Gamma = 0.0*diag(ones(1,nn*c.dimd));
Sigma = 0.0*diag(ones(1,nn));
Delta = 0.0*diag(ones(1,nn));
Upsilon = 0.0*ones(nn,1); %zeros(nn,1); %1/(1+4/c.dt)*(ones(nn,1));

%Construct the matrix E (diagonal with columns of ones, each of length dimd)
E = zeros(c.dimd*nn,nn);
ci = 1; %column index
ri = 1; %row index
for cj = 1:nn
    E(ri:ri+c.dimd-1,ci) = 1;
    ci = ci +1;
    ri = ri+c.dimd;
end

%% Dynamics
for ii=1:N
    % ---------------- STEP 1: Compute system matrices at mid point ---------------- %

    %Obtain the linear and angular velocity at tA
    tic
    vA = BV_AB(:,ii);
    vACross = [];
    cnt = 0;
    for kk = 1:Nobj
        cnt = cnt+1;
        if obj{kk}.dynamics == false
            vtemp = zeros(6,1);
        else
            vtemp = BV_AB((6*(cnt-1)+1:6*(cnt-1)+6),ii); %velocity of body kk
        end
        %Kinematics: Compute the configuration at the mid-time (tM)
        AH_Bm(:,:,kk) = AH_B(:,:,kk,ii)*expm(0.5*c.dt*hat(vtemp));
        AR_Bm(:,:,kk) = AH_Bm(1:3,1:3,kk);

        %Compute the wrench at tM
        if obj{kk}.dynamics
            B_fM((6*(cnt-1)+1:6*(cnt-1)+6),1) = ([AR_Bm(:,:,kk) zeros(3); zeros(3) AR_Bm(:,:,kk)])'*BA_f((6*(cnt-1)+1:6*(cnt-1)+6),1);
        else
            B_fM((6*(cnt-1)+1:6*(cnt-1)+6),1) = zeros(6,1);
        end

        vACross = [vACross zeros(length(vACross),6); zeros(6,length(vACross))  [hat(vtemp(4:6)), zeros(3); hat(vtemp(1:3)), hat(vtemp(4:6))]];

    end


    % ---------------- STEP 2: Perform SAT test to find closed contacts ---------------- %
    [contact,WNA,WTA,WNM,WTM,cp] = CollisionCheckLCP(AH_B(:,:,:,ii),AH_Bm,obj,c);


    % ---------------- STEP 3: Solve the contact problem for closed contacts ---------------- %
    if  sum(contact>0)
        nn = length(cp(1,:));

        %Construct the diagonal friction matrix
        Mu = diag(repmat(c.mu,1,nn));

        term1 = B_fM*c.dt - vACross*B_M_B*vA*c.dt;

        %Define the normal and tangential velocities at tA and tM. We
        %substranct from the relative velocity the velocity of the
        %contact surface.
        gammaNA=[]; gammaNE=[]; gammaTA=[]; gammaTE=[]; conv = [];
        % for jj=1:length(surf)
            gammaNA = WNA'*vA;
            gammaTA = WTA'*vA;
            % speed =surface{surf(jj)}.transform(1:3,1:3)*surface{surf(jj)}.speed;
            % conv = [conv; repmat(D*speed,sum(indx(:,surf(jj))),1)];
        % end

        %Setup the matrices for the LCP (See Claude's thesis, page 128)
        C =[Gamma(1:nn*c.dimd,1:nn*c.dimd), zeros(nn*c.dimd,nn),  E(1:nn*c.dimd,1:nn);
            zeros(nn,nn*c.dimd),            Sigma(1:nn,1:nn),     zeros(nn,nn);
            -E(1:nn*c.dimd,1:nn)',          Mu,                   Delta(1:nn,1:nn)];

        G =  [WTM'; WNM'; zeros(nn,Nobj*6)];

        %In the vector below, the restitution parameters appear (not sure this is correct)
        % v = [-c.eT*gammaTA+conv; -c.eN*gammaNA; zeros(nn,1)];
        v = [-c.eT*gammaTA; -c.eN*gammaNA; Upsilon(1:nn,1)];
        u = B_M_B*vA-term1;

        A = C + G/B_M_B*G';
        b = G/B_M_B*u-v;

        %Solve the LCP formulation
        part1(ii) = toc;
        tic
        [x,~] = LCP(A,b);
        LCPtime(ii) = toc;

        PT = x(1:c.dimd*nn);
        PN = x(c.dimd*nn+1:c.dimd*nn+nn);
        % PN = x(c.dimd*nn+nn+1:c.dimd*nn+nn+nn);

        vE = vA + B_M_B\(B_fM*c.dt - vACross*B_M_B*vA*c.dt + WNM*PN + WTM*PT);
        BV_AB(:,ii+1) = vE;
    else
        %Update the velocity to the next time step using configuration at tM
        vE = B_M_B\(B_fM*c.dt - vACross*B_M_B*vA*c.dt) + vA;
        BV_AB(:,ii+1) = vE;
        PN = 0;
        PT = [0;0];
    end

    %Complete the time step
    for kk=1:Nobj
        if obj{kk}.dynamics == false
            AH_B(:,:,kk,ii+1) = AH_B(:,:,kk,ii);
        else
            vtemp = BV_AB((6*(kk-1)+1:6*(kk-1)+6),ii+1);
            AH_B(:,:,kk,ii+1)  = AH_Bm(:,:,kk)*expm(0.5*c.dt*hat(vtemp));
        end
    end

    FN(1:length(PN),ii) = PN;
    FT(1:length(PT),ii) = PT;
end

%%
%Plot the trajectory of the box
%Plot the simulation results into a figure
figure(Position=[232,146,1060,620]);

for Iplot=1:20:length(AH_B(1,1,1,:))
    if Iplot ==1
        pause;
    end
    for ib = 1:length(obj)
        plotBox(AH_B(:,:,ib,Iplot),obj{ib},[0.7 0.7 0.7]);
    hold on
    end
    % plotBox(AH_B(:,:,1,Iplot),obj{1},[0.7 0.7 0.7]); hold on;
    % plotBox(AH_B(:,:,2,Iplot),obj{2},[0 0 1]); 
    % plotBox(AH_B(:,:,3,Iplot),obj{3},[0 0 1]); 
    % plotBox(AH_B(:,:,4,Iplot),obj{4},[0.5 0.5 0.5]); 
    % plotBox(AH_B{5}(:,:,Iplot),obj{5},[0 0 1]); 
    % plotBox(AH_B{6}(:,:,Iplot),obj{6},[0 0 1]); 
    % plotBox(AH_B{7}(:,:,Iplot),obj{7},[0 0 1]); 
    grid on; axis equal
    % axis([-2 2 -1 3 -2 2]);   
    axis([-5 5 -5 5 -1 3]);   
    
    view(-44,11);
    % view(36,31)
    % zoom(2)
    
    ax = gca;
    ax.Clipping = "off";
    drawnow

    hold off

    
    % drawnow
    % pause;
    
    
end
