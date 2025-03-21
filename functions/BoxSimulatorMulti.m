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
N = 590;
c.dt = 1/1000;
c.eN = 0.4;
c.eT = 0;
c.mu = 0.4;
c.a = 0.0001;
c.tol = 1e-7;
% Bv_AB = x.releaseLinVel;                      %Linear velocity at t0 of B wrt frame A  [m/s]
% Bomg_AB = x.releaseAngVel;                    %Angular velocity at t0 of B wrt frame A [m/s]
PNfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
LNfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
lnfull = ones(5000,1);       %Initial guess for momenta PN            [N*s]
PTfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
LTfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
ltfull(1:5000,1) = {[1;1]};   %initial guess for momenta PT            [N*s]
% N = round(c.endtime/c.dt);                           %Run to this frame                       [-]

%% Preallocate memory to speed up the process
% FT = NaN(length(box.vertices),N);
% FN = NaN(length(box.vertices),N); 
% AH_B = NaN(4,4,N);
% E  = NaN(1,N);

%% Retrieve info of box struct
% B_M_B = box.B_M_B;
                           
%% Wrench acting on body B: expressed in B with orientation of A:
% BA_fo  = [0; 0; -box.mass*g;];
% BA_Tau = [0; 0; 0];
% BA_f   = [BA_fo; BA_Tau];

%% For testing purposes
for ii = 1:3 %3 objects
    obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
    obj{ii}.mass = 1;
    obj{ii}.initAH_B = [eye(3) [0; 0; 0.6*ii;]; zeros(1,3),1];
    obj{ii}.initBV_AB = zeros(6,1);
    obj{ii}.dim = [0.1*ii;0.2*ii;0.1];
    % obj{ii}.dim = [0.3+0.0001*ii;0.2+0.0001*ii;0.1];
    % obj{ii}.dim = [0.3-0.05*ii;0.2-0.05*ii;0.1];

    %Go over the surfaces
    obj{ii}.surface{1}.transform = [eye(3), [0; 0; obj{ii}.dim(3)/2]; zeros(1,3),1 ];
    obj{ii}.surface{2}.transform = [Ry(90), [obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
    obj{ii}.surface{3}.transform = [Ry(180), [0; 0; -obj{ii}.dim(3)/2]; zeros(1,3),1 ];
    obj{ii}.surface{4}.transform = [Ry(-90), [-obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
    obj{ii}.surface{5}.transform = [Rx(90), [0; -obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
    obj{ii}.surface{6}.transform = [-Rx(90), [0; obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
    obj{ii}.surface{1}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
    obj{ii}.surface{2}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
    obj{ii}.surface{3}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
    obj{ii}.surface{4}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
    obj{ii}.surface{5}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];
    obj{ii}.surface{6}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];    
    obj{ii}.surface{1}.speed = [0; 0; 0];
    obj{ii}.surface{2}.speed = [0; 0; 0];
    obj{ii}.surface{3}.speed = [0; 0; 0];
    obj{ii}.surface{4}.speed = [0; 0; 0];
    obj{ii}.surface{5}.speed = [0; 0; 0];
    obj{ii}.surface{6}.speed = [0; 0; 0];

    %Contact points
    Ndisc=2;
    [X,Y,Z]=meshgrid(linspace(-obj{ii}.dim(1)/2,obj{ii}.dim(1)/2,Ndisc),linspace(-obj{ii}.dim(2)/2,obj{ii}.dim(2)/2,Ndisc),linspace(-obj{ii}.dim(3)/2,obj{ii}.dim(3)/2,Ndisc));
    pbool = (abs(X(:))==obj{ii}.dim(1)/2) | (abs(Y(:))==obj{ii}.dim(2)/2) | (abs(Z(:))==obj{ii}.dim(3)/2);
    obj{ii}.vertices= [X(pbool)';Y(pbool)';Z(pbool)'];

    %Determine if the object has dynamics
    obj{ii}.dynamics = true;
end
% obj{1}.dynamics = false; %first object fixed in space
for ii = 1:length(obj{2}.surface)
    obj{2}.surface{ii}.speed = [1; 1; 0];
end
obj{1}.initBV_AB = [0;0;5;0;0;0];
obj{1}.initAH_B = obj{1}.initAH_B*[Rx(15)*Ry(2) zeros(3,1); zeros(1,3),1];
obj{2}.initAH_B = obj{2}.initAH_B*[Rx(3)*Ry(5) zeros(3,1); zeros(1,3),1];
obj{3}.initAH_B = obj{3}.initAH_B*[Rx(-10)*Ry(-5) zeros(3,1); zeros(1,3),1];

% obj{2}.surface{2}.speed = [1; 1; 0];
% obj{2}.surface{3}.speed = [1; 1; 0];
% obj{2}.surface{4}.speed = [1; 1; 0];
% obj{2}.surface{5}.speed = [1; 1; 0];
% obj{2}.surface{6}.speed = [1; 1; 0];
% obj{3}.dynamics = false; %first object fixed in space

% figure; 
% hold on;
% for ii = 1:3
%     plotBox(obj{ii}.initAH_B,obj{ii},[0.7 0.7 0.7])
%     hold on
%     for jj = 1:6
%         AH_C = obj{ii}.initAH_B*obj{ii}.surface{jj}.transform;
%         tip = [AH_C(1:3,4)+0.1*AH_C(1:3,1) AH_C(1:3,4)+0.1*AH_C(1:3,2) AH_C(1:3,4)+0.1*AH_C(1:3,3)];
%         plot3([AH_C(1,4) tip(1,1)],[AH_C(2,4) tip(2,1)],[AH_C(3,4) tip(3,1)],'r'); hold on
%         plot3([AH_C(1,4) tip(1,2)],[AH_C(2,4) tip(2,2)],[AH_C(3,4) tip(3,2)],'g');
%         plot3([AH_C(1,4) tip(1,3)],[AH_C(2,4) tip(2,3)],[AH_C(3,4) tip(3,3)],'b');
%     end
% end
% axis equal

%define ground plane
ii = 4;

obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
obj{ii}.mass = 1;
obj{ii}.initAH_B = eye(4);
obj{ii}.initBV_AB = zeros(6,1);
obj{ii}.dim = [10;10;0.01];
% obj{ii}.dim = [0.3+0.0001*ii;0.2+0.0001*ii;0.1];
% obj{ii}.dim = [0.3-0.05*ii;0.2-0.05*ii;0.1];

%Go over the surfaces
obj{ii}.surface{1}.transform = [eye(3), [0; 0; obj{ii}.dim(3)/2]; zeros(1,3),1 ];
obj{ii}.surface{2}.transform = [Ry(90), [obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
obj{ii}.surface{3}.transform = [Ry(180), [0; 0; -obj{ii}.dim(3)/2]; zeros(1,3),1 ];
obj{ii}.surface{4}.transform = [Ry(-90), [-obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
obj{ii}.surface{5}.transform = [Rx(90), [0; -obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
obj{ii}.surface{6}.transform = [-Rx(90), [0; obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
obj{ii}.surface{1}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
obj{ii}.surface{2}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
obj{ii}.surface{3}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
obj{ii}.surface{4}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
obj{ii}.surface{5}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];
obj{ii}.surface{6}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];

%Contact points
Ndisc=2;
[X,Y,Z]=meshgrid(linspace(-obj{ii}.dim(1)/2,obj{ii}.dim(1)/2,Ndisc),linspace(-obj{ii}.dim(2)/2,obj{ii}.dim(2)/2,Ndisc),linspace(-obj{ii}.dim(3)/2,obj{ii}.dim(3)/2,Ndisc));
pbool = (abs(X(:))==obj{ii}.dim(1)/2) | (abs(Y(:))==obj{ii}.dim(2)/2) | (abs(Z(:))==obj{ii}.dim(3)/2);
obj{ii}.vertices= [X(pbool)';Y(pbool)';Z(pbool)'];

%Determine if the object has dynamics
obj{ii}.dynamics = false;

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



%% Dynamics
for ii=1:N
    % ---------------- STEP 1: Compute system matrices at mid point ---------------- %
    
    %Obtain the linear and angular velocity at tA
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

    if ii == 568 || ii == 591 || ii == 592 || ii == 593
        hoi = 1;
    end
   
    % ---------------- STEP 2: Perform SAT test to find closed contacts ---------------- %
    [contact,WNA,WTA,WNM,WTM,cp] = CollisionCheck(AH_B(:,:,:,ii),AH_Bm,obj);

    % ---------------- STEP 3: Solve the contact problem for closed contacts ---------------- %
    if  sum(contact>0)
        cpi = length(WNM(1,:));
        
        %Give an initial guess for the normal and tangential momenta
        PN=PNfull(1:cpi);
        PT=cell2mat(PTfull(1:cpi));
        term1 = B_fM*c.dt - vACross*B_M_B*vA*c.dt;
        converged = 0;
        tel = 1;
        while converged==0
            %Discrete time dynamics: Equations of motion
            vE = vA + B_M_B\(term1 + WNM*PN + WTM*PT);
            
            %Define the normal and tangential velocities at tA and tM. We substranct from the relative velocity the velocity of the contact surface.
            gammaNA=[]; gammaNE=[]; gammaTA=[]; gammaTE=[];

            % speed = obj{cidx(2,surf(jj))}.surface{cidx(3,surf(jj))}.speed; %Speed of the plane with which we are in contact
            % speed = zeros(3,1);
            gammaNA = WNA'*vA; %Because we consider only 1 contact point at the time, we can substract directly the velocity without repmat-ing it to the amount of contact points that make contact with that plane
            gammaNE = WNM'*vE;

            %Define the tangential velocities at tA and tM
            gammaTA = WTA'*vA;
            gammaTE = WTM'*vE;
            
            %Newtons restitution law
            xiN = gammaNE+c.eN*gammaNA;
            xiT = gammaTE+c.eT*gammaTA;
            
            %Using prox functions: project PN and PT
            PNold = PN;
            PTold = PT;
            PN = proxCN(PN-c.a*xiN);
            PT = proxCT(PT-c.a*xiT,c.mu*PN);

            %Compute the error
            error = norm(PN-PNold)+norm(PT-PTold);
            
            %Check for convergence
            converged = error<c.tol;
            tel = tel+1;
        end
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
for Iplot=530:1:length(AH_B(1,1,1,:))
    if Iplot ==1
        pause;
    end
    plotBox(AH_B(:,:,1,Iplot),obj{1},[0 0 1]); hold on;
    plotBox(AH_B(:,:,2,Iplot),obj{2},[0 0 1]); 
    plotBox(AH_B(:,:,3,Iplot),obj{3},[0 0 1]); 
    plotBox(AH_B(:,:,4,Iplot),obj{4},[0.5 0.5 0.5]); 
    % plotBox(AH_B{5}(:,:,Iplot),obj{5},[0 0 1]); 
    % plotBox(AH_B{6}(:,:,Iplot),obj{6},[0 0 1]); 
    % plotBox(AH_B{7}(:,:,Iplot),obj{7},[0 0 1]); 
    grid on; axis equal
    axis([-3 3 -3 3 0 3]);    
    hold off
    view(-44,11);
    drawnow
    % pause;
    ax = gca;
    ax.Clipping = "off";
    
end
