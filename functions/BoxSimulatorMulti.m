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
%                                  obj{ii}.B_M_B : 6x6 double intertia tensor of the box  []  
%                                  obj{ii}.mass  : 1x1 double mass of the box             [kg]
%                                  obj{ii}.cp    : 3xN double position of the contact     [m]
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
%% Constants and settings
close all;clearvars;clc;
g     = 9.81;                                 %Gravitational acceleration              [m/s^2]
N = 1500;
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
    obj{ii}.cp= [X(pbool)';Y(pbool)';Z(pbool)'];

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
obj{ii}.cp= [X(pbool)';Y(pbool)';Z(pbool)'];

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
    AH_B{ii}(:,:,1) = obj{ii}.initAH_B;

    % if obj{cnt}.dynamics == false
    %     continue;
    % end
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
        AH_Bm{kk} = AH_B{kk}(:,:,ii)*expm(0.5*c.dt*hat(vtemp));
        AR_Bm{kk} = AH_Bm{kk}(1:3,1:3);
        Ao_Bm{kk} = AH_Bm{kk}(1:3,4);
        AR_B{kk}  = AH_B{kk}(1:3,1:3,ii);


        %Compute the wrench at tM
        if obj{kk}.dynamics
            B_fM((6*(cnt-1)+1:6*(cnt-1)+6),1) = ([AR_Bm{kk} zeros(3); zeros(3) AR_Bm{kk}])'*BA_f((6*(cnt-1)+1:6*(cnt-1)+6),1);
        else
            B_fM((6*(cnt-1)+1:6*(cnt-1)+6),1) = zeros(6,1);
        end

        vACross = [vACross zeros(length(vACross),6); zeros(6,length(vACross))  [hat(vtemp(4:6)), zeros(3); hat(vtemp(1:3)), hat(vtemp(4:6))]];

    end
   
    % ---------------- STEP 3: Compute for closed contacts ---------------- %
    %If a gap function has been closed, contact has been made
    %
    [contact,WNA,WTA,WNM,WTM,cp] = CollisionCheck(AH_B,AH_Bm,obj,ii);
    %
    if ii == 569
        hoi=1;
    end
    if  sum(contact>0)
        cpi = length(WNM(1,:));
        
        %Give an initial guess for the normal and tangential momenta
        PN=PNfull(1:cpi);
        PT=cell2mat(PTfull(1:cpi));
        LN = LNfull(1:cpi); 
        ln = lnfull(1:cpi); 
        LT=cell2mat(LTfull(1:cpi));
        lt=cell2mat(ltfull(1:cpi));
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
        LN = 0;
        ln = 0;
        PT = [0;0];
        LT = [0;0];
        lt = [0;0];
    end
    %Give estimate for PN and PT for next timestep (speeds up convergence
    %in case of persistant contact)
    % if ~sum(contact>0)
    % % if contact
    %     PNfull(cpi)=PN;
    %     LNfull(cpi)=LN;
    %     lnfull(cpi)=ln;
    %     cnt=1;
    %     for jj = length(cpi)
    %         PTfull(cpi(jj)) = {PT(cnt:cnt+1)};
    %         LTfull(cpi(jj)) = {LT(cnt:cnt+1)};
    %         ltfull(cpi(jj)) = {lt(cnt:cnt+1)};
    %         cnt=cnt+2;
    %     end
    % end
    
    %Complete the time step
    for kk=1:Nobj
        if obj{kk}.dynamics == false
            AH_B{kk}(:,:,ii+1) = AH_B{kk}(:,:,ii);
        else
            vtemp = BV_AB((6*(kk-1)+1:6*(kk-1)+6),ii+1);
            AH_B{kk}(:,:,ii+1)  = AH_Bm{kk}*expm(0.5*c.dt*hat(vtemp));
        end        
    end
    
    %Kinetic energy
    % E(ii) = 0.5*BV_AB(:,ii)'*B_M_B*BV_AB(:,ii);  %Kinetic energy

    FN(1:length(PN),ii) = PN;
    FT(1:length(PT),ii) = PT;
end

%%
%Plot the trajectory of the box
%Plot the simulation results into a figure
figure(Position=[232,146,1060,620]);
for Iplot=1:16:length(AH_B{1})
    if Iplot ==1
        pause;
    end
    plotBox(AH_B{1}(:,:,Iplot),obj{1},[0 0 1]); hold on;
    plotBox(AH_B{2}(:,:,Iplot),obj{2},[0 0 1]); 
    plotBox(AH_B{3}(:,:,Iplot),obj{3},[0 0 1]); 
    plotBox(AH_B{4}(:,:,Iplot),obj{4},[0.5 0.5 0.5]); 
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
%%

% %Check if body 1 is in contact with body 2
% close all; clear cp
% 

% Test case 1
AH_Bm{1} =   [1.0000         0         0         0
                   0   -1.0000   -0.0000         0
                   0    0.0000   -1.0000    1.0000
                   0         0         0    1.0000];
AH_Bm{2} = [    1.0000         0         0         0
                     0    1.0000         0         0
                     0         0    1.0000    1.0994
                     0         0         0    1.0000];

AH_Bm{3} =   [1.0000         0         0         0
                   0   -1.0000   -0.0000         0
                   0    0.0000   -1.0000    4.0000
                   0         0         0    1.0000];
% 
% 
% %Test case 2
AH_Bm{2} = AH_Bm{2}*[Rx(3)*Ry(3), [0;0;0.001]; zeros(1,3),1];

% %Test case 3
% AH_Bm{2} = AH_Bm{2}*[Rx(30)*Ry(3), [0;0;+0.055]; zeros(1,3),1];
% % 
%Test case 3
AH_Bm{1} = AH_Bm{1}*[Rx(180), zeros(3,1); zeros(1,3),1];
AH_Bm{2} = AH_Bm{2}*[Rx(5)*Ry(3), [0;0;0.012]; zeros(1,3),1];
% 
% Test case 4: Edge-edge
% AH_Bm{2} = AH_Bm{2}*[Rx(5)*Ry(3), [0;0.15;0.015]; zeros(1,3),1];
% 
figure;
% plotBox(AH_B{1}(:,:,ii),obj{1},[0 0 1]); hold on;
% plotBox(AH_B{2}(:,:,ii),obj{2},[0 0 1]); 
% plotBox(AH_B{3}(:,:,ii),obj{3},[0 0 1]);  
plotBox(AH_Bm{1}(:,:),obj{1},[0 0 1]); hold on;
plotBox(AH_Bm{2}(:,:),obj{2},[0 0 1]); 
% plotBox(AH_Bm{3}(:,:),obj{3},[0 0 1]); 
grid on; axis equal
axis off
axis([-0.4 0.5 -1 1 0 3.5]);    
hold off
view(-2,3);
drawnow
ax = gca;
ax.Clipping = "off";
% 
tic
[contact,WNA,WTA,WNM,WTM,cp] = CollisionCheck(AH_B,AH_Bm,obj,ii);
toc
% 
% %Contact points in world frame used for plotting
ibref = 1;
Acp = cp;
hold on; plot3(Acp(1,:),Acp(2,:),Acp(3,:),'.','color',[1 0 0],'markersize',30)

%%
function [contact,WNA,WTA,WNM,WTM,cp] = CollisionCheck(AH_B,AH_Bm,obj,JJ)
%Function to check if there is a collision between bodies. If there is, it will compute the matrices
%containing the force directions based on a set of computed contact points. 

Nobj = length(obj);     %Amount of objects in the scene

WNM = zeros(6*Nobj,1); %Preallocate: each contact point is one column with 6xNobj rows
WTM = zeros(6*Nobj,2); %Preallocate: each contact point is two columns with 6xNobj rows
WNA = zeros(6*Nobj,1); %Preallocate: each contact point is one column with 6xNobj rows
WTA = zeros(6*Nobj,2); %Preallocate: each contact point is two columns with 6xNobj rows
cp = [];

%If we have only 1 body in the scene, nothing will happen. 
if Nobj <2    
    return;
end

%Create the vector of objects we need to check contact for
%each body needs to be checked with each other body only once
vecObj = 1:Nobj;
vecCheck = [];
for ii = 1:length(vecObj)-1
    vecCheck = [vecCheck [repmat(vecObj(ii),1,vecObj(end-ii)); vecObj(ii+1:end)]]; 
end

%Now we are going to build the matrices of the force direction
cntr = 1; %To count the contact points (and hence the size of the matrices containing the force directions)
for II = 1:length(vecCheck) 
    contact(II) = false;

    %Get the body indexes we are going to check for
    ib1 = vecCheck(1,II); %Body 1 index
    ib2 = vecCheck(2,II); %Body 2 index

    %The 15 vectors we need to check for SAT test
    b1v1 = AH_Bm{ib1}(1:3,1); b1v2 = AH_Bm{ib1}(1:3,2); b1v3 = AH_Bm{ib1}(1:3,3);
    b2v1 = AH_Bm{ib2}(1:3,1); b2v2 = AH_Bm{ib2}(1:3,2); b2v3 = AH_Bm{ib2}(1:3,3);
    
    SATvec = [b1v1/norm(b1v1) b1v2/norm(b1v2) b1v3/norm(b1v3) b2v1/norm(b2v1) b2v2/norm(b2v2) b2v3/norm(b2v3) cross(b1v1,b2v1)/norm(cross(b1v1,b2v1)) cross(b1v1,b2v2)/norm(cross(b1v1,b2v2)) cross(b1v1,b2v3)/norm(cross(b1v1,b2v3))...
        cross(b1v2,b2v1)/norm(cross(b1v2,b2v1)) cross(b1v2,b2v2)/norm(cross(b1v2,b2v2)) cross(b1v2,b2v3)/norm(cross(b1v2,b2v3)) cross(b1v3,b2v1)/norm(cross(b1v3,b2v1)) cross(b1v3,b2v2)/norm(cross(b1v3,b2v2)) cross(b1v3,b2v3)/norm(cross(b1v3,b2v3))];
    SATidx = [1 2 3 0 0 0 1 1 1 2 2 2 3 3 3; 0 0 0 1 2 3 1 2 3 1 2 3 1 2 3]; %Tracking which vector is used for each body (1=x, 2=y, 3=z)
    
    
    ha = AH_Bm{ib1}(1:3,1:3)*(obj{ib1}.dim/2); %According to Erleben book this should be defined in local frame (half width of the object in local frame)
    hb = AH_Bm{ib2}(1:3,1:3)*(obj{ib2}.dim/2);
    dis = (AH_Bm{ib2}(1:3,4)-AH_Bm{ib1}(1:3,4));
    
    for ii = 1:length(SATvec)
        s1(ii) = abs(dis'*SATvec(:,ii)) - ( abs((ha/norm(ha))'*SATvec(:,ii)) + abs((hb/norm(hb))'*SATvec(:,ii)) ); %Equation 86 of Erleben book
        s(ii) = abs(dis'*SATvec(:,ii)) - ( abs(min((AH_Bm{ib1}(1:3,1:3)*obj{ib1}.cp)'*SATvec(:,ii))) + abs(min((AH_Bm{ib2}(1:3,1:3)*obj{ib2}.cp)'*SATvec(:,ii))) );
    end
    
    
    %Check if we have contact. In that case, they all must show overlap (there exist no axis that shows separation)
    contact(II) = sum(s<=0)==15;
    
    if contact(II)
        %First need to compensate for a scenario with equal orientation
        % if sum(abs(SATvec(:,7)))<1e-10 || sum(abs(SATvec(:,11)))<1e-10 || sum(abs(SATvec(:,15)))<1e-10 %If there is a (numerically) zero vector in the SATvec, at least one of the cross products must be zero
        if vecnorm(SATvec(:,7))<0.01 || vecnorm(SATvec(:,11))<0.01 || vecnorm(SATvec(:,15))<0.01 %If there is a (numerically) zero vector in the SATvec, at least one of the cross products must be zero
            %In this case we only consider the first six vectors, as we are dealing with a face-face collision
            isat = find(s(1:6)==max(s(1:6))); %Take the max of s in the first 6 vectors
            isat = isat(1); %If orientations are exactly the same, we choose the first body as the reference body
        else
            %Vector with minimum overlap is the one where the impact is happening
            isat = find(s==max(s)); %Take the max of s in all vectors
        end    
    
        %If the collision along the first 3 axis, the first body has the reference face
        %If the collision is along the next 3 axis, the second body has the reference face
        %If the collision is along one of the other 9 axes, we have a edge-edge collision, and we need to retreive the axis
        if isat<4
            %Body 1 has the reference face
            ibref = ib1;
            ibinc = ib2;
        elseif isat <7
            %Body 2 has the reference face
            ibref = ib2;
            ibinc = ib1;
        else
            %We have an edge-edge collision
            v1 = AH_Bm{ib1}(1:3,SATidx(1,isat)); %Using the index vector to get back which vector was used to create the cross product that is now giving the smallest overlap
            v2 = AH_Bm{ib2}(1:3,SATidx(2,isat));
            % normal = cross(v1,v2)/norm(cross(v1,v2)); % cross(v1,v2) .. I guess? no.. need to take the edge vectors and map them to A.. but edge vectors would be the same as the frame axis vectors.. already expressed in A.. so should be good!
        end
        normal = SATvec(:,isat);  
        
        if isat <7
            %Here we have a face intersection with a feature from the other box
            
            %Here we have a face intersection with a feature from the other box
            %We first need to get the normal of the surface of the reference box that is in contact
            % This is the surface on the box that is closest to the other box, along the normal that we used to compute it    
            dis1 = (AH_Bm{ibref}(1:3,4)-AH_Bm{ibinc}(1:3,4))'*SATvec(:,isat); %distance from incident body to the reference body projected on separating axis
            for ii = 1:length(obj{ibref}.surface)                
                dis2 = (AH_Bm{ibref}(1:3,1:3)*obj{ibref}.surface{ii}.transform(1:3,4))'*SATvec(:,isat); %distance from the reference body frame to its surfaces projected on the separating axes
                f(ii) = dis1+dis2; %face that is closest to the incident body is the reference face
            end

            fidx = find(abs(f)==min(abs(f))); %Face index of reference box that is closest to incident box
            AH_C = AH_Bm{ibref}*obj{ibref}.surface{fidx}.transform; %The transform of that surface is our surface transform
    
            %Get the pose of incident body expressed in reference body
            B1H_B2 = AH_Bm{ibref}\AH_Bm{ibinc};

            %Get the pose of incident body expressed in reference face frame
            CH_B2 = AH_C\AH_Bm{ibinc};
    
            %Locations of the vertices of incident body expressed in reference body
            B1p2 = B1H_B2(1:3,1:3)*obj{ibinc}.cp +B1H_B2(1:3,4);
            Cp2  = CH_B2(1:3,1:3)*obj{ibinc}.cp +CH_B2(1:3,4);
            B1p1 = obj{ibref}.cp;
    
            %Find the vertices that lie below reference plane
            % p_i = find((abs(B1p2(3,:))<abs(B1p1(3)))==1);
            p_i = find(Cp2(3,:)<0); %Vertices that lie below reference plane
    
            for ii = 1:length(p_i)
                for id = 1:2
                    if abs(B1p2(id,p_i(ii)))<(obj{ibref}.dim(id)/2)
                        B1cp(id,ii) = B1p2(id,p_i(ii));
                        B1cp(3,ii) = B1p2(3,p_i(ii));
                    else
                        index = find(vecnorm(B1p2(:,p_i(ii))-B1p1)==min(vecnorm(B1p2(:,p_i(ii))-B1p1))); %Find the contact point on body 1 closest to body 2
                        B1cp(id,ii) = B1p1(id,index);
                        B1cp(3,ii) = B1p1(3,index);
                    end
                end
            end
            %Express the contact points in incident body
            cp = AH_Bm{ibref}(1:3,1:3)*B1cp + AH_Bm{ibref}(1:3,4);
            B2cp = AH_Bm{ibinc}(1:3,1:3)\(cp - AH_Bm{ibinc}(1:3,4)); %B2H_B1(1:3,1:3)*B1cp + B2H_B1(1:3,4); %its the incident body's contact points that we consider

            %Now we need to create the matrices containing the force directions
            for icp = 1:length(B1cp(1,:)) %Count over the contact points between these two bodies
                % AR_C = [null(normal') normal]; % We assume the convention of having the normal of the reference face as positive
                AR_C = AH_C(1:3,1:3);
                wB1 = (AR_C'*[-AH_Bm{ibref}(1:3,1:3) AH_Bm{ibref}(1:3,1:3)*hat(B1cp(:,icp))])';     %This is the reaction force, in oposite direction 
                wB1a= (AR_C'*[-AH_B{ibref}(1:3,1:3,JJ) AH_B{ibref}(1:3,1:3,JJ)*hat(B1cp(:,icp))])'; %This is the reaction force, in oposite direction 
                wB2 = (AR_C'*[AH_Bm{ibinc}(1:3,1:3) -AH_Bm{ibinc}(1:3,1:3)*hat(B2cp(:,icp))])';     %Starndard convention for positive normal
                wB2a= (AR_C'*[AH_B{ibinc}(1:3,1:3,JJ) -AH_B{ibinc}(1:3,1:3,JJ)*hat(B2cp(:,icp))])'; %Starndard convention for positive normal

                if obj{ibref}.dynamics
                    WNM((6*(ibref-1)+1:6*(ibref-1)+6),cntr) = wB1(:,3); %For the ibref body we compute the force directions
                    WTM((6*(ibref-1)+1:6*(ibref-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1(:,1:2); %Assuming this works the same way in the tangential direction

                    WNA((6*(ibref-1)+1:6*(ibref-1)+6),cntr) = wB1a(:,3); %For the ibref body we compute the force directions
                    WTA((6*(ibref-1)+1:6*(ibref-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1a(:,1:2); %Assuming this works the same way in the tangential direction
                end

                if obj{ibinc}.dynamics
                    WNM((6*(ibinc-1)+1:6*(ibinc-1)+6),cntr) = wB2(:,3); %For the ibinc body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                    WTM((6*(ibinc-1)+1:6*(ibinc-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2(:,1:2);

                    WNA((6*(ibinc-1)+1:6*(ibinc-1)+6),cntr) = wB2a(:,3); %For the ibinc body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                    WTA((6*(ibinc-1)+1:6*(ibinc-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2a(:,1:2);
                end

                cntr = cntr+1;
            end
        else
            %Here we have an edge-edge collision and only a single contact point
            %Find the edges that are colliding
            for jj = 1:2 %loop over the two bodies
                pbool = (abs(obj{vecCheck(jj,II)}.cp(1,:))==max(abs(obj{vecCheck(jj,II)}.cp(1,:))))&(abs(obj{vecCheck(jj,II)}.cp(2,:))==max(abs(obj{vecCheck(jj,II)}.cp(2,:))))&(abs(obj{vecCheck(jj,II)}.cp(3,:))==max(abs(obj{vecCheck(jj,II)}.cp(3,:))));
                % Ap = AH_Bm{vecCheck(jj,II)}(1:3,1:3)*obj{vecCheck(jj,II)}.cp(:,pbool); %+AH_Bm{vecCheck(jj,II)}(1:3,4);
                Ap = obj{vecCheck(jj,II)}.cp(:,pbool); %So now we are just considering edges in body frame

                E(:,:,1) = [Ap(:,1) Ap(:,2)];
                E(:,:,2) = [Ap(:,5) Ap(:,6)];
                E(:,:,3) = [Ap(:,7) Ap(:,8)];
                E(:,:,4) = [Ap(:,3) Ap(:,4)];
                E(:,:,5) = [Ap(:,1) Ap(:,3)];%
                E(:,:,6) = [Ap(:,2) Ap(:,4)];%
                E(:,:,7) = [Ap(:,5) Ap(:,7)];
                E(:,:,8) = [Ap(:,6) Ap(:,8)];
                E(:,:,9) = [Ap(:,1) Ap(:,5)];
                E(:,:,10) = [Ap(:,2) Ap(:,6)];
                E(:,:,11) = [Ap(:,3) Ap(:,7)];
                E(:,:,12) = [Ap(:,4) Ap(:,8)];

                edges = [E(:,2,1)-E(:,1,1),E(:,2,2)-E(:,1,2),E(:,2,3)-E(:,1,3),E(:,2,4)-E(:,1,4),E(:,2,5)-E(:,1,5),E(:,2,6)-E(:,1,6),E(:,2,7)-E(:,1,7),E(:,2,8)-E(:,1,8),E(:,2,9)-E(:,1,9),E(:,2,10)-E(:,1,10),E(:,2,11)-E(:,1,11),E(:,2,12)-E(:,1,12) ];

                %Find the edges that are parallel to the normal used for the cross product (v1 or v2, their cross product with the normal will be zero)

                if jj ==1
                    Aedges = AH_Bm{ib1}(1:3,1:3)*(edges./vecnorm(edges)); %Edges expressed in world frame and normalized
                    for ii = 1:length(Aedges)
                        ed{jj}(ii) = norm(cross(v1,Aedges(:,ii)));
                    end
                    % idx = find(~round(ed{jj})); %Indices of edges that are parallel to the axis used for computing axis of separation
                    idx = find(ed{jj}==min(ed{jj})); %Safer way of computing it
                    B2H_B1 = AH_Bm{ib2}\AH_Bm{ib1}; %Transformation from Body 1 to Body 2
                    B2p1 = B2H_B1(1:3,1:3)*reshape(E(:,:,idx),3,8) +B2H_B1(1:3,4); %Expressing the edges of body 1 in body 2
                    dis2 = [mean(B2p1(:,1:2),2), mean(B2p1(:,3:4),2), mean(B2p1(:,5:6),2), mean(B2p1(:,7:8),2)]; %Vectors from body 2 to the edges of body 1
                    dis3 = AH_Bm{ib2}(1:3,1:3)*dis2; %Expressing these vectors in frame A
                    dis4 = abs(dis3'*SATvec(:,isat)); %These vectors projected on the separating axis
                    idx2 = find(dis4==min(dis4)); %Index of parallel edge closest to other body
                    eA = AH_Bm{ib1}(1:3,1:3)*edges(:,idx(idx2)); %Edge eA from body A closest to body B
                    P1 = AH_Bm{ib1}(1:3,1:3)*E(:,1,idx(idx2))+AH_Bm{ib1}(1:3,4);
                    V1 = AH_Bm{ib1}(1:3,1:3)*edges(:,idx(idx2));
                elseif jj==2
                    Aedges = AH_Bm{ib2}(1:3,1:3)*(edges./vecnorm(edges)); %Edges expressed in world frame and normalized
                    for ii = 1:length(Aedges)
                        ed{jj}(ii) = norm(cross(v2,Aedges(:,ii)));
                    end
                    % idx = find(~round(ed{jj})); %Indices of edges that are parallel to the axis used for computing axis of separation
                    idx = find(ed{jj}==min(ed{jj})); %Safer way of computing it
                    B1H_B2 = AH_Bm{ib1}\AH_Bm{ib2}; %Transformation from Body 2 to Body 1
                    B1p2 = B1H_B2(1:3,1:3)*reshape(E(:,:,idx),3,8) +B1H_B2(1:3,4); %Expressing the edges of body 2 in body 1
                    dis2 = [mean(B1p2(:,1:2),2), mean(B1p2(:,3:4),2), mean(B1p2(:,5:6),2), mean(B1p2(:,7:8),2)]; %Distance from body 1 to the edges of body 2
                    dis3 = AH_Bm{ib1}(1:3,1:3)*dis2; %Expressing these vectors in frame A
                    dis4 = abs(dis3'*SATvec(:,isat)); %These vectors projected on the separating axis
                    idx2 = find(dis4==min(dis4)); %Index of parallel edge closest to other body
                    eB = AH_Bm{ib2}(1:3,1:3)*edges(:,idx(idx2)); %Edge eA from body A closest to body B
                    P2 = AH_Bm{ib2}(1:3,1:3)*E(:,1,idx(idx2))+AH_Bm{ib2}(1:3,4);
                    V2 = AH_Bm{ib2}(1:3,1:3)*edges(:,idx(idx2));
                end
            end
            normal = cross(eA,eB)/norm(cross(eA,eB)); %Erleben eq. 87
            %Alternative, we know that the normal is the same as the separating axis SATvec(:,isat) ... 
            % normal = SATvec(:,isat);
            b = s(isat)*SATvec(:,isat)-P2+P1;
            A = [-V1(:,1) V2(:,1)]; %Taking first column. If there are more columns, their values must be the same, as it means the distance to those edges are the same.
            t = linsolve(A,b);
            L1 = P1+t(1)*V1; %Contact point on body 1 expressed in world frame
            L2 = P2+t(2)*V2; %Contact point on body 1 expressed in world frame
            cp1 = 0.5*(L1+L2); %Single contact point in midpoint of overlapping region expressed in world frame.
            B1cp = AH_Bm{ib1}(1:3,1:3)'*(cp1 - AH_Bm{ib1}(1:3,4)); %Single contact point exprssed in body 1
            cp = AH_Bm{ib1}(1:3,1:3)*B1cp + AH_Bm{ib1}(1:3,4);
            B2cp = AH_Bm{ib2}(1:3,1:3)\(cp - AH_Bm{ib2}(1:3,4)); %its the incident body's contact points that we consider

            %If the projected contact point distance in body 1 is negative on the axis of separation, we need to flip the axis
            if (AH_Bm{ib1}(1:3,1:3)*B1cp)'*normal <0
                normal = -normal;
            end

            %Now we need to create the matrices containing the force directions
            for icp = 1:length(B1cp(1,:)) %Count over the contact points between these two bodies
                AR_C = [null(normal') normal]; % We assume the convention of having the normal of the reference face as positive
                wB1 = (AR_C'*[-AH_Bm{ib1}(1:3,1:3) AH_Bm{ib1}(1:3,1:3)*hat(B1cp(:,icp))])';     %This is the reaction force, in oposite direction 
                wB1a= (AR_C'*[-AH_B{ib1}(1:3,1:3,JJ) AH_B{ib1}(1:3,1:3,JJ)*hat(B1cp(:,icp))])'; %This is the reaction force, in oposite direction 
                wB2 = (AR_C'*[AH_Bm{ib2}(1:3,1:3) -AH_Bm{ib2}(1:3,1:3)*hat(B2cp(:,icp))])';     %Starndard convention for positive normal
                wB2a= (AR_C'*[AH_B{ib2}(1:3,1:3,JJ) -AH_B{ib2}(1:3,1:3,JJ)*hat(B2cp(:,icp))])'; %Starndard convention for positive normal
                
                if obj{ib1}.dynamics
                    WNM((6*(ib1-1)+1:6*(ib1-1)+6),cntr) = wB1(:,3); %For the ib1 body we compute the force directions
                    WTM((6*(ib1-1)+1:6*(ib1-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1(:,1:2); %Assuming this works the same way in the tangential direction                
    
                    WNA((6*(ib1-1)+1:6*(ib1-1)+6),cntr) = wB1a(:,3); %For the ib1 body we compute the force directions
                    WTA((6*(ib1-1)+1:6*(ib1-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1a(:,1:2); %Assuming this works the same way in the tangential direction  
                end

                if obj{ib2}.dynamics
                    WNM((6*(ib2-1)+1:6*(ib2-1)+6),cntr) = wB2(:,3); %For the ib2 body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                    WTM((6*(ib2-1)+1:6*(ib2-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2(:,1:2);

                    WNA((6*(ib2-1)+1:6*(ib2-1)+6),cntr) = wB2a(:,3); %For the ib2 body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                    WTA((6*(ib2-1)+1:6*(ib2-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2a(:,1:2);
                end

                cntr = cntr+1;
            end
        end
    end
end
end





%% Matrix with force directions
function [WN,WT] = CompW(AR_B,AR_C,cp,cidx,Nobj,obj)
% Compute the matrix containing the tangential force directions.
%Knowning the total amount of object Nobj, we can construct the column of WN and WT
WN = zeros(6*Nobj,1);
WT = zeros(6*Nobj,2);

tel = 1;
for ii = 1:length(cp(1,:)) %For now we leave this, but essentially this should be 1 (because we call this function for every contact point now)
    w = (AR_C'*[AR_B -AR_B*hat(cp(:,ii))])';
    % if obj{cidx(1)}.dynamics
        %Only if it has dynamics
        WN((6*(cidx(1)-1)+1:6*(cidx(1)-1)+6),ii) = w(:,3); %For the cidx(1) body we compute the force directions
        WT((6*(cidx(1)-1)+1:6*(cidx(1)-1)+6),tel:tel+1) = w(:,1:2); %Assuming this works the same way in the tangential direction
    % end

    % if obj{cidx(2)}.dynamics
        WN((6*(cidx(2)-1)+1:6*(cidx(2)-1)+6),ii) = -w(:,3); %For the cidx(2) body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)    
        WT((6*(cidx(2)-1)+1:6*(cidx(2)-1)+6),tel:tel+1) = -w(:,1:2);
    % end
    tel = tel+2;
end
end

%% Proximal point Normal
function y=proxCN(x)
% Proximal point formulation for CN. See thesis for reference.
% prox_CN(x) = 0 for x=< 0
%            = x for x > 0
y=max(x,0);
end

%% Proximal point Tangential
function y=proxCT(x,a)
% Proximal point formulation for CT. See thesis for reference.
%
% prox_CT(x) = x           for ||x|| =< a
%            = a*(x/||x||) for ||x||  > a
% for CT = {x in R^n| ||x|| =< a}
cnt = 1;
for ii = 1:length(a) %For each point in contact
    if norm(x(cnt:cnt+1)) <= a(ii)
        y(cnt:cnt+1,1) = x(cnt:cnt+1); %Stick
    else
        y(cnt:cnt+1,1) = a(ii)*x(cnt:cnt+1)/norm(x(cnt:cnt+1)); %Slip
    end
    cnt = cnt+2;
end
end

function res = hat(vec)
% Take the 3- or 6-vector representing an isomorphism of so(3) or se(3) and
% writes this as element of so(3) or se(3). 
%
% INPUTS:    vec     : 3- or 6-vector. Isomorphism of so(3) or se(3)
%
% OUTPUTS:   res     : element of so(3) or se(3)
%
% Hat operator
if length(vec) == 3
    res = [0, -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];
elseif length(vec) == 6
    skew = [0, -vec(6), vec(5); vec(6), 0, -vec(4); -vec(5), vec(4), 0];
    v = [vec(1);vec(2);vec(3)];
    res = [skew, v; zeros(1,4)];
end
end

function res = vee(mat)
% Takes an element of so(3) or se(3) and returns its isomorphism in R^n.
%
% INPUTS:    mat     : element of so(3) or se(3)
%
% OUTPUTS:   res     : 3- or 6-vector. Isomorphism of so(3) or se(3)
%
% Vee operator

xi1 = (mat(3,2)-mat(2,3))/2;
xi2 = (mat(1,3)-mat(3,1))/2;
xi3 = (mat(2,1)-mat(1,2))/2;

if length(mat) == 3
   res = [xi1; xi2; xi3];
elseif length(mat) == 4
   res = [mat(1:3,4);xi1;xi2;xi3]; 
end
end

function R = Rx(th)
%Rotate around x with th degrees;
R = [1 0 0; 0 cos(deg2rad(th)) -sin(deg2rad(th)); 0 sin(deg2rad(th)) cos(deg2rad(th))];
end

function R = Ry(th)
%Rotate around y with th degrees;
R = [cos(deg2rad(th)) 0 sin(deg2rad(th)); 0 1 0;  -sin(deg2rad(th)) 0 cos(deg2rad(th))];
end

function R = Rz(th)
%Rotate around z with th degrees;
R = [cos(deg2rad(th)) -sin(deg2rad(th)) 0; sin(deg2rad(th)) cos(deg2rad(th)) 0; 0 0 1];
end

function Bplot = plotBox(AH_B,obj,color) 
% Box-simulator-FixedPoint:
%This script is used to plot the box given a certain state
%
% INPUTS:    AH_B  : 4x4 double, pose of the box
%            box   : struct, with fields of box properties as
%                    box.B_M_B  : 6x6 double intertia tensor of the box
%                    box.mass    : 1x1 double mass of the box
%                    box.vertices: 3x8 double position of the vertices of 
%                                  the box w.r.t body-fixed frame
%            color : 3x1 double, rgb color of the box
%
% OUTPUTS:   Bplot : Plot of the box
% Plot the box
        AR_B = AH_B(1:3,1:3);
        %Output the position of the current time step for plotting purposes
        q(:,1)  = AH_B(1:3,4);
        R1(:,1) = AH_B(1:3,1);
        R2(:,1) = AH_B(1:3,2);
        R3(:,1) = AH_B(1:3,3);
        
        %Plot the origin of the box with its unit vectors
        tip = [q(:,1)+ 0.3*R1(:,1) q(:,1)+ 0.3*R2(:,1) q(:,1)+ 0.3*R3(:,1)];
        plot3([q(1,1) tip(1,1)],[q(2,1) tip(2,1)],[q(3,1) tip(3,1)],'r'); hold on
        plot3([q(1,1) tip(1,2)],[q(2,1) tip(2,2)],[q(3,1) tip(3,2)],'g');
        plot3([q(1,1) tip(1,3)],[q(2,1) tip(2,3)],[q(3,1) tip(3,3)],'b');
        
        %Create the box
        pbool = (abs(obj.cp(1,:))==max(abs(obj.cp(1,:))))&(abs(obj.cp(2,:))==max(abs(obj.cp(2,:))))&(abs(obj.cp(3,:))==max(abs(obj.cp(3,:))));
        Ap = AR_B*obj.cp(:,pbool)+AH_B(1:3,4);
%         Ap = AR_B*obj.cp+AH_B(1:3,4);
        Ap_1 = Ap(:,1);
        Ap_2 = Ap(:,2);
        Ap_3 = Ap(:,6);
        Ap_4 = Ap(:,5);
        Ap_5 = Ap(:,3);
        Ap_6 = Ap(:,4);
        Ap_7 = Ap(:,8);
        Ap_8 = Ap(:,7);
        
        plot3([Ap_1(1) Ap_2(1)],[Ap_1(2) Ap_2(2)],[Ap_1(3) Ap_2(3)],'k');hold on;
        plot3([Ap_2(1) Ap_3(1)],[Ap_2(2) Ap_3(2)],[Ap_2(3) Ap_3(3)],'k');%
        plot3([Ap_3(1) Ap_4(1)],[Ap_3(2) Ap_4(2)],[Ap_3(3) Ap_4(3)],'k');
        plot3([Ap_4(1) Ap_1(1)],[Ap_4(2) Ap_1(2)],[Ap_4(3) Ap_1(3)],'k');
        plot3([Ap_5(1) Ap_6(1)],[Ap_5(2) Ap_6(2)],[Ap_5(3) Ap_6(3)],'k');%
        plot3([Ap_6(1) Ap_7(1)],[Ap_6(2) Ap_7(2)],[Ap_6(3) Ap_7(3)],'k');%
        plot3([Ap_7(1) Ap_8(1)],[Ap_7(2) Ap_8(2)],[Ap_7(3) Ap_8(3)],'k');
        plot3([Ap_8(1) Ap_5(1)],[Ap_8(2) Ap_5(2)],[Ap_8(3) Ap_5(3)],'k');
        plot3([Ap_1(1) Ap_5(1)],[Ap_1(2) Ap_5(2)],[Ap_1(3) Ap_5(3)],'k');
        plot3([Ap_2(1) Ap_6(1)],[Ap_2(2) Ap_6(2)],[Ap_2(3) Ap_6(3)],'k');
        plot3([Ap_3(1) Ap_7(1)],[Ap_3(2) Ap_7(2)],[Ap_3(3) Ap_7(3)],'k');
        plot3([Ap_4(1) Ap_8(1)],[Ap_4(2) Ap_8(2)],[Ap_4(3) Ap_8(3)],'k');
        
        %Color the surfaces of the obj
        Bplot = fill3([Ap_1(1) Ap_2(1) Ap_6(1) Ap_5(1)],[Ap_1(2) Ap_2(2) Ap_6(2) Ap_5(2)],[Ap_1(3) Ap_2(3) Ap_6(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face F
        fill3([Ap_1(1) Ap_2(1) Ap_3(1) Ap_4(1)],[Ap_1(2) Ap_2(2) Ap_3(2) Ap_4(2)],[Ap_1(3) Ap_2(3) Ap_3(3) Ap_4(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face A
        fill3([Ap_8(1) Ap_7(1) Ap_6(1) Ap_5(1)],[Ap_8(2) Ap_7(2) Ap_6(2) Ap_5(2)],[Ap_8(3) Ap_7(3) Ap_6(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face C
        fill3([Ap_8(1) Ap_7(1) Ap_3(1) Ap_4(1)],[Ap_8(2) Ap_7(2) Ap_3(2) Ap_4(2)],[Ap_8(3) Ap_7(3) Ap_3(3) Ap_4(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face E
        fill3([Ap_1(1) Ap_4(1) Ap_8(1) Ap_5(1)],[Ap_1(2) Ap_4(2) Ap_8(2) Ap_5(2)],[Ap_1(3) Ap_4(3) Ap_8(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face D
        fill3([Ap_2(1) Ap_3(1) Ap_7(1) Ap_6(1)],[Ap_2(2) Ap_3(2) Ap_7(2) Ap_6(2)],[Ap_2(3) Ap_3(3) Ap_7(3) Ap_6(3)],1,'FaceColor',color,'FaceAlpha',0.8);%Face B

        %Plot all contact points
        % vertices = AR_B*obj.cp+AH_B(1:3,4);
        % plot3(vertices(1,:),vertices(2,:),vertices(3,:),'.','MarkerSize',15,'color',[1 0 0])
        
end