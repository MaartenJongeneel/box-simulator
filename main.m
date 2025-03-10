close all; clc; warning ('off','all'); addpath('readyaml'); addpath('scenes'); addpath('functions');
clearvars;
%This script is used to simulate a box that is tossed on a surface. In the
%settings below, one can set different settings, such as the size of the
%box or the coefficients of friction, tangential and normal restitution.
%% General settings
dosave             = false;         %Save the trajectory (AH_B) to a .mat file
doPlot             = true;          %Show the trajectory of the box
MakeVideo          = false;         %Save the simulation result to video
%% Read the scene that you want to run
scenefile = "IFZ.yml";
% scenefile = "DoubleConveyor.yml";
data = readyaml(scenefile);
%% Parameters for input
c.a                  = 0.001;           %Prox point auxilary parameter             [-]
c.tol                = 1e-7;            %Error tol for fixed-point                 [-]
c.m                  = 1;               %Mass of the box                           [kg]  
c.endtime            = 5;             %Runtime of the simulation                 [s]
c.dt                 = 1/100;          %Timestep at which the simulator runs      [s]
c.dimd               = 4;              %Discretization of the friction cone
step                 = ceil(1/c.dt/50); %Number of discrete points we skip per shown frame
%% Read the scene data
x.releaseOrientation = data.box.release.orientation;  %Release orientation of the box            [deg]
x.releasePosition    = data.box.release.position';    %Release position of the box               [m]
x.releaseLinVel      = data.box.release.linVel';      %Release linear velocity (expressed in B)  [m/s]
x.releaseAngVel      = data.box.release.angVel';      %Release angular velocity (expressed in B) [rad/s]
c.eN                 = data.box.parameters.eN;        %Normal coefficient of restitution         [-]
c.eT                 = data.box.parameters.eT;        %Tangential coefficient of restitution     [-]s
c.mu                 = data.box.parameters.mu;        %Coefficient of friction                   [-]
box                  = data.box;                      %Obtain the box struct
box.B_M_B            = data.box.inertia_tensor;       %Rewrite inertia tensor
surface              = data.surface;                  %Obtain the surfaces 
%% Create the box struct
%Discretization of the box vertices
Ndisc=data.box.discretization;
[X,Y,Z]=meshgrid(linspace(-box.dimensions(1)/2,box.dimensions(1)/2,Ndisc),linspace(-box.dimensions(2)/2,box.dimensions(2)/2,Ndisc),linspace(-box.dimensions(3)/2,box.dimensions(3)/2,Ndisc));
pbool = (abs(X(:))==box.dimensions(1)/2) | (abs(Y(:))==box.dimensions(2)/2) | (abs(Z(:))==box.dimensions(3)/2);
box.vertices= [X(pbool)';Y(pbool)';Z(pbool)'];
% box.vertices = box.vertices;


%% Define the impact planes
for jj = 1:length(surface)
    surface{jj}.speed = surface{jj}.speed';
end

%% Create IFZ scene
clear surface;

%MultiBelt Junction 60deg
nbelts = 8; %number of multibelts
Length = 1.050; %longest belt length
width = 0.760; %total multibelt width
deg = 60;      %multibelt angle
vel = 2.5;  
transform = [eye(3), [0; 2.552/2; 0]; zeros(1,3), 1];

%Create the multibelt
surface = multibelt(nbelts,Length,width,deg,vel,transform,1); 

%Throwbelt
surface{9}.dim = [0.760 2.550];
surface{9}.speed = [0;1.5;0];
surface{9}.transform = [eye(4)];
surface{9}.id = "TB01";

%Catchbelt
surface{10}.dim = [1.450 2.100];
surface{10}.speed = [0;1.5;0];
surface{10}.transform = [Rz(60), [0.30; 2.85; 0]; zeros(1,3),1];
surface{10}.id = "CB01";

%Shortbelt 1
surface{11}.dim = [1.450 0.6];
surface{11}.speed = [0; 1.75; 0];
surface{11}.transform = surface{10}.transform * [eye(3) [0; 1.350; 0]; zeros(1,3), 1];
surface{11}.id = "SB01";

%Shortbelt 2
surface{12}.dim = [1.450 0.6];
surface{12}.speed = [0; 2.0; 0];
surface{12}.transform = surface{10}.transform * [eye(3) [0; 1.950; 0]; zeros(1,3), 1];
surface{12}.id = "SB02";

%Shortbelt 3
surface{13}.dim = [1.450 0.6];
surface{13}.speed = [0; 2.25; 0];
surface{13}.transform = surface{10}.transform * [eye(3) [0; 2.550; 0]; zeros(1,3), 1];
surface{13}.id = "SB03";

%Shortbelt 4
surface{14}.dim = [1.450 0.6];
surface{14}.speed = [0; 2.5; 0];
surface{14}.transform = surface{10}.transform * [eye(3) [0; 3.150; 0]; zeros(1,3), 1];
surface{14}.id = "SB04";

%Shortbelt 5
surface{15}.dim = [1.450 0.6];
surface{15}.speed = [0; 2.75; 0];
surface{15}.transform = surface{10}.transform * [eye(3) [0; 3.750; 0]; zeros(1,3), 1];
surface{15}.id = "SB05";

%MultiBelt Junction 60deg
nbelts = 14; %number of multibelts
Length = 2.8; %longest belt length
width = 1.45; %total multibelt width
deg = 30;      %multibelt angle
vel = 2.89; %Based on crossorter speed of 2.5  
transform = surface{10}.transform * [eye(3) [0; 4.050; 0]; zeros(1,3), 1];

%Create the multibelt
mbelt = multibelt(nbelts,Length,width,deg,vel,transform,2); 

for ii = 1:length(mbelt)
    surface{15+ii} = mbelt{ii};
end

%Crossorter deck SCS700 
% surface{15+ii+1}.dim =[0.700 0.485];
% surface{15+ii+1}.speed =[1.44; 2.5; 0];
% surface{15+ii+1}.transform= [Rz(90), [5; 6.1; 0]; zeros(1,3),1];
% surface{15+ii+1}.id = "CS01";
% pitch = 0.525;

%Crossorter deck SCS1200
surface{15+ii+1}.dim =[1.200 0.520];
surface{15+ii+1}.speed =[1.44; 2.5; 0];
surface{15+ii+1}.transform= [Rz(90), [5; 6.35; 0]; zeros(1,3),1];
surface{15+ii+1}.id = "CS01";
pitch = 0.7;

nconv = length(surface); %check how many conveyors we already have
for jj = 1:10 %create 10 additional decks
    surface{nconv+jj}.dim = surface{nconv}.dim; 
    surface{nconv+jj}.speed = surface{nconv}.speed;
    surface{nconv+jj}.transform= surface{nconv}.transform*[eye(3), [0; jj*-pitch; 0]; zeros(1,3),1];
    surface{nconv+jj}.id = sprintf('CS%0.2d',jj+1);
end

% plotEnvironment(surface,c.dt,ii);

%% Run the simulation
tic
[AH_B, BV_AB, ~, ~] = BoxSimulatorLCP(x,c,box,surface);
toc

% tic
% [AH_B_fp,BV_AB_fp,PN,PT] = BoxSimulator(x,c,box,surface);
% toc

%% Figures
%Set plots to use LaTeX interface
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%Extract some values
dt =c.dt;
runtime=c.endtime;

%Plot the trajectory of the box
if doPlot

    %Plot the simulation results into a figure
    figure(Position=[232,246,1560,820]);
    pause;
    for ii=1:(1/dt)/(runtime*4):(runtime/dt+1)
        try; plotBox(AH_B_fp(:,:,ii),box,[0 0 1]); catch; end; 
        try; plotBox(AH_B(:,:,ii),box,[1 0 0]); catch; end; hold on;
        plotEnvironment(surface,dt,ii);
    end
end
if dosave ==1
    save('AH_B.mat','AH_B');
end
%%
if MakeVideo
    close all;
    video = VideoWriter('static/Boxsimulator.avi'); %create the video object
    video.FrameRate=0.5/c.dt/step;
    open(video); %open the file for writing

    figure(1);
    for ii=1:step:length(AH_B)
        plotBox(AH_B(:,:,ii),box,[0 0 1]);        
        
        %Plot the inclined table C
        for jj=1:length(surface)
            table3 = fill3(spoints{jj}(1,1:4),spoints{jj}(2,1:4),spoints{jj}(3,1:4),1);hold on;
            set(table3,'FaceColor',0.8*[1 1 1],'FaceAlpha',1);


            %Plot the origin of the contact surface with its unit vectors
            % tip = [surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,1) surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,2) surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,3)];
            % plot3([surface{jj}.transform(1,4) tip(1,1)],[surface{jj}.transform(2,4) tip(2,1)],[surface{jj}.transform(3,4) tip(3,1)],'r'); hold on
            % plot3([surface{jj}.transform(1,4) tip(1,2)],[surface{jj}.transform(2,4) tip(2,2)],[surface{jj}.transform(3,4) tip(3,2)],'g');
            % plot3([surface{jj}.transform(1,4) tip(1,3)],[surface{jj}.transform(2,4) tip(2,3)],[surface{jj}.transform(3,4) tip(3,3)],'b');


            % %Draw the velocity of the contact plane
            % temp = (vecveltemp+surface{jj}.speed*(dt*(ii-1))); %Move the grid according to the conveyor speed
            % pbool = temp(1,:)>(-0.5*surface{jj}.dim(1))&temp(1,:)<(0.5*surface{jj}.dim(1))&temp(2,:)>(-0.5*surface{jj}.dim(2))&temp(2,:)<(0.5*surface{jj}.dim(2)); %Select the grid points inside the surface area            
            % vecvel = surface{jj}.transform(1:3,1:3)*temp+surface{jj}.transform(1:3,4); %Rotate and translate those points according to surface pose
            % speed = surface{jj}.speed/norm(surface{jj}.speed); %Get the normalized velocity vector
            % vecvel2 = surface{jj}.transform(1:3,1:3)*repmat(0.15*speed,1,length(vecvel)); %Get the end points of the velocity vector
            % 
            % quiver3(vecvel(1,pbool),vecvel(2,pbool),vecvel(3,pbool),vecvel2(1,pbool),vecvel2(2,pbool),vecvel2(3,pbool),'off','color',[0 0.4470 0.7410]);
        end
        
        %Plot the origin of the world coordinate frame
        tip = [0.3*[1;0;0] 0.3*[0;1;0] 0.3*[0;0;1]];
        plot3([0 tip(1,1)],[0 tip(2,1)],[0 tip(3,1)],'r'); hold on
        plot3([0 tip(1,2)],[0 tip(2,2)],[0 tip(3,2)],'g');
        plot3([0 tip(1,3)],[0 tip(2,3)],[0 tip(3,3)],'b');

        grid on;axis equal; 
        axis([-1.5 1.5 -0.7 3.5 -0.3 0.7]);
        axis([-0.7 0.7 -0.7 2 -0.3 0.7]);
        axis([-1 1 -0.7 2 -0.3 0.7]);
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        view(-35,31);
%         axis off
        drawnow
        hold off

        fname = ['static/image']; % full name of image
        print('-djpeg','-r500',fname)     % save image with '-r200' resolution
        I = imread([fname '.jpg']);       % read saved image
        frame = im2frame(I);              % convert image to frame
        writeVideo(video,frame); %write the image to file

        hold off;
    end
    close(video); %close the file
end

%% Plot Position, Velocity, Force
% AH_Bm = cat(3,AH_B{:});
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
% close all;
time=0:dt:runtime;

figure('pos',[95,530,560,420]); 
plot(time(1:end-1),flipud(PN)');
% axis([0 10 0.002451 0.002454]);
grid on;
legend('Contact point 1','Contact point 2','Contact point 3','Contact point 4','Contact point 5','Contact point 6','Contact point 7','Contact point 8');
ylabel('$P_N$ [N]');
xlabel('time [s]');
title('Normal force of the contact points for a box of 1 kg');

figure('pos',[670,530,560,420]); 
sgtitle('Velocity of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(time(1:end-1),BV_AB_fp(1,1:end-1)); hold on;
try; plot(time(1:end-1),BV_AB(1,1:end-1)); catch; end;
ylabel('Velocity [m/s]');
xlabel('Time [s]');
title('x');
legend('FP','LCP')

subplot(1,3,2);
plot(time(1:end-1),BV_AB_fp(2,1:end-1)); hold on;
try; plot(time(1:end-1),BV_AB(2,1:end-1)); catch; end;
xlabel('Time [s]');
title('y');

subplot(1,3,3);
plot(time(1:end-1),BV_AB_fp(3,1:end-1)); hold on;
try; plot(time(1:end-1),BV_AB(3,1:end-1)); catch; end;
xlabel('Time [s]');
title('z');


figure('pos',[1245,530,560,420]); 
sgtitle('Position of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(time(1:end-1),squeeze(AH_B_fp(1,4,1:end-1))); hold on;
try; plot(time(1:end-1),squeeze(AH_B(1,4,1:end-1))); catch; end;
title('x');
xlabel('Time [s]');
ylabel('Position [m]');
legend('FP','LCP')

subplot(1,3,2);
plot(time(1:end-1),squeeze(AH_B_fp(2,4,1:end-1))); hold on;
try; plot(time(1:end-1),squeeze(AH_B(2,4,1:end-1))); catch; end;
title('y');
xlabel('Time [s]');

subplot(1,3,3);
plot(time(1:end-1),squeeze(AH_B_fp(3,4,1:end-1))); hold on;
try; plot(time(1:end-1),squeeze(AH_B(3,4,1:end-1))); catch; end;
title('z');
xlabel('Time [s]');


%% Functions
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


function Bplot = plotBox(AH_B,box,color) 
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
%         tip = [q(:,1)+ 0.3*R1(:,1) q(:,1)+ 0.3*R2(:,1) q(:,1)+ 0.3*R3(:,1)];
%         plot3([q(1,1) tip(1,1)],[q(2,1) tip(2,1)],[q(3,1) tip(3,1)],'r'); hold on
%         plot3([q(1,1) tip(1,2)],[q(2,1) tip(2,2)],[q(3,1) tip(3,2)],'g');
%         plot3([q(1,1) tip(1,3)],[q(2,1) tip(2,3)],[q(3,1) tip(3,3)],'b');
        
        %Create the box
        pbool = (abs(box.vertices(1,:))==max(abs(box.vertices(1,:))))&(abs(box.vertices(2,:))==max(abs(box.vertices(2,:))))&(abs(box.vertices(3,:))==max(abs(box.vertices(3,:))));
        Ap = AR_B*box.vertices(:,pbool)+AH_B(1:3,4);
%         Ap = AR_B*box.vertices+AH_B(1:3,4);
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
        
        %Color the surfaces of the box
        Bplot = fill3([Ap_1(1) Ap_2(1) Ap_6(1) Ap_5(1)],[Ap_1(2) Ap_2(2) Ap_6(2) Ap_5(2)],[Ap_1(3) Ap_2(3) Ap_6(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',1);%Face F
        fill3([Ap_1(1) Ap_2(1) Ap_3(1) Ap_4(1)],[Ap_1(2) Ap_2(2) Ap_3(2) Ap_4(2)],[Ap_1(3) Ap_2(3) Ap_3(3) Ap_4(3)],1,'FaceColor',color,'FaceAlpha',1);%Face A
        fill3([Ap_8(1) Ap_7(1) Ap_6(1) Ap_5(1)],[Ap_8(2) Ap_7(2) Ap_6(2) Ap_5(2)],[Ap_8(3) Ap_7(3) Ap_6(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',1);%Face C
        fill3([Ap_8(1) Ap_7(1) Ap_3(1) Ap_4(1)],[Ap_8(2) Ap_7(2) Ap_3(2) Ap_4(2)],[Ap_8(3) Ap_7(3) Ap_3(3) Ap_4(3)],1,'FaceColor',color,'FaceAlpha',1);%Face E
        fill3([Ap_1(1) Ap_4(1) Ap_8(1) Ap_5(1)],[Ap_1(2) Ap_4(2) Ap_8(2) Ap_5(2)],[Ap_1(3) Ap_4(3) Ap_8(3) Ap_5(3)],1,'FaceColor',color,'FaceAlpha',1);%Face D
        fill3([Ap_2(1) Ap_3(1) Ap_7(1) Ap_6(1)],[Ap_2(2) Ap_3(2) Ap_7(2) Ap_6(2)],[Ap_2(3) Ap_3(3) Ap_7(3) Ap_6(3)],1,'FaceColor',color,'FaceAlpha',1);%Face B

        %Plot all contact points
%         vertices = AR_B*box.vertices+AH_B(1:3,4);
%         plot3(vertices(1,:),vertices(2,:),vertices(3,:),'.',MarkerSize=1)
        
end

function plotEnvironment(surface,dt,ii)

%Crossorter movement
for jj = 1:11
    surface{29+jj}.transform = surface{29+jj}.transform*expm(hat([0;2.5;0;0;0;0]*dt*ii));
end

% Plotting options
%For plotting the contact surface
for jj=1:length(surface)
    ws    = surface{jj}.dim(1);                 %With of the contact surface             [m]
    ls    = surface{jj}.dim(2);               %Length of the contact surface           [m]
    surfacepoints = [0.5*ws -0.5*ws -0.5*ws 0.5*ws 0.5*ws; -0.5*ls -0.5*ls 0.5*ls 0.5*ls -0.5*ls; 0 0 0 0 0;];
    spoints{jj} = surface{jj}.transform(1:3,1:3)*surfacepoints + surface{jj}.transform(1:3,4); %Transform the vertices according to position/orientation of the surface
end

%Vector field
X = linspace(-20,20,205);
Y = linspace(-20,20,205);
Z = 0;
tel = 1;
vecveltemp =[];
for xx =1:length(X)
    for yy = 1:length(Y)
        vecveltemp(:,tel) = [X(xx);Y(yy);Z]';
        tel=tel+1;
    end
end

% Plot the origin of the world coordinate frame
tip = [0.3*[1;0;0] 0.3*[0;1;0] 0.3*[0;0;1]];
plot3([0 tip(1,1)],[0 tip(2,1)],[0 tip(3,1)],'r'); hold on
plot3([0 tip(1,2)],[0 tip(2,2)],[0 tip(3,2)],'g');
plot3([0 tip(1,3)],[0 tip(2,3)],[0 tip(3,3)],'b');

%Plot the inclined table C
for jj=1:length(surface)
    table3 = fill3(spoints{jj}(1,1:4),spoints{jj}(2,1:4),spoints{jj}(3,1:4),1);hold on;
    set(table3,'FaceColor',0.8*[1 1 1],'FaceAlpha',1);

    % %Plot the origin of the contact surface with its unit vectors
    % tip = [surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,1) surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,2) surface{jj}.transform(1:3,4)+0.3*surface{jj}.transform(1:3,3)];
    % plot3([surface{jj}.transform(1,4) tip(1,1)],[surface{jj}.transform(2,4) tip(2,1)],[surface{jj}.transform(3,4) tip(3,1)],'r'); hold on
    % plot3([surface{jj}.transform(1,4) tip(1,2)],[surface{jj}.transform(2,4) tip(2,2)],[surface{jj}.transform(3,4) tip(3,2)],'g');
    % plot3([surface{jj}.transform(1,4) tip(1,3)],[surface{jj}.transform(2,4) tip(2,3)],[surface{jj}.transform(3,4) tip(3,3)],'b');

    % % Draw the velocity of the contact plane
    % temp = (vecveltemp+surface{jj}.speed*(dt*(ii-1))); %Move the grid according to the conveyor speed
    % pbool = temp(1,:)>(-0.5*surface{jj}.dim(1))&temp(1,:)<(0.5*surface{jj}.dim(1))&temp(2,:)>(-0.5*surface{jj}.dim(2))&temp(2,:)<(0.5*surface{jj}.dim(2)); %Select the grid points inside the surface area
    % vecvel = surface{jj}.transform(1:3,1:3)*temp+surface{jj}.transform(1:3,4); %Rotate and translate those points according to surface pose
    % speed = surface{jj}.speed/norm(surface{jj}.speed); %Get the normalized velocity vector
    % vecvel2 = surface{jj}.transform(1:3,1:3)*repmat(0.15*speed,1,length(vecvel)); %Get the end points of the velocity vector
    % 
    % quiver3(vecvel(1,pbool),vecvel(2,pbool),vecvel(3,pbool),vecvel2(1,pbool),vecvel2(2,pbool),vecvel2(3,pbool),'off','color',[0 0.4470 0.7410]);
end

grid on;axis equal;
% axis([-1 1 -0.7 2 -0.3 0.7]);
axis([-8 2 -2 7 -0.1 0.7]);
% xlabel('x [m]');
% ylabel('y [m]');
% zlabel('z [m]');
view(-35,31);
ax = gca;
ax.Clipping = "off";
axis off;
hold off
drawnow
end


function mbelt = multibelt(nbelts,length,width,deg,vel,transform,id)

%MultiBelt Junction 60deg
% nbelts = 8; %number of multibelts
% length = 1.050; %longest belt length
% width = 0.760; %total multibelt width
% deg = 60;      %multibelt angle
% vel = 2.5;     %multibelt velocity in m/s

x = width/nbelts; %multibelts width
slength = length -(width-x)*tan(deg2rad(90-deg)); %Shortest multibelt length
y = linspace(length,slength,nbelts); %multibelt lengths


for ii = 1:nbelts
    mbelt{ii}.dim = [x y(ii)];
    mbelt{ii}.speed = [0; vel; 0];
    mbelt{ii}.transform = [eye(3) [(x/2)+(ii-1)*x; (ii-1)*(mean(diff(y))/2); 0]; zeros(1,3), 1];
    mbelt{ii}.id = "MB" + num2str(id) + "_"+num2str(ii);

    %transform origin of multibelt to middle in x, and zero in y
    mbelt{ii}.transform = mbelt{ii}.transform*[eye(3), [-0.5*width; 0.5*length; 0]; zeros(1,3),1];

    %transform multibelt to desired pose
    mbelt{ii}.transform = transform*mbelt{ii}.transform;
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
%% Hat operator
if length(vec) == 3
    res = [0, -vec(3), vec(2); vec(3), 0, -vec(1); -vec(2), vec(1), 0];
elseif length(vec) == 6
    skew = [0, -vec(6), vec(5); vec(6), 0, -vec(4); -vec(5), vec(4), 0];
    v = [vec(1);vec(2);vec(3)];
    res = [skew, v; zeros(1,4)];
end
end