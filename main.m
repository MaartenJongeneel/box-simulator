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
% SingleConveyor;
DoubleConveyor;
% IFZ41;
%% Parameters for input
c.a                  = 0.001;           %Prox point auxilary parameter             [-]
c.tol                = 1e-7;            %Error tol for fixed-point                 [-]
c.m                  = 1;               %Mass of the box                           [kg]  
c.endtime            = 5;             %Runtime of the simulation                 [s]
c.dt                 = 1/100;          %Timestep at which the simulator runs      [s]
c.dimd               = 4;              %Discretization of the friction cone
step                 = ceil(1/c.dt/50); %Number of discrete points we skip per shown frame
%% Read the scene data
x.releaseOrientation = box.release.orientation;  %Release orientation of the box            [deg]
x.releasePosition    = box.release.position';    %Release position of the box               [m]
x.releaseLinVel      = box.release.linVel';      %Release linear velocity (expressed in B)  [m/s]
x.releaseAngVel      = box.release.angVel';      %Release angular velocity (expressed in B) [rad/s]
c.eN                 = box.parameters.eN;        %Normal coefficient of restitution         [-]
c.eT                 = box.parameters.eT;        %Tangential coefficient of restitution     [-]s
c.mu                 = box.parameters.mu;        %Coefficient of friction                   [-]
box.B_M_B            = box.inertia_tensor;       %Rewrite inertia tensor
%% Create the box struct
%Discretization of the box vertices
Ndisc=box.discretization;
[X,Y,Z]=meshgrid(linspace(-box.dimensions(1)/2,box.dimensions(1)/2,Ndisc),linspace(-box.dimensions(2)/2,box.dimensions(2)/2,Ndisc),linspace(-box.dimensions(3)/2,box.dimensions(3)/2,Ndisc));
pbool = (abs(X(:))==box.dimensions(1)/2) | (abs(Y(:))==box.dimensions(2)/2) | (abs(Z(:))==box.dimensions(3)/2);
box.vertices= [X(pbool)';Y(pbool)';Z(pbool)'];

%% Run the simulation
tic
[AH_B, BV_AB, PN, PT] = BoxSimulatorLCP(x,c,box,surface);
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
try; plot(time(1:end-1),BV_AB_fp(1,1:end-1)); hold on; catch; end;
try; plot(time(1:end-1),BV_AB(1,1:end-1)); catch; end;
ylabel('Velocity [m/s]');
xlabel('Time [s]');
title('x');
legend('FP','LCP')

subplot(1,3,2);
try; plot(time(1:end-1),BV_AB_fp(2,1:end-1)); hold on; catch; end;
try; plot(time(1:end-1),BV_AB(2,1:end-1)); catch; end;
xlabel('Time [s]');
title('y');

subplot(1,3,3);
try;  plot(time(1:end-1),BV_AB_fp(3,1:end-1)); hold on; catch; end;
try; plot(time(1:end-1),BV_AB(3,1:end-1)); catch; end;
xlabel('Time [s]');
title('z');


figure('pos',[1245,530,560,420]); 
sgtitle('Position of the COM of the box in x,y,z-direction');
subplot(1,3,1);
try; plot(time(1:end-1),squeeze(AH_B_fp(1,4,1:end-1))); hold on; catch; end;
try; plot(time(1:end-1),squeeze(AH_B(1,4,1:end-1))); catch; end;
title('x');
xlabel('Time [s]');
ylabel('Position [m]');
legend('FP','LCP')

subplot(1,3,2);
try; plot(time(1:end-1),squeeze(AH_B_fp(2,4,1:end-1))); hold on; catch; end;
try; plot(time(1:end-1),squeeze(AH_B(2,4,1:end-1))); catch; end;
title('y');
xlabel('Time [s]');

subplot(1,3,3);
try; plot(time(1:end-1),squeeze(AH_B_fp(3,4,1:end-1))); hold on; catch; end;
try; plot(time(1:end-1),squeeze(AH_B(3,4,1:end-1))); catch; end;
title('z');
xlabel('Time [s]');
