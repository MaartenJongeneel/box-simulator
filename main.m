clearvars; close all; clc; warning ('off','all'); addpath('Functions');
%This script is used to simulate a box that is tossed on a surface. In the
%settings below, one can set different settings, such as the size of the
%box or the coefficients of friction, tangential and normal restitution.
%% General settings
dosave             = false;         %Save the trajectory (AH_B) to a .mat file
doPlot             = true;          %Show the trajectory of the box
%% Parameters for input
releaseOrientation = Rx(0);         %Release orientation of the box            [deg]
releasePosition    = [0; 0; 0.0251];%Release position of the box               [m]
releaseLinVel      = [0; 0; 0];     %Release linear velocity (expressed in B)  [m/s]
releaseAngVel      = [0; 0; 0];     %Release angular velocity (expressed in B) [rad/s]
eN                 = 0.4;           %Normal coefficient of restitution         [-]
eT                 = 0.4;           %Tangential coefficient of restitution     [-]
mu                 = 0.4;           %Coefficient of friction                   [-]
l                  = 0.1;           %Length of the box                         [m]
w                  = 0.15;          %Width of the box                          [m]
h                  = 0.05;          %Height of the box                         [m]
m                  = 1;             %Mass of the box                           [kg]  
AR_C               = Rx(0);         %Orientation of the contact plane          [deg]
Ao_C               = [0; 0.6; 0];   %Position of the contact plane             [m]
runtime            = 2;             %Runtime of the simulation                 [s]
dt                 = 1/1000;        %Timestep at which the simulator runs      [s]

%% Create the box struct
%Mass matrix of the box
Ml = m*eye(3);

%Inertia matrix of the box
I  = [(m/12)*(w^2+h^2),                 0,                  0;
                     0,  (m/12)*(l^2+h^2),                  0;
                     0,                 0,   (m/12)*(l^2+w^2);];
%Inertia tensor
box.B_M_B = [Ml zeros(3,3); zeros(3,3) I];

%Vertices
box.vertices = BoxVertices(l,w,h); 

%Mass
box.mass = m;

%% Run the dynamics
[AH_B,BV_AB,FN,FT] = BoxSimulator(releasePosition,releaseOrientation,releaseLinVel,releaseAngVel,eN,eT,mu,box,AR_C,Ao_C,dt,runtime);

restPosition = AH_B{end}(1:3,4);      %Rest position of the box
restOrientation = AH_B{end}(1:3,1:3); %Rest orientation of the box

%% Figures
%Set plots to use LaTeX interface
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Plotting options 
%For plotting the contact surface
ws    = 1;                 %With of the contact surface             [m]
ls    = 1.6;               %Length of the contact surface           [m]
surfacepoints = [0.5*ws -0.5*ws -0.5*ws 0.5*ws 0.5*ws; -0.5*ls -0.5*ls 0.5*ls 0.5*ls -0.5*ls; 0 0 0 0 0;];
spoints = AR_C*surfacepoints +Ao_C; %Transform the vertices according to position/orientation of the surface

%Plot the trajectory of the box
if doPlot
    figure;
    for ii=1:15:length(AH_B)
        plotBox(AH_B{ii},box,[0 0 1]);
        
        
        %Plot the inclined table C
        table3 = fill3(spoints(1,1:4),spoints(2,1:4),spoints(3,1:4),1);hold on;
        set(table3,'FaceColor',0.8*[1 1 1],'FaceAlpha',1);
        
        %Plot the origin of the contact surface with its unit vectors
        tip = [Ao_C+0.3*AR_C(:,1) Ao_C+0.3*AR_C(:,2) Ao_C+0.3*AR_C(:,3)];
        plot3([Ao_C(1) tip(1,1)],[Ao_C(2) tip(2,1)],[Ao_C(3) tip(3,1)],'r'); hold on
        plot3([Ao_C(1) tip(1,2)],[Ao_C(2) tip(2,2)],[Ao_C(3) tip(3,2)],'g');
        plot3([Ao_C(1) tip(1,3)],[Ao_C(2) tip(2,3)],[Ao_C(3) tip(3,3)],'b');
        
        %Plot the origin of the world coordinate frame
        tip = [0.3*[1;0;0] 0.3*[0;1;0] 0.3*[0;0;1]];
        plot3([0 tip(1,1)],[0 tip(2,1)],[0 tip(3,1)],'r'); hold on
        plot3([0 tip(1,2)],[0 tip(2,2)],[0 tip(3,2)],'g');
        plot3([0 tip(1,3)],[0 tip(2,3)],[0 tip(3,3)],'b');
        
        
        grid on;axis equal;
        axis([-0.5 0.5 -0.5 2 -0.3 0.6]);
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        view(-35,31);
        drawnow
        hold off
    end
end
if dosave ==1
    save('AH_B.mat','AH_B');
end


% Plot Position, Velocity, Force
AH_Bm = cat(3,AH_B{:});
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
% close all;

figure; 
plot(0:dt:runtime,flipud(FN)');
% axis([0 10 0.002451 0.002454]);
grid on;
legend('Contact point 1','Contact point 2','Contact point 3','Contact point 4');
ylabel('$F_N$ [N]');
xlabel('time [s]');
title('Normal force of the contact points for a box of 1 kg');

figure;
sgtitle('Velocity of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(0:dt:runtime,BV_AB(1,1:end-1));
ylabel('Velocity [m/s]');
xlabel('Time [s]');
title('x');

subplot(1,3,2);
plot(0:dt:runtime,BV_AB(2,1:end-1));
xlabel('Time [s]');
title('y');

subplot(1,3,3);
plot(0:dt:runtime,BV_AB(3,1:end-1));
xlabel('Time [s]');
title('z');
ylim([-8e-6 8e-6]);

figure;
sgtitle('Position of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(0:dt:runtime,squeeze(AH_Bm(1,4,1:end-1)));
title('x');
xlabel('Time [s]');
ylabel('Position [m]');

subplot(1,3,2);
plot(0:dt:runtime,squeeze(AH_Bm(2,4,1:end-1)))
title('y');
xlabel('Time [s]');

subplot(1,3,3);
plot(0:dt:runtime,squeeze(AH_Bm(3,4,1:end-1)))
title('z');
xlabel('Time [s]');