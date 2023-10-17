clearvars; close all; clc; warning ('off','all'); addpath('Functions');
%This script is used to simulate a box that is tossed on a surface. In the
%settings below, one can set different settings, such as the size of the
%box or the coefficients of friction, tangential and normal restitution.
%% General settings
dosave             = false;         %Save the trajectory (AH_B) to a .mat file
doPlot             = true;          %Show the trajectory of the box
%% Parameters for input
x.releaseOrientation = Rx(0);         %Release orientation of the box            [deg]
x.releasePosition    = [0; 0; 0.3];   %Release position of the box               [m]
x.releaseLinVel      = [0; 0; 0];     %Release linear velocity (expressed in B)  [m/s]
x.releaseAngVel      = [0; 0; 0];     %Release angular velocity (expressed in B) [rad/s]
c.eN                 = 0.4;           %Normal coefficient of restitution         [-]
c.eT                 = 0.0;           %Tangential coefficient of restitution     [-]
c.mu                 = 0.6;           %Coefficient of friction                   [-]
l                    = 0.1;           %Length of the box                         [m]
w                    = 0.15;          %Width of the box                          [m]
h                    = 0.05;          %Height of the box                         [m]
c.a                  = 0.01;          %Prox point auxilary parameter             [-]
c.tol                = 1e-9;          %Error tol for fixed-point                 [-]
c.m                  = 1;             %Mass of the box                           [kg]  
c.AR_C               = Rx(0);         %Orientation of the contact plane          [deg]
c.Ao_C               = [0; 0.6; 0];   %Position of the contact plane             [m]
c.Cv_AC              = [0; 0; 0];     %Linear velocity of the contact plane      [m/s]
c.endtime            = 1.5;           %Runtime of the simulation                 [s]
c.dt                 = 1/1000;        %Timestep at which the simulator runs      [s]

%% Create the box struct
%Mass matrix of the box
Ml = c.m*eye(3);

%Inertia matrix of the box
I  = [(c.m/12)*(w^2+h^2),                 0,                  0;
                     0,  (c.m/12)*(l^2+h^2),                  0;
                     0,                 0,   (c.m/12)*(l^2+w^2);];
%Inertia tensor
box.B_M_B = [Ml zeros(3,3); zeros(3,3) I];

%Vertices of the box expressed in body fixed frame B (frame B at COM)
box.vertices(:,1) = [ l/2; -w/2; -h/2];      box.vertices(:,5) = [ l/2; -w/2;  h/2];
box.vertices(:,2) = [-l/2; -w/2; -h/2];      box.vertices(:,6) = [-l/2; -w/2;  h/2];
box.vertices(:,3) = [-l/2;  w/2; -h/2];      box.vertices(:,7) = [-l/2;  w/2;  h/2];
box.vertices(:,4) = [ l/2;  w/2; -h/2];      box.vertices(:,8) = [ l/2;  w/2;  h/2];

%Mass
box.mass = c.m;

%% Run the dynamics
[AH_B,BV_AB,FN,FT] = BoxSimulator(x,c,box);

restPosition = AH_B{end}(1:3,4);      %Rest position of the box
restOrientation = AH_B{end}(1:3,1:3); %Rest orientation of the box

%% Figures
%Set plots to use LaTeX interface
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%Extract some values
AR_C=c.AR_C;
Ao_C=c.Ao_C;
dt =c.dt;
runtime=c.endtime;


% Plotting options 
%For plotting the contact surface
ws    = 1;                 %With of the contact surface             [m]
ls    = 1.6;               %Length of the contact surface           [m]
surfacepoints = [0.5*ws -0.5*ws -0.5*ws 0.5*ws 0.5*ws; -0.5*ls -0.5*ls 0.5*ls 0.5*ls -0.5*ls; 0 0 0 0 0;];
spoints = AR_C*surfacepoints +Ao_C; %Transform the vertices according to position/orientation of the surface
step = 25; %Number of discrete points we skip per shown frame

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


%Plot the trajectory of the box
if doPlot
    figure;
    for ii=1:step:length(AH_B)
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

        %Draw the velocity vector of the
%         plot3(vecvel(1,:),vecvel(2,:),vecvel(3,:),'o'); hold on;
%         annotation('textarrow',[vecvel(1,:)' vecvel2(1,:)'],[vecvel(2,:)' vecvel2(2,:)'],[vecvel(3,:)' vecvel2(3,:)'])

        vecvel = AR_C*(vecveltemp+c.Cv_AC*(dt*(ii-1)))+Ao_C;
        speed = c.Cv_AC/norm(c.Cv_AC);
        vecvel2 = AR_C*repmat(0.15*speed,1,length(vecvel));
        pbool = (vecvel(1,:)>=(min(spoints(1,:))))&(vecvel(1,:)<=(max(spoints(1,:))))&(vecvel(2,:)>=(min(spoints(2,:))))&(vecvel(2,:)<=(max(spoints(2,:))))&(vecvel(3,:)>=(min(spoints(3,:))))&(vecvel(3,:)<=(max(spoints(3,:))));
        quiver3(vecvel(1,pbool),vecvel(2,pbool),vecvel(3,pbool),vecvel2(1,pbool),vecvel2(2,pbool),vecvel2(3,pbool),'off');

        
        
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


% % Plot Position, Velocity, Force
% AH_Bm = cat(3,AH_B{:});
% set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
% % close all;
% time=0:dt:runtime;
% 
% figure; 
% plot(time(1:end-1),flipud(FN)');
% % axis([0 10 0.002451 0.002454]);
% grid on;
% legend('Contact point 1','Contact point 2','Contact point 3','Contact point 4');
% ylabel('$F_N$ [N]');
% xlabel('time [s]');
% title('Normal force of the contact points for a box of 1 kg');
% 
% figure;
% sgtitle('Velocity of the COM of the box in x,y,z-direction');
% subplot(1,3,1);
% plot(time(1:end-1),BV_AB(1,1:end-1));
% ylabel('Velocity [m/s]');
% xlabel('Time [s]');
% title('x');
% 
% subplot(1,3,2);
% plot(time(1:end-1),BV_AB(2,1:end-1));
% xlabel('Time [s]');
% title('y');
% 
% subplot(1,3,3);
% plot(time(1:end-1),BV_AB(3,1:end-1));
% xlabel('Time [s]');
% title('z');
% ylim([-8e-6 8e-6]);
% 
% figure;
% sgtitle('Position of the COM of the box in x,y,z-direction');
% subplot(1,3,1);
% plot(time(1:end-1),squeeze(AH_Bm(1,4,1:end-1)));
% title('x');
% xlabel('Time [s]');
% ylabel('Position [m]');
% 
% subplot(1,3,2);
% plot(time(1:end-1),squeeze(AH_Bm(2,4,1:end-1)))
% title('y');
% xlabel('Time [s]');
% 
% subplot(1,3,3);
% plot(time(1:end-1),squeeze(AH_Bm(3,4,1:end-1)))
% title('z');
% xlabel('Time [s]');