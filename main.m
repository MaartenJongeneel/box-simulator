clearvars; close all; clc; warning ('off','all'); addpath('Functions');
%This script is used to simulate a box that is tossed on a surface. In the
%settings below, one can set different settings, such as the size of the
%box or the coefficients of friction, tangential and normal restitution.
%% General settings
dosave             = false;         %Save the trajectory (AH_B) to a .mat file
doPlot             = true;          %Show the trajectory of the box
MakeVideo          = false;         %Save the simulation result to video
%% Parameters for input
x.releaseOrientation = Rx(0);         %Release orientation of the box            [deg]
x.releasePosition    = [0; 0; 0.3];   %Release position of the box               [m]
x.releaseLinVel      = [0; 0; 0];     %Release linear velocity (expressed in B)  [m/s]
x.releaseAngVel      = [3; 1; 0];     %Release angular velocity (expressed in B) [rad/s]
c.eN                 = 0.4;           %Normal coefficient of restitution         [-]
c.eT                 = 0.0;           %Tangential coefficient of restitution     [-]
c.mu                 = 0.6;           %Coefficient of friction                   [-]
l                    = 0.1;           %Length of the box                         [m]
w                    = 0.15;          %Width of the box                          [m]
h                    = 0.05;          %Height of the box                         [m]
c.a                  = 0.001;         %Prox point auxilary parameter             [-]
c.tol                = 1e-7;          %Error tol for fixed-point                 [-]
c.m                  = 1;             %Mass of the box                           [kg]  
c.AR_C               = Rx(10)*Ry(5);  %Orientation of the contact plane          [deg]
c.Ao_C               = [0; 0.6; 0];   %Position of the contact plane             [m]
c.Cv_AC              = [0; 1; 0];     %Linear velocity of the contact plane      [m/s]
c.endtime            = 1.5;           %Runtime of the simulation                 [s]
c.dt                 = 1/1000;        %Timestep at which the simulator runs      [s]
step                 = 1/c.dt/100;    %Number of discrete points we skip per shown frame

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

        %Draw the velocity of the contact plane
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

%%
if MakeVideo
    video = VideoWriter('static/Boxsimulator.avi'); %create the video object
    video.FrameRate=0.5/c.dt/step;
    open(video); %open the file for writing

    figure(1);
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

        %Draw the velocity of the contact plane
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
AH_Bm = cat(3,AH_B{:});
set(groot,'defaulttextinterpreter','latex'); set(groot,'defaultAxesTickLabelInterpreter','latex'); set(groot,'defaultLegendInterpreter','latex');
% close all;
time=0:dt:runtime;

figure; 
plot(time(1:end-1),flipud(FN)');
% axis([0 10 0.002451 0.002454]);
grid on;
legend('Contact point 1','Contact point 2','Contact point 3','Contact point 4');
ylabel('$F_N$ [N]');
xlabel('time [s]');
title('Normal force of the contact points for a box of 1 kg');

figure;
sgtitle('Velocity of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(time(1:end-1),BV_AB(1,1:end-1));
ylabel('Velocity [m/s]');
xlabel('Time [s]');
title('x');

subplot(1,3,2);
plot(time(1:end-1),BV_AB(2,1:end-1));
xlabel('Time [s]');
title('y');

subplot(1,3,3);
plot(time(1:end-1),BV_AB(3,1:end-1));
xlabel('Time [s]');
title('z');

figure;
sgtitle('Position of the COM of the box in x,y,z-direction');
subplot(1,3,1);
plot(time(1:end-1),squeeze(AH_Bm(1,4,1:end-1)));
title('x');
xlabel('Time [s]');
ylabel('Position [m]');

subplot(1,3,2);
plot(time(1:end-1),squeeze(AH_Bm(2,4,1:end-1)))
title('y');
xlabel('Time [s]');

subplot(1,3,3);
plot(time(1:end-1),squeeze(AH_Bm(3,4,1:end-1)))
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
%% Box-simulator-FixedPoint:
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
%% Plot the box
        AR_B = AH_B(1:3,1:3);
        %Output the position of the current time step for plotting purposes
        ii =1;
        q(:,ii)  = AH_B(1:3,4);
        Ao_B1{ii} = AH_B(1:3,4);
        R1(:,ii) = AH_B(1:3,1);
        R2(:,ii) = AH_B(1:3,2);
        R3(:,ii) = AH_B(1:3,3);
        
        %Plot the origin of the box with its unit vectors
        tip = [q(:,ii)+ 0.3*R1(:,ii) q(:,ii)+ 0.3*R2(:,ii) q(:,ii)+ 0.3*R3(:,ii)];
        plot3([q(1,ii) tip(1,1)],[q(2,ii) tip(2,1)],[q(3,ii) tip(3,1)],'r'); hold on
        plot3([q(1,ii) tip(1,2)],[q(2,ii) tip(2,2)],[q(3,ii) tip(3,2)],'g');
        plot3([q(1,ii) tip(1,3)],[q(2,ii) tip(2,3)],[q(3,ii) tip(3,3)],'b');
        
        %Create the box
        Ap_1 = AR_B*box.vertices(:,1) + Ao_B1{ii};
        Ap_2 = AR_B*box.vertices(:,2) + Ao_B1{ii};
        Ap_3 = AR_B*box.vertices(:,3) + Ao_B1{ii};
        Ap_4 = AR_B*box.vertices(:,4) + Ao_B1{ii};
        Ap_5 = AR_B*box.vertices(:,5) + Ao_B1{ii};
        Ap_6 = AR_B*box.vertices(:,6) + Ao_B1{ii};
        Ap_7 = AR_B*box.vertices(:,7) + Ao_B1{ii};
        Ap_8 = AR_B*box.vertices(:,8) + Ao_B1{ii};
        
        plot3([Ap_1(1) Ap_2(1)],[Ap_1(2) Ap_2(2)],[Ap_1(3) Ap_2(3)],'k');%
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
        
end