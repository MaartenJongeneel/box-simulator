clearvars; close all; clc;
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

c.dimd = 8; %discretization of the friction cone

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

%% %Check if body 1 is in contact with body 2
% Test case 1
AH_Bm(:,:,1) =   [1.0000         0         0         0
                   0   -1.0000   -0.0000         0
                   0    0.0000   -1.0000    1.0000
                   0         0         0    1.0000];
AH_Bm(:,:,2) = [    1.0000         0         0         0
                     0    1.0000         0         0
                     0         0    1.0000    1.0999999
                     0         0         0    1.0000];

AH_Bm(:,:,3) =   [1.0000         0         0         0
                   0   -1.0000   -0.0000         0
                   0    0.0000   -1.0000    4.0000
                   0         0         0    1.0000];
% 
% 
% %Test case 2
% AH_Bm(:,:,2) = AH_Bm(:,:,2)*[Rx(3)*Ry(3), [0;0;0.001]; zeros(1,3),1];

% %Test case 3
% AH_Bm(:,:,2) = AH_Bm(:,:,2)*[Rx(30)*Ry(3), [0;0;+0.055]; zeros(1,3),1];
% % 
%Test case 3
% AH_Bm(:,:,1) = AH_Bm(:,:,1)*[Rx(180), zeros(3,1); zeros(1,3),1];
% AH_Bm(:,:,2) = AH_Bm(:,:,2)*[Rx(5)*Ry(3), [0;0;0.012]; zeros(1,3),1];
% 
% Test case 4: Edge-edge
AH_Bm(:,:,2) = AH_Bm(:,:,2)*[Rx(5)*Ry(3), [0;0.15;0.015]; zeros(1,3),1];

%Test case 5: same orientation
% AH_Bm(:,:,2) = AH_Bm(:,:,2)*[Rx(0), [0;0;-0.001]; zeros(1,3),1];


% tic
% [contact,WNA,WTA,WNM,WTM,cp] = CollisionCheck(AH_Bm,AH_Bm,obj);
% toc

tic
[contact,WNA,WTA,WNM,WTM,cp] = CollisionCheckLCP(AH_Bm,AH_Bm,obj,c);
toc

figure;
% plotBox(AH_B{1}(:,:,ii),obj{1},[0 0 1]); hold on;
% plotBox(AH_B{2}(:,:,ii),obj{2},[0 0 1]); 
% plotBox(AH_B{3}(:,:,ii),obj{3},[0 0 1]);  
plotBox(AH_Bm(:,:,1),obj{1},[0 0 1]); hold on;
plotBox(AH_Bm(:,:,2),obj{2},[0 0 1]); 
plotBox(AH_Bm(:,:,3),obj{3},[0 0 1]); 
grid on; axis equal
axis off
axis([-0.4 0.5 -1 1 0 3.5]);    
hold off
view(-2,3);
drawnow
ax = gca;
ax.Clipping = "off";

% %Contact points in world frame used for plotting
ibref = 1;
Acp = cp;
hold on; plot3(Acp(1,:),Acp(2,:),Acp(3,:),'.','color',[1 0 0],'markersize',30)