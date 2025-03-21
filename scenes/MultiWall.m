%% For testing purposes
for ii = 1:10 %3 objects
    obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
    obj{ii}.mass = 0.1;
    obj{ii}.initAH_B = [eye(3) [0; 0.2*ii; 0.05]; zeros(1,3),1]; %[eye(3) [0; 0.2*ii+0.01; 0.05]; zeros(1,3),1];
    obj{ii}.initBV_AB = zeros(6,1);
    obj{ii}.dim = [0.1;0.19;0.1];
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

for ii = 6:10
    obj{ii}.initAH_B = [eye(3) [0; 0.1+0.2*(ii-5); 0.155]; zeros(1,3),1];
end

% for ii = 11:15
%     obj{ii}.initAH_B = [eye(3) [0; 0.2*(ii-10); 0.26]; zeros(1,3),1];
% end
% 
% for ii = 31:40
%     obj{ii}.initAH_B = [eye(3) [0; 0.1+0.2*(ii-30); 0.35]; zeros(1,3),1];
% end


%define ground plane
ii = ii+1;

obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
obj{ii}.mass = 1;
obj{ii}.initAH_B = [eye(3), [0;0;-0.05]; zeros(1,3),1];
obj{ii}.initBV_AB = zeros(6,1);
obj{ii}.dim = [10;10;0.1];
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

% %% Define inflying object
% %define ground plane
% ii = ii+1;
% 
% obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
% obj{ii}.mass = 2;
% obj{ii}.initAH_B = [eye(3) [1.5; 0.65; 0.5]; zeros(1,3),1];
% obj{ii}.initBV_AB = [-6;0;0;2;0;0];
% obj{ii}.dim = [0.1;0.2;0.1];
% 
% 
% %Go over the surfaces
% obj{ii}.surface{1}.transform = [eye(3), [0; 0; obj{ii}.dim(3)/2]; zeros(1,3),1 ];
% obj{ii}.surface{2}.transform = [Ry(90), [obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
% obj{ii}.surface{3}.transform = [Ry(180), [0; 0; -obj{ii}.dim(3)/2]; zeros(1,3),1 ];
% obj{ii}.surface{4}.transform = [Ry(-90), [-obj{ii}.dim(1)/2; 0; 0]; zeros(1,3),1 ];
% obj{ii}.surface{5}.transform = [Rx(90), [0; -obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
% obj{ii}.surface{6}.transform = [-Rx(90), [0; obj{ii}.dim(2)/2; 0]; zeros(1,3),1 ];
% obj{ii}.surface{1}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
% obj{ii}.surface{2}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
% obj{ii}.surface{3}.dim=[obj{ii}.dim(1) obj{ii}.dim(2)];
% obj{ii}.surface{4}.dim=[obj{ii}.dim(3) obj{ii}.dim(2)];
% obj{ii}.surface{5}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];
% obj{ii}.surface{6}.dim=[obj{ii}.dim(1) obj{ii}.dim(3)];
% 
% %Contact points
% Ndisc=2;
% [X,Y,Z]=meshgrid(linspace(-obj{ii}.dim(1)/2,obj{ii}.dim(1)/2,Ndisc),linspace(-obj{ii}.dim(2)/2,obj{ii}.dim(2)/2,Ndisc),linspace(-obj{ii}.dim(3)/2,obj{ii}.dim(3)/2,Ndisc));
% pbool = (abs(X(:))==obj{ii}.dim(1)/2) | (abs(Y(:))==obj{ii}.dim(2)/2) | (abs(Z(:))==obj{ii}.dim(3)/2);
% obj{ii}.vertices= [X(pbool)';Y(pbool)';Z(pbool)'];
% 
% %Determine if the object has dynamics
% obj{ii}.dynamics = true;