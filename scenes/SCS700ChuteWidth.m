%Now we have the base objects, we're going to adjust them to create the scene
% OBJ 1: sorter deck
obj{1}.dynamics = false; 
obj{1}.dim = [1.5 0.588 0.01];
obj{1}.initAH_B = eye(4);
obj{1}.initBV_AB = [0; 2.5; 0; 0; 0; 0];
% obj{1}.surface{1}.speed = [-1; 0; 0]; %Sorting speed

% OBJ 2: bridge plate
obj{2}.dynamics = false;
obj{2}.dim = [0.5; 6; 0.01];
obj{2}.initAH_B = [Ry(-8), [-1; 3; (sin(deg2rad(-8))*0.5) / 2]; zeros(1,3),1];
obj{2}.initBV_AB = zeros(6,1);

% OBJ 3: chute
obj{3}.dynamics = false;
obj{3}.dim = [5; 6; 0.01];
obj{3}.initAH_B = [Ry(-25), [-3.5; 3; obj{2}.initAH_B(3,4) + (sin(deg2rad(-25))*5) / 2]; zeros(1,3),1];
obj{3}.initBV_AB = zeros(6,1);

% OBJ 4: left guarding
obj{4}.dynamics = false;
obj{4}.dim = [5; 0.02; 0.3];
obj{4}.initAH_B = [Ry(-25), [-3.5; 6; obj{2}.initAH_B(3,4) + (sin(deg2rad(-25))*5) / 2]; zeros(1,3),1]  * [eye(3), [0;0; obj{4}.dim(3)/2]; zeros(1,3),1];
obj{4}.initBV_AB = zeros(6,1);

% OBJ 5: right guarding
obj{5}.dynamics = false;
obj{5}.dim = [5; 0.02; 0.3];
obj{5}.initAH_B = [Ry(-25), [-3.5; 0; obj{2}.initAH_B(3,4) + (sin(deg2rad(-25))*5) / 2]; zeros(1,3),1] * [eye(3), [0;0; obj{5}.dim(3)/2]; zeros(1,3),1];
obj{5}.initBV_AB = zeros(6,1);

% OBJ 6: Box
obj{6}.dynamics = true;
obj{6}.dim = [1; 0.4; 0.3];
obj{6}.initAH_B = [Rz(10), [0; 0; 0.25]; zeros(1,3), 1];
obj{6}.initBV_AB = zeros(6,1);

%% For testing purposes
for ii = 1:6 %3 objects
    obj{ii}.B_M_B = diag([1,1,1,0.0021,0.001,0.0027]);
    obj{ii}.mass = 1;
    % obj{ii}.initAH_B = [eye(3) [0; 0; 0;]; zeros(1,3),1];
    % obj{ii}.initBV_AB = zeros(6,1);
    % obj{ii}.dim = [1 1 1];

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
    % obj{ii}.dynamics = true;
end