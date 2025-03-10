%% Box release position
box.release.orientation = eye(3);  %Release orientation of the box            [deg]
box.release.position = [0 0 0.2];  %Release position of the box               [m]
box.release.linVel = [0 0 0];      %Release linear velocity (expressed in B)  [m/s]
box.release.angVel = [3 1 0];      %Release angular velocity (expressed in B) [rad/s]
box.parameters.eN = 0.4;           %Normal coefficient of restitution         [-]
box.parameters.eT = 0.0;           %Tangential coefficient of restitution     [-]s
box.parameters.mu = 0.5;           %Coefficient of friction                   [-]
box.mass = 1;                      %Box mass                                  [kg]
box.dimensions = [0.1 0.15 0.05];  %Box dimensions                            [m]
box.inertia_tensor = [eye(3), zeros(3,3); zeros(3,3), [0.0021, 0, 0; 0, 0.001, 0; 0, 0, 0.0027]];  %Box inertia tensor
box.discretization = 4;            %Box discretization of the contact points  [-]

%% Environment
surface{1}.dim = [1 2];            %Dimension of the surface                  [m]
surface{1}.speed = [0; 1; 0];      %Speed of the surface                      [m/s]
surface{1}.transform = [Rx(10)*Ry(5), [0; 0.5; 0]; zeros(1,3), 1]; %4x4 Transformation matrix 

