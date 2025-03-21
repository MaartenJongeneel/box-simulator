function [contact,WNA,WTA,WNM,WTM,CP] = CollisionCheck(AH_B,AH_Bm,obj)
%Function to check if there is a collision between bodies. If there is, it will compute the matrices
%containing the force directions based on a set of computed contact points. 
%
% INPUTS:    AH_B                : 4x4xN double, transformation matrices of the N bodies at tA   [-]
%            AH_Bm               : 4x4xN double, transformation matrices of the N bodies at tM   [-]
%            obj                 : 1xN struct array, with fields of box properties as
%                                  obj{1}.B_M_B    : 6x6 double intertia tensor of the box       []  
%                                  obj{1}.mass     : 1x1 double mass of the box                  [kg]
%                                  obj{1}.initAH_B : 4x4 double, initial pose of the object      [-]
%                                  obj{1}.initBV_AB: 6x1 double, initial velocity of object      [-]
%                                  obj{1}.dim      : 3x1 double, dimensions of the object        [-]
%                                  obj{1}.surface  : 6x1 struct array with surface properties    [-]
%                                  obj{1}.vertices : 3xN double position of the contact          [m]
%                                                points w.r.t the body-fixed frame
%                                  obj{1}.dynamics : boolean, indicates if object has dynamics   [-]
%
% OUTPUTS:   contact             : vector of booleans indicating contact between 
%                                  the different bodies                                          [-] 
%            WNA                 : 6Nx1 double, matrix containing normal force directions at tA  [-]
%            WTA                 : 6Nx2 double, matrix containing tangent force directions at tA [-]
%            WNM                 : 6Nx1 double, matrix containing normal force directions at tM  [-]
%            WTM                 : 6Nx2 double, matrix containing tangent force directions at tM [-]
%            cp                  : 3xNcp double, matrix of contact points expressed in frame A   [m]
%
% Copyright (c) 2025, Maarten Jongeneel
% All rights reserved.
%% SAT test

Nobj = length(obj);     %Amount of objects in the scene

WNM = zeros(6*Nobj,1); %Preallocate: each contact point is one column with 6xNobj rows
WTM = zeros(6*Nobj,2); %Preallocate: each contact point is two columns with 6xNobj rows
WNA = zeros(6*Nobj,1); %Preallocate: each contact point is one column with 6xNobj rows
WTA = zeros(6*Nobj,2); %Preallocate: each contact point is two columns with 6xNobj rows
CP = []; %All contact points resulting from the total SAT test

%If we have only 1 body in the scene, nothing will happen. 
if Nobj <2    
    return;
end

%Create the vector of objects we need to check contact for
%each body needs to be checked with each other body only once
vecObj = 1:Nobj;
Nc     = sum(1:(Nobj-1)); %Number of contact checks we do
vecCheck = [];
for ii = 1:length(vecObj)-1
    vecCheck = [vecCheck [repmat(vecObj(ii),1,vecObj(end-ii)); vecObj(ii+1:end)]]; 
end

%Now we are going to build the matrices of the force direction
cntr = 1; %To count the contact points (and hence the size of the matrices containing the force directions)
contact = zeros(1,Nc); %Preallocate to speed up
for II = 1:length(vecCheck) %Loop through each body-to-body check
    contact(II) = false;

    %Get the body indexes we are going to check for
    ib1 = vecCheck(1,II); %Body 1 index
    ib2 = vecCheck(2,II); %Body 2 index

    %The 15 vectors we need to check for SAT test are based on principle axis of each object
    b1v1 = AH_Bm(1:3,1,ib1); b1v2 = AH_Bm(1:3,2,ib1); b1v3 = AH_Bm(1:3,3,ib1);
    b2v1 = AH_Bm(1:3,1,ib2); b2v2 = AH_Bm(1:3,2,ib2); b2v3 = AH_Bm(1:3,3,ib2);
    
    %The first six are the principle axis of each object, the remaining 9 are the cross products of these vectors
    SATvec = [b1v1/norm(b1v1) b1v2/norm(b1v2) b1v3/norm(b1v3) b2v1/norm(b2v1) b2v2/norm(b2v2) b2v3/norm(b2v3) cross(b1v1,b2v1)/norm(cross(b1v1,b2v1)) cross(b1v1,b2v2)/norm(cross(b1v1,b2v2)) cross(b1v1,b2v3)/norm(cross(b1v1,b2v3))...
        cross(b1v2,b2v1)/norm(cross(b1v2,b2v1)) cross(b1v2,b2v2)/norm(cross(b1v2,b2v2)) cross(b1v2,b2v3)/norm(cross(b1v2,b2v3)) cross(b1v3,b2v1)/norm(cross(b1v3,b2v1)) cross(b1v3,b2v2)/norm(cross(b1v3,b2v2)) cross(b1v3,b2v3)/norm(cross(b1v3,b2v3))];
    %Tracking which vector is used for each body (1=x, 2=y, 3=z vector)
    SATidx = [1 2 3 0 0 0 1 1 1 2 2 2 3 3 3; 0 0 0 1 2 3 1 2 3 1 2 3 1 2 3];     
    
    %Distance between the centers of the bodies in world frame
    dis = (AH_Bm(1:3,4,ib2)-AH_Bm(1:3,4,ib1));
    
    %Different from Erleben, we should instead use the full vertices positions (instead of half of the dimension) to compute seperation s:
    for ii = 1:length(SATvec)
        s(ii) = abs(dis'*SATvec(:,ii)) - ( abs(min((AH_Bm(1:3,1:3,ib1)*obj{ib1}.vertices)'*SATvec(:,ii))) + abs(min((AH_Bm(1:3,1:3,ib2)*obj{ib2}.vertices)'*SATvec(:,ii))) );
    end
    
    %Check if we have contact. In that case, they all must show overlap (there exist no axis that shows separation)
    contact(II) = sum(s(~isnan(s))<0) == sum(~isnan(s));  %all values that are not NaN should be smaller than zero
    cp = []; %Reset contact points within one check
    f = [];
    B1cp = [];
    
    if contact(II)   
        %Find the axis of seperation. We ignore the ones that are NaN (occurs in case of exact same orientation)
        isat = find(s(~isnan(s))==max(s(~isnan(s))),1);
        
        %If the collision along the first 3 axis, the first body has the reference face. If the collision is along the next 3 axis, the second body has the reference face. Otherwise, we have a edge-edge collision, and we need to retreive the axis
        if isat <7
            if isat<4 %Body 1 has the reference face                
                ibref = ib1;    ibinc = ib2;
            else      %Body 2 has the reference face
                ibref = ib2;    ibinc = ib1;
            end
            %The surface of the reference box that is in contact, is the one that closest to the indicent body. So let's see which face of ibref is closest to ibinc
            for ii = 1:length(obj{ibref}.surface)                
                f(ii) = (AH_Bm(1:3,4,ibref) + AH_Bm(1:3,1:3,ibref)*obj{ibref}.surface{ii}.transform(1:3,4) - AH_Bm(1:3,4,ibinc))' * SATvec(:,isat);
            end
            fidx = abs(f)==min(abs(f)); %Face index of reference body that is closest to incident body
            AH_C = AH_Bm(:,:,ibref)*obj{ibref}.surface{fidx}.transform; %The transform of that surface is our surface transform
    
            %Locations of the vertices of incident body expressed in reference body
            Cp   = AH_C(1:3,1:3)'*(AH_Bm(1:3,4,ibinc) + AH_Bm(1:3,1:3,ibinc)*obj{ibinc}.vertices-AH_C(1:3,4)); %Gap function 
            B1p2 = obj{ibref}.surface{fidx}.transform(1:3,1:3) * Cp + obj{ibref}.surface{fidx}.transform(1:3,4);
            B1p1 = obj{ibref}.vertices;

            %Find the vertices that lie below reference plane
            p_i = find(Cp(3,:)<0); %Vertices that lie below reference plane

            %Here we do some clipping of the contact points to the surfaces of interest
            for ii = 1:length(p_i) %loop over the contact points
                for id = 1:2 %loop over the two dimensions of the surface
                    if abs(B1p2(id,p_i(ii)))<(obj{ibref}.dim(id)/2) %If the point is inside the dimension of the object
                        B1cp(id,ii) = B1p2(id,p_i(ii)); %Then we use its coordinates
                    else
                        index = find(vecnorm(B1p2(:,p_i(ii))-B1p1)==min(vecnorm(B1p2(:,p_i(ii))-B1p1))); %Find the contact point on reference body 1 closest to body 2
                        B1cp(id,ii) = B1p1(id,index); %If not, we use the reference body vertex
                    end
                end
                B1cp(3,ii) = B1p2(3,p_i(ii)); %Their z-value is projected to the reference plane
            end

            %Express the contact points in world frame and incident body
            cp = AH_Bm(1:3,1:3,ibref)*B1cp + AH_Bm(1:3,4,ibref);
            AR_C = AH_C(1:3,1:3);
            B2cp = AH_Bm(1:3,1:3,ibinc)\(cp - AH_Bm(1:3,4,ibinc)); %B2H_B1(1:3,1:3)*B1cp + B2H_B1(1:3,4); %its the incident body's contact points that we consider
            bi1 = ibref; %Put the body index
            bi2 = ibinc; %Put the body index          
        else
            %Here we have an edge-edge collision and only a single contact point
            %Find the edges that are colliding
            v1 = AH_Bm(1:3,SATidx(1,isat),ib1); %Using the index vector to get back which vector was used to create the cross product that is now giving the smallest overlap
            v2 = AH_Bm(1:3,SATidx(2,isat),ib2);
            for jj = 1:2 %loop over the two bodies
                pbool = (abs(obj{vecCheck(jj,II)}.vertices(1,:))==max(abs(obj{vecCheck(jj,II)}.vertices(1,:))))&(abs(obj{vecCheck(jj,II)}.vertices(2,:))==max(abs(obj{vecCheck(jj,II)}.vertices(2,:))))&(abs(obj{vecCheck(jj,II)}.vertices(3,:))==max(abs(obj{vecCheck(jj,II)}.vertices(3,:))));
                % Ap = AH_Bm{vecCheck(jj,II)}(1:3,1:3)*obj{vecCheck(jj,II)}.vertices(:,pbool); %+AH_Bm{vecCheck(jj,II)}(1:3,4);
                Ap = obj{vecCheck(jj,II)}.vertices(:,pbool); %So now we are just considering edges in body frame

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
                    Aedges = AH_Bm(1:3,1:3,ib1)*(edges./vecnorm(edges)); %Edges expressed in world frame and normalized
                    for ii = 1:length(Aedges)
                        ed{jj}(ii) = norm(cross(v1,Aedges(:,ii)));
                    end
                    % idx = find(~round(ed{jj})); %Indices of edges that are parallel to the axis used for computing axis of separation
                    idx = find(ed{jj}==min(ed{jj})); %Safer way of computing it
                    B2H_B1 = AH_Bm(:,:,ib2)\AH_Bm(:,:,ib1); %Transformation from Body 1 to Body 2
                    B2p1 = B2H_B1(1:3,1:3)*reshape(E(:,:,idx),3,8) +B2H_B1(1:3,4); %Expressing the edges of body 1 in body 2
                    dis2 = [mean(B2p1(:,1:2),2), mean(B2p1(:,3:4),2), mean(B2p1(:,5:6),2), mean(B2p1(:,7:8),2)]; %Vectors from body 2 to the edges of body 1
                    dis3 = AH_Bm(1:3,1:3,ib2)*dis2; %Expressing these vectors in frame A
                    dis4 = abs(dis3'*SATvec(:,isat)); %These vectors projected on the separating axis
                    idx2 = find(dis4==min(dis4)); %Index of parallel edge closest to other body
                    % eA = AH_Bm(1:3,1:3,ib1)*edges(:,idx(idx2)); %Edge eA from body A closest to body B
                    P1 = AH_Bm(1:3,1:3,ib1)*E(:,1,idx(idx2))+AH_Bm(1:3,4,ib1);
                    V1 = AH_Bm(1:3,1:3,ib1)*edges(:,idx(idx2));
                elseif jj==2
                    Aedges = AH_Bm(1:3,1:3,ib2)*(edges./vecnorm(edges)); %Edges expressed in world frame and normalized
                    for ii = 1:length(Aedges)
                        ed{jj}(ii) = norm(cross(v2,Aedges(:,ii)));
                    end
                    % idx = find(~round(ed{jj})); %Indices of edges that are parallel to the axis used for computing axis of separation
                    idx = find(ed{jj}==min(ed{jj})); %Safer way of computing it
                    B1H_B2 = AH_Bm(:,:,ib1)\AH_Bm(:,:,ib2); %Transformation from Body 2 to Body 1
                    B1p2 = B1H_B2(1:3,1:3)*reshape(E(:,:,idx),3,8) +B1H_B2(1:3,4); %Expressing the edges of body 2 in body 1
                    dis2 = [mean(B1p2(:,1:2),2), mean(B1p2(:,3:4),2), mean(B1p2(:,5:6),2), mean(B1p2(:,7:8),2)]; %Distance from body 1 to the edges of body 2
                    dis3 = AH_Bm(1:3,1:3,ib1)*dis2; %Expressing these vectors in frame A
                    dis4 = abs(dis3'*SATvec(:,isat)); %These vectors projected on the separating axis
                    idx2 = find(dis4==min(dis4)); %Index of parallel edge closest to other body
                    % eB = AH_Bm(1:3,1:3,ib2)*edges(:,idx(idx2)); %Edge eA from body A closest to body B
                    P2 = AH_Bm(1:3,1:3,ib2)*E(:,1,idx(idx2))+AH_Bm(1:3,4,ib2);
                    V2 = AH_Bm(1:3,1:3,ib2)*edges(:,idx(idx2));
                end
            end
            % normal = cross(eA,eB)/norm(cross(eA,eB)); %Erleben eq. 87
            %Alternative, we know that the normal is the same as the separating axis SATvec(:,isat) ... 
            normal = SATvec(:,isat);
            b = s(isat)*SATvec(:,isat)-P2+P1;
            A = [-V1(:,1) V2(:,1)]; %Taking first column. If there are more columns, their values must be the same, as it means the distance to those edges are the same.
            t = linsolve(A,b);
            L1 = P1+t(1)*V1; %Contact point on body 1 expressed in world frame
            L2 = P2+t(2)*V2; %Contact point on body 1 expressed in world frame
            cp1 = 0.5*(L1+L2); %Single contact point in midpoint of overlapping region expressed in world frame.
            B1cp = AH_Bm(1:3,1:3,ib1)'*(cp1 - AH_Bm(1:3,4,ib1)); %Single contact point exprssed in body 1
            cp = AH_Bm(1:3,1:3,ib1)*B1cp + AH_Bm(1:3,4,ib1);
            B2cp = AH_Bm(1:3,1:3,ib2)\(cp - AH_Bm(1:3,4,ib2)); %its the incident body's contact points that we consider

            %If the projected contact point distance in body 1 is negative on the axis of separation, we need to flip the axis
            if (AH_Bm(1:3,1:3,ib1)*B1cp)'*normal <0
                normal = -normal;
            end
            AR_C = [null(normal') normal]; % We assume the convention of having the normal of the reference face as positive
            bi1 = ib1; %Put the body index
            bi2 = ib2; %Put the body index
        end

        %Now we need to create the matrices containing the force directions
        for icp = 1:length(B1cp(1,:)) %Count over the contact points between these two bodies

            wB1 = (AR_C'*[-AH_Bm(1:3,1:3,bi1) AH_Bm(1:3,1:3,bi1)*hat(B1cp(:,icp))])';     %This is the reaction force, in oposite direction
            wB1a= (AR_C'*[-AH_B(1:3,1:3,bi1) AH_B(1:3,1:3,bi1)*hat(B1cp(:,icp))])'; %This is the reaction force, in oposite direction
            wB2 = (AR_C'*[AH_Bm(1:3,1:3,bi2) -AH_Bm(1:3,1:3,bi2)*hat(B2cp(:,icp))])';     %Starndard convention for positive normal
            wB2a= (AR_C'*[AH_B(1:3,1:3,bi2) -AH_B(1:3,1:3,bi2)*hat(B2cp(:,icp))])'; %Starndard convention for positive normal

            if obj{bi1}.dynamics
                WNM((6*(bi1-1)+1:6*(bi1-1)+6),cntr) = wB1(:,3); %For the ibref body we compute the force directions
                WTM((6*(bi1-1)+1:6*(bi1-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1(:,1:2); %Assuming this works the same way in the tangential direction

                WNA((6*(bi1-1)+1:6*(bi1-1)+6),cntr) = wB1a(:,3); %For the ibref body we compute the force directions
                WTA((6*(bi1-1)+1:6*(bi1-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB1a(:,1:2); %Assuming this works the same way in the tangential direction
            end

            if obj{bi2}.dynamics
                WNM((6*(bi2-1)+1:6*(bi2-1)+6),cntr) = wB2(:,3); %For the ibinc body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                WTM((6*(bi2-1)+1:6*(bi2-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2(:,1:2);

                WNA((6*(bi2-1)+1:6*(bi2-1)+6),cntr) = wB2a(:,3); %For the ibinc body it is in contact with, the force direction is oposite (see also eq 17 in K.Erleben, Contact and Friction Simulation for Computer Graphics.)
                WTA((6*(bi2-1)+1:6*(bi2-1)+6),(2*(cntr-1)+1:2*(cntr-1)+2)) = wB2a(:,1:2);
            end

            cntr = cntr+1;
        end
    end
    CP = [CP cp];
end
end