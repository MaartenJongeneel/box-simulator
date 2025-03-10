function plotEnvironment(surface,dt,ii)

%Crossorter movement
% for jj = 1:11
%     surface{29+jj}.transform = surface{29+jj}.transform*expm(hat([0;2.5;0;0;0;0]*dt*ii));
% end

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