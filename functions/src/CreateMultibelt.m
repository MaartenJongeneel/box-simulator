%% Functions
function mbelt = CreateMultibelt(nbelts,length,width,deg,vel,transform,id)

%MultiBelt Junction 60deg
% nbelts = 8; %number of multibelts
% length = 1.050; %longest belt length
% width = 0.760; %total multibelt width
% deg = 60;      %multibelt angle
% vel = 2.5;     %multibelt velocity in m/s

x = width/nbelts; %multibelts width
slength = length -(width-x)*tan(deg2rad(90-deg)); %Shortest multibelt length
y = linspace(length,slength,nbelts); %multibelt lengths


for ii = 1:nbelts
    mbelt{ii}.dim = [x y(ii)];
    mbelt{ii}.speed = [0; vel; 0];
    mbelt{ii}.transform = [eye(3) [(x/2)+(ii-1)*x; (ii-1)*(mean(diff(y))/2); 0]; zeros(1,3), 1];
    mbelt{ii}.id = "MB" + num2str(id) + "_"+num2str(ii);

    %transform origin of multibelt to middle in x, and zero in y
    mbelt{ii}.transform = mbelt{ii}.transform*[eye(3), [-0.5*width; 0.5*length; 0]; zeros(1,3),1];

    %transform multibelt to desired pose
    mbelt{ii}.transform = transform*mbelt{ii}.transform;
end
end