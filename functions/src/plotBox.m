function Bplot = plotBox(AH_B,box,color) 
% Box-simulator-FixedPoint:
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
%
% Copyright (c) 2021, Maarten Jongeneel
% All rights reserved.
%% Plot the box
        AR_B = AH_B(1:3,1:3);
        %Output the position of the current time step for plotting purposes
        q(:,1)  = AH_B(1:3,4);
        R1(:,1) = AH_B(1:3,1);
        R2(:,1) = AH_B(1:3,2);
        R3(:,1) = AH_B(1:3,3);
        
        %Plot the origin of the box with its unit vectors
%         tip = [q(:,1)+ 0.3*R1(:,1) q(:,1)+ 0.3*R2(:,1) q(:,1)+ 0.3*R3(:,1)];
%         plot3([q(1,1) tip(1,1)],[q(2,1) tip(2,1)],[q(3,1) tip(3,1)],'r'); hold on
%         plot3([q(1,1) tip(1,2)],[q(2,1) tip(2,2)],[q(3,1) tip(3,2)],'g');
%         plot3([q(1,1) tip(1,3)],[q(2,1) tip(2,3)],[q(3,1) tip(3,3)],'b');
        
        %Create the box
        pbool = (abs(box.vertices(1,:))==max(abs(box.vertices(1,:))))&(abs(box.vertices(2,:))==max(abs(box.vertices(2,:))))&(abs(box.vertices(3,:))==max(abs(box.vertices(3,:))));
        Ap = AR_B*box.vertices(:,pbool)+AH_B(1:3,4);
%         Ap = AR_B*box.vertices+AH_B(1:3,4);
        Ap_1 = Ap(:,1);
        Ap_2 = Ap(:,2);
        Ap_3 = Ap(:,6);
        Ap_4 = Ap(:,5);
        Ap_5 = Ap(:,3);
        Ap_6 = Ap(:,4);
        Ap_7 = Ap(:,8);
        Ap_8 = Ap(:,7);
        
        plot3([Ap_1(1) Ap_2(1)],[Ap_1(2) Ap_2(2)],[Ap_1(3) Ap_2(3)],'k');hold on;
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

        %Plot all contact points
%         vertices = AR_B*box.vertices+AH_B(1:3,4);
%         plot3(vertices(1,:),vertices(2,:),vertices(3,:),'.',MarkerSize=1)
        
end