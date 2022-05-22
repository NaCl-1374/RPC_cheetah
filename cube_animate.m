%% ���ӻ�����
function cube_animate(X,time,foot_pos,swingstate,pre_contact,feetforce_used,fig,X_ref_mpc,path,mpcTra,Xtra,Polygon,View,dt,footSwingTrajectories,CLF)
%���ӻ����� �ο�ϵΪ����ϵ �Ӵ���Ϊ��ɫ���ڶ�Ϊ��ɫ
%����
%   X:��Ҫ���Ƶ�״̬����
%   time:��ʾ��ͼ�ϵ�ʱ��
%   foot_pos:�����λ��
%   swingstate:�ڶ�״̬
%   pre_contact:��ǰ�Ӵ�״̬
%   feetforce_used:������
%   fig:�������
%   X_ref_mpc:�ο��켣
%   path:����·��
%���
%   ��
figure(fig);
if CLF
clf; 
end
hold on; grid on; axis equal;
if ~isempty(View)
    view(View(1), View(2))
else
    view(-0, 0);%�ӽǿ��ƣ�ĿǰΪ��������˺󷽣��� 45 ���򿴻����� X(3)*180/pi
end
%���Ʋο��켣
if ~isempty(X_ref_mpc)
    plot_temp=[];
    for i=1:length(X_ref_mpc)/13
        plot_temp=[plot_temp,X_ref_mpc((i-1)*13+4:(i-1)*13+6)];
    end
    plot3(plot_temp(1,:),plot_temp(2,:),plot_temp(3,:),'*');%���Ʋο��켣,'linewidth',6
end
%%���Ƹ���·��
if ~isempty(path)
    plot3(path(:,1),path(:,2),ones(length(path(:,2)),1)*0.28);
end
if ~isempty(Polygon)
    plot3(Polygon(1),Polygon(2),0,'o');
end

%����mpcԤ��Ĺ켣
if ~isempty(mpcTra)
    plot_temp2=[];
    for i=1:length(mpcTra)/12
        plot_temp2=[plot_temp2,full(mpcTra((i-1)*12+4:(i-1)*12+6))];
    end
    plot3(plot_temp2(1,:),plot_temp2(2,:),plot_temp2(3,:),'o');
end

if ~isempty(Xtra)
    plot3(Xtra(4,:),Xtra(5,:),Xtra(6,:));
end

R = rotz(X(3))*roty(X(2))*rotx(X(1));%��ת����

plot3([X(4),X(4)+0.1*X(7)],[X(5),X(5)+0.1*X(8)],[X(6),X(6)+0.1*X(9)]);%��ǰ���ٶ�����
plot3([X(4),X(4)+0.1*X(10)],[X(5),X(5)+0.1*X(11)],[X(6),X(6)+0.1*X(12)]);%��ǰ�ٶ�����

plot_cube(R,0.38,0.22,0.10,[X(4),X(5),X(6)]');%���ƻ���

str=['pos: ',num2str(X(4:6)')];
str1=['rad: ',num2str(X(1:3)')];
if CLF
text(X(4)+0.25,X(5)+0.25,X(6)-0.25,str);%��עPOS/m
text(X(4)-0.25,X(5)-0.25,X(6)+0.25,str1);%��עrpy/rad
text(X(4),X(5)+0.5,X(6),num2str(time*dt(1)));%��עʱ��/s
end
%������ŵ�λ�á�r������������
for i=1:4
    x_indx=3*(i-1)+1;
    y_indx=3*(i-1)+2;
    z_indx=3*(i-1)+3;
    if ~isempty(footSwingTrajectories)
%         plot3(footSwingTrajectories(i).startPoint(1),footSwingTrajectories(i).startPoint(2),footSwingTrajectories(i).startPoint(3),'*','Color','r');
        plot3(footSwingTrajectories(i).endPoint(1),footSwingTrajectories(i).endPoint(2),footSwingTrajectories(i).endPoint(3),'o','Color','b');
        
    end
    if  swingstate(i)>0 && pre_contact(i)~=1
        
%         text(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),num2str(i),'Color','red');%������ע����� �ڶ��ú�ɫ
        %         plot3(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),'*');
        plot3([foot_pos(x_indx),foot_pos(x_indx)+0.01*feetforce_used(x_indx)],...
            [foot_pos(y_indx),foot_pos(y_indx)+0.01*feetforce_used(y_indx)],...
            [foot_pos(z_indx),foot_pos(z_indx)+0.01*feetforce_used(z_indx)]);
    else
        text(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),num2str(i),'Color','green');%������ע����� �Ӵ�����ɫ
        plot3([foot_pos(x_indx),foot_pos(x_indx)+0.01*feetforce_used(x_indx)],...
            [foot_pos(y_indx),foot_pos(y_indx)+0.01*feetforce_used(y_indx)],...
            [foot_pos(z_indx),foot_pos(z_indx)+0.01*feetforce_used(z_indx)]);
        plot3([X(4),foot_pos(x_indx)],[X(5),foot_pos(y_indx)],[X(6),foot_pos(z_indx)]);%ri
        
    end
end

axis([X(4)-1,X(4)+1,X(5)-1,X(5)+1,-0.1,1]);%����ϵ��Χ �Ե�ǰ����Ϊ����
xlabel('x');ylabel('y');zlabel('z');%�淶x,y,z������̶ȷ�Χ�����ڸ����������ϱ�ע��ĸx,y,z
drawnow;
end

function plot_cube(R, a, b, c,pos)
x = [-1,  1,  1, -1, -1,  1,  1, -1] * a/2;
y = [-1, -1,  1,  1, -1, -1,  1,  1] * b/2;
z = [-1, -1, -1, -1,  1,  1,  1,  1] * c/2;
P = R * [x; y; z]+pos;
order = [1, 2, 3, 4, 1, 5, 6, 7, 8, 5, 6, 2, 3, 7, 8, 4];
plot3(P(1, order), P(2, order), P(3, order), 'b');
end
%% ���ߺ���
function rotxm=rotx(theta)
s=sin(theta);
c=cos(theta);
% rotxm=[1,0,0;
%     0,c,s
%     0,-s c]';
rotxm=[1,0,0;
    0,c,-s
    0,s c];
end

function rotym=roty(theta)
s=sin(theta);
c=cos(theta);
% rotym =[c,0,-s;
%     0,1,0;
%     s,0,c]';
rotym =[c,0,s;
    0,1,0;
    -s,0,c];
end

function rotzm=rotz(theta)
s=sin(theta);
c=cos(theta);

% rotzm=[c,s,0;
%     -s,c,0;
%     0,0,1]';
rotzm=[c,-s,0;
    s,c,0;
    0,0,1];
end
%Rsb
function R=rotsb(theta)
% R=rotx(theta(1))*roty(theta(2))*rotz(theta(3));
R=rotz(theta(3))*roty(theta(2))*rotx(theta(1));

end

function s=Skew(in)
s=zeros(3,3);
s(1,2)=-in(3);
s(1,3)=in(2);
s(2,3)=-in(1);
s(2,1)=in(3);
s(3,1)=-in(2);
s(3,2)=in(1);
% s = [0 -in(3) in(2);
%     in(3) 0 -in(1);
%     -in(2) in(1) 0];
end

