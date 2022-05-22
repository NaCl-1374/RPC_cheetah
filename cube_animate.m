%% 可视化函数
function cube_animate(X,time,foot_pos,swingstate,pre_contact,feetforce_used,fig,X_ref_mpc,path,mpcTra,Xtra,Polygon,View,dt,footSwingTrajectories,CLF)
%可视化函数 参考系为世界系 接触脚为绿色，摆动为红色
%输入
%   X:需要绘制的状态向量
%   time:显示在图上的时间
%   foot_pos:落足点位置
%   swingstate:摆动状态
%   pre_contact:提前接触状态
%   feetforce_used:控制力
%   fig:窗口序号
%   X_ref_mpc:参考轨迹
%   path:跟踪路径
%输出
%   无
figure(fig);
if CLF
clf; 
end
hold on; grid on; axis equal;
if ~isempty(View)
    view(View(1), View(2))
else
    view(-0, 0);%视角控制，目前为跟随机器人后方，从 45 方向看机器人 X(3)*180/pi
end
%绘制参考轨迹
if ~isempty(X_ref_mpc)
    plot_temp=[];
    for i=1:length(X_ref_mpc)/13
        plot_temp=[plot_temp,X_ref_mpc((i-1)*13+4:(i-1)*13+6)];
    end
    plot3(plot_temp(1,:),plot_temp(2,:),plot_temp(3,:),'*');%绘制参考轨迹,'linewidth',6
end
%%绘制跟踪路径
if ~isempty(path)
    plot3(path(:,1),path(:,2),ones(length(path(:,2)),1)*0.28);
end
if ~isempty(Polygon)
    plot3(Polygon(1),Polygon(2),0,'o');
end

%绘制mpc预测的轨迹
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

R = rotz(X(3))*roty(X(2))*rotx(X(1));%旋转矩阵

plot3([X(4),X(4)+0.1*X(7)],[X(5),X(5)+0.1*X(8)],[X(6),X(6)+0.1*X(9)]);%当前角速度向量
plot3([X(4),X(4)+0.1*X(10)],[X(5),X(5)+0.1*X(11)],[X(6),X(6)+0.1*X(12)]);%当前速度向量

plot_cube(R,0.38,0.22,0.10,[X(4),X(5),X(6)]');%绘制机身

str=['pos: ',num2str(X(4:6)')];
str1=['rad: ',num2str(X(1:3)')];
if CLF
text(X(4)+0.25,X(5)+0.25,X(6)-0.25,str);%标注POS/m
text(X(4)-0.25,X(5)-0.25,X(6)+0.25,str1);%标注rpy/rad
text(X(4),X(5)+0.5,X(6),num2str(time*dt(1)));%标注时间/s
end
%绘制落脚点位置、r向量、力向量
for i=1:4
    x_indx=3*(i-1)+1;
    y_indx=3*(i-1)+2;
    z_indx=3*(i-1)+3;
    if ~isempty(footSwingTrajectories)
%         plot3(footSwingTrajectories(i).startPoint(1),footSwingTrajectories(i).startPoint(2),footSwingTrajectories(i).startPoint(3),'*','Color','r');
        plot3(footSwingTrajectories(i).endPoint(1),footSwingTrajectories(i).endPoint(2),footSwingTrajectories(i).endPoint(3),'o','Color','b');
        
    end
    if  swingstate(i)>0 && pre_contact(i)~=1
        
%         text(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),num2str(i),'Color','red');%在落点标注脚序号 摆动用红色
        %         plot3(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),'*');
        plot3([foot_pos(x_indx),foot_pos(x_indx)+0.01*feetforce_used(x_indx)],...
            [foot_pos(y_indx),foot_pos(y_indx)+0.01*feetforce_used(y_indx)],...
            [foot_pos(z_indx),foot_pos(z_indx)+0.01*feetforce_used(z_indx)]);
    else
        text(foot_pos(x_indx),foot_pos(y_indx),foot_pos(z_indx),num2str(i),'Color','green');%在落点标注脚序号 接触用绿色
        plot3([foot_pos(x_indx),foot_pos(x_indx)+0.01*feetforce_used(x_indx)],...
            [foot_pos(y_indx),foot_pos(y_indx)+0.01*feetforce_used(y_indx)],...
            [foot_pos(z_indx),foot_pos(z_indx)+0.01*feetforce_used(z_indx)]);
        plot3([X(4),foot_pos(x_indx)],[X(5),foot_pos(y_indx)],[X(6),foot_pos(z_indx)]);%ri
        
    end
end

axis([X(4)-1,X(4)+1,X(5)-1,X(5)+1,-0.1,1]);%坐标系范围 以当前机身为中心
xlabel('x');ylabel('y');zlabel('z');%规范x,y,z坐标轴刻度范围，及在各自坐标轴上标注字母x,y,z
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
%% 工具函数
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

