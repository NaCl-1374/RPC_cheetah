%% 用来测试非线性模型预测控制
clc;
clear;
warning off
addpath('E:\Matlab\casadi-matlab')
import casadi.*
%%
mpcHz=60;%mpc频率
delta_t = 0.002;%仿真时间间隔
IterMPCRate=fix(1/mpcHz/delta_t);%mpc频率控制
N=8;
runtime=14;
%% %物理参数

Body.m = 9;%机器人质量
%机身惯量
Body.Ib = diag([0.07,0.26,0.242]);%转动惯量矩阵
Body.length_body=0.38;
Body.width_body=0.22;
Body.hipPos=diag([1 1 1, 1 -1 1, -1 1 1, -1 -1 1])*repmat([0.2 0.1 0],1,4)';
world.g = 9.8;%重力加速度
world.mu=0.4;%摩擦系数

%%	步态参数
%一个步态片段时长
Gait.dt=0.05;
%一个步态片段在主循环中的iter数
Gait.iterationsBetween=fix(Gait.dt/delta_t);
dt_val = repmat(Gait.dt,1,N)';
%步态相位
Gait.phsae=0;
%步态片段数
Gait.nIterations=8;
Gait.Tp=Gait.nIterations*Gait.dt;
Gait.iterLength=Gait.iterationsBetween*Gait.nIterations;
%步态相位差片段计数形式
% offsets=[3,0,8,5]*5;
Gait.offsets=[0,4,4,0];
%步态相位差小数形式
Gait.offsetFloat=Gait.offsets/Gait.nIterations;
%每个的支撑持续时间时间片段计数
Gait.durations=ones(1,4)*Gait.nIterations*0.5;
%每个的支撑持续时间主循环时间片段计数
Gait.durationsFloat=Gait.durations/Gait.nIterations;
%步态片段计数时间片段计数
Gait.gait_iteration=0;
%支撑时间片段计数
Gait.standtime=Gait.durations(1);
%摆动时间片段计数
Gait.swingtime=Gait.nIterations-Gait.standtime;
%支撑时间 s
Gait.standTime=Gait.standtime*Gait.dt;%步态分段时间
%摆动时间 s
Gait.swingTime=Gait.swingtime*Gait.dt;
%%	落足点计算需要

firstSwing=ones(1,4);%开始摆动标志
swingTimeRemaining=ones(4,1)*Gait.swingTime;%剩余摆动时间
side_signy=[1,-1,1,-1];%落脚点偏执符号
side_signx=[1,1,-1,-1];%落脚点偏执符号
pre_contact=[0,0,0,0];%提前接触标志
footRef=Body.hipPos;
% 摆动轨迹起始点
footSwingTrajectories(1:4)=struct;
for i=1:4
    footSwingTrajectories(i).startPoint=zeros(3,1);
    footSwingTrajectories(i).endPoint=zeros(3,1);
end

%% %状态的权重

weight.QX = [5000 1000 2500, 0 0 5000, 2.5 2.5 300, 2500 2500 2500 ]';
weight.QN = [5000 1000 2500, 0 0 5000, 2.5 2.5 300, 2500 2500 2500 ]';
weight.Qc = [350 350 350]';
weight.Qf = [0.25 0.25 0.0001]';
%% 构造微分方程
X_=SX.sym('X_', 12, 1);
n_state=size(X_,1);
F_=SX.sym('F_', 12, 1);
n_F=size(F_,1);
R_=SX.sym('R_', 12, 1);
n_r=size(R_,1);
%% 计算微分方程
I3=eye(3);
Rbody=rotsb(X_(1:3));

cy = cos(X_(3));
sy = sin(X_(3));
cp = cos(X_(2));
sp = sin(X_(2));

R_yaw =[cy sy 0;
    -sy cy 0;
    0 0 1];%世界到机身
R_w=[cy/cp,sy/cp,0;
    -sy,cy,0;
    cy*sp/cp,sy*sp/cp,1];
Ig = Rbody*Body.Ib*Rbody';
Ig_inv=Ig\I3;
%hip在世界系

A = [zeros(3) zeros(3) R_yaw zeros(3)  ;
    zeros(3) zeros(3) zeros(3) I3 ;
    zeros(3) zeros(3) zeros(3) zeros(3);
    zeros(3) zeros(3) zeros(3) zeros(3) ;
    ];%状态矩阵
AA=A;
AA(1:3,7:9)=R_w;
B=[zeros(3)           zeros(3)           zeros(3)            zeros(3);
    zeros(3)           zeros(3)           zeros(3)            zeros(3);
    Ig_inv*Skew(R_(1:3)) Ig_inv*Skew(R_(4:6)) Ig_inv*Skew(R_(7:9))  Ig_inv*Skew(R_(10:12));
    I3/Body.m   I3/Body.m   I3/Body.m    I3/Body.m;];%控制矩阵
grav=zeros(12,1);
grav(12)=-world.g;

dotX=A*X_+B*F_+grav;
Dotx=AA*X_+B*F_+grav;%
%% 定义微分方程函数
% 简化
f=Function('f',{X_;F_;R_},{dotX},{'input_states','control_inputs','foot_input'},{'dotX'});
% 精确
fa=Function('fa',{X_;F_;R_},{Dotx},{'input_states','control_inputs','foot_input'},{'dotX'});

%%  构造代价和约束 变量定义
X = SX.sym('X', n_state, N+1); % N+1步状态
F = SX.sym('F', n_F, N); % N步内的控制
r = SX.sym('r', n_r, N); % N步内的控制

RefX = SX.sym('RefX', n_state, N+1); % N步内的控制输出
RefF = SX.sym('RefF', n_F, N); % N步内的控制输出
Refr = SX.sym('Refr', n_r, N); % N步内的控制输出

ContactState=SX.sym('ConState', 4, N);
obj=0;
% g=[];
%%  构造代价和约束 摩擦和运动学限制
mu_inv = 1.0/world.mu;
%摩擦约束
f_block =[ mu_inv, 0,  -1.0;
    -mu_inv, 0,  -1.0;
    0,  mu_inv, -1.0;
    0, -mu_inv, -1.0;];
%运动学约束
kin_box_x = 0.05;
kin_box_y = 0.05;
kin_box_z = 0.3;

Kin_block =[ 1, 0,  0,-kin_box_x;
    -1, 0,  0,-kin_box_x;
    0, 1, 0,-kin_box_y;
    0, -1, 0,-kin_box_y;
    0, 0, 1,0.06;
    0, 0, -1,-kin_box_z];
%足端换算用
Phip_swing=Body.hipPos;
Phip_swing([3,6,9,12])=[-0.2,-0.2,-0.2,-0.2];
%%  约束暂存变量定义 %初态约束

defect_init=X(:,1)-RefX(:,1);%12*1 eq

defect_state=SX.zeros(12*(N+1),1);%12(N+1)*1 eq
defect_FootOnGround=SX.zeros(4*(N),1);%4(N)*1 eq
defect_footStance=SX.zeros(12*(N),1);%(3*4)(N)*1 eq
n_equa_c=12+12*(N+1)+4*(N)+12*(N);%12+
%共
defect_legLimits=SX.zeros(24*(N),1);%(4*6)(N)*1
defect_footforce=SX.zeros(16*(N),1);%(4*4)(N)*1 xy摩擦约束4个
defect_ForceNormal=SX.zeros(N,1);% (N)*1
defect_footswing=SX.zeros(4*(N),1);%4(N)*1
n_inequa_c=24*(N)+16*(N)+N+4*(N);
%%	约束和代价计算
for k = 1:N
    %%	代价计算
    Xk=X(:,k);
    Fk=F(:,k);
    rk=r(:,k);
    Pk=repmat(Xk(4:6),4,1)+rk;
    ContactStatek=ContactState(:,k);
    dtk=dt_val(k);
    X_err = Xk - RefX(:,k);                                         % 基座状态误差
    pf_err = repmat(Xk(4:6),4,1) + Phip_swing - Pk;                      %  悬空时约束foot位置
    r_err=rk-Refr(:,k);
    U_err = Fk - RefF(:,k);                                         % GRF 误差
    obj = obj + (X_err'*diag(weight.QX)*X_err+...                     % 误差求和
        r_err'*diag(repmat(weight.Qc,4,1))*r_err+...
        U_err'*diag(repmat(weight.Qf,4,1))*U_err)*dtk;
    %%	约束计算
    %状态约束
    %% runge kutta method
    %     k1 = f(Xk,Fk,Pk);   % new
    %     k2 = f(Xk + dtk/2*k1,Fk,rk); % new
    %     k3 = f(Xk + dtk/2*k2,Fk,rk); % new
    %     k4 = f(Xk + dtk*k3,Fk,rk); % new
    %     st_next_RK4=Xk +dtk/6*(k1+2*k2+2*k3+k4); % new
    %     defect_state((k-1)*12+1:(k-1)*12+12)=X(:,k+1)-(st_next_RK4);
    %%
    defect_state((k-1)*12+1:(k-1)*12+12)=X(:,k+1)-(Xk+f(Xk,Fk,rk)*dtk);
    %法向力大于0 不等式
    defect_ForceNormal(k)=-dot(Fk,repmat([0;0;1],4,1));
    %结合法向力大于0，摩擦约束来约束摆动中力为0 和最大力 不等式
    defect_footswing((k-1)*4+1:(k-1)*4+4)=Fk([3,6,9,12])-ContactStatek.*repmat(1000,4,1);
    for leg=1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        %脚在地上约束 0是此时地面高度等式
        defect_FootOnGround((k-1)*4+leg)=ContactStatek(leg)*Pk(3*(leg-1)+3);
        %限制腿长 限制范围不等式
        Rbody=rotsb(Xk(1:3));
        Phip=Rbody*Body.hipPos(xyz_idx)+Xk(4:6);
        p_rel = (Pk(xyz_idx) - Phip);%hip->足端
        defect_legLimits((k-1)*24+(leg-1)*6+1:(k-1)*24+(leg-1)*6+6)= Kin_block*[p_rel;1];
        %接触中脚不滑动
        if (k < N)
            Pk1=repmat(X(4:6,k+1),4,1)+r(:,k+1);
            defect_footStance((k-1)*12+(leg-1)*3+1:(k-1)*12+(leg-1)*3+3)=ContactStatek(leg)*(Pk1(xyz_idx)-Pk(xyz_idx));%ContactState(leg,k+1)*
        end
        %摩擦约束 不等式
        defect_footforce((k-1)*16+(leg-1)*4+1:(k-1)*16+(leg-1)*4+4)=f_block*Fk(xyz_idx);
    end
end
%%	约束生成 defect_init;
g=[defect_init;defect_state;defect_FootOnGround;defect_footStance;...
    defect_legLimits;defect_footforce;defect_ForceNormal;defect_footswing];
display_str=['等式约束数量',num2str(n_equa_c),'   不等式约束数量',num2str(n_inequa_c)];
disp(display_str);
%%	终端 cost
X_err = X(:,end)-RefX(:,end);
obj = obj + X_err'*diag(weight.QN)*X_err;
%%	构造问题和问题变量
OPT_variables = [reshape(X,n_state*(N+1),1);reshape(F,n_F*N,1);reshape(r,n_r*N,1)];
OPT_Param = [reshape(RefX,n_state*(N+1),1);reshape(RefF,n_F*N,1);reshape(Refr,n_r*N,1);reshape(ContactState,4*N,1)];
nlp_prob =struct('f', obj, 'x',OPT_variables,'p',OPT_Param, 'g',g);
%%  优化设置
opts_setting=struct;
opts_setting.ipopt.max_iter=100;
opts_setting.ipopt.tol=1e-4;
opts_setting.ipopt.print_level=0;
opts_setting.ipopt.acceptable_tol=1e-4;
opts_setting.ipopt.constr_viol_tol=1e-3;
opts_setting.ipopt.acceptable_iter= 5;
opts_setting.ipopt.acceptable_obj_change_tol=1e-6;
opts_setting.ipopt.nlp_scaling_method='gradient-based';
opts_setting.ipopt.nlp_scaling_max_gradient=50; % (100), % 50 works well
opts_setting.ipopt.bound_relax_factor= 1e-6; % (1e-8), % 1e-6 works well
opts_setting.ipopt.fixed_variable_treatment='relax_bounds'; % {'make_parameter','make_constraint','relax_bounds'}; % relax bounds works well
opts_setting.ipopt.bound_frac=5e-1; % (1e-2), 5e-1 works well
opts_setting.ipopt.bound_push=5e-1; % (1e-2), 5e-1 works well
opts_setting.ipopt.mu_strategy='adaptive'; % {'monotone','adaptive'}; % adaptive works very well
opts_setting.ipopt.mu_oracle='probing'; % {'quality-function','probing','loqo'}; % probing works very well
opts_setting.ipopt.fixed_mu_oracle='probing'; % {'average_compl','quality-function','probing','loqo'}; % probing decent
opts_setting.ipopt.adaptive_mu_globalization='obj-constr-filter'; % {'obj-constr-filter','kkt-error','never-monotone-mode'};
opts_setting.ipopt.mu_init=1e-1; % [1e-1 1e-2 1]
opts_setting.ipopt.alpha_for_y='bound-mult'; % {'primal','bound-mult','min','max','full','min-dual-infeas','safer-min-dual-infeas','primal-and-full'}; % primal or bound-mult seems best
opts_setting.ipopt.alpha_for_y_tol=1e1; % (1e1)
opts_setting.ipopt.recalc_y='no'; % {'no','yes'};
opts_setting.ipopt.max_soc=4; % (4)
opts_setting.ipopt.accept_every_trial_step='no'; % {'no','yes'}
opts_setting.ipopt.linear_solver='mumps'; % {'ma27','mumps','ma57','ma77','ma86'} % ma57 seems to work well
opts_setting.ipopt.linear_system_scaling='slack-based'; {'mc19','none','slack-based'}; % Slack-based
opts_setting.ipopt.linear_scaling_on_demand='yes'; % {'yes','no'};
opts_setting.ipopt.max_refinement_steps=10; % (10)
opts_setting.ipopt.min_refinement_steps=1; % (1)
opts_setting.ipopt.warm_start_init_point= 'no'; % (no)

%% 构造求解器
solver = nlpsol('solver', 'ipopt', nlp_prob,opts_setting);
%%	约束上下界面 args
args.lbg(1:n_equa_c) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_equa_c) = 0;  % 1e-20   % Equality constraints

args.lbg(n_equa_c+1 : n_equa_c+ n_inequa_c) = -inf; % inequality constraints
args.ubg(n_equa_c+1 : n_equa_c+ n_inequa_c) = 0; % inequality constraints

%%  绘图用矩阵
contact_list=[];%绘制接触相位图
swing_list=[];%绘制摆动相位图

%%  输入速度
des_vel=[0.0,0,0]';%机身下速度
v_des_world=[0.6,0,1]';%世界下速度 最后一个必须为1
yaw_turn_rate=0.0;%yaw速度期望
yaw_des=0;%yaw期望
onmiMode=true;%是否无头模式，1 使用世界速度控制机器人，0 使用本体速度控制

%%  仿真变量
X_init = [0;0.0;0; 0.0;0.0;0.28 ;0;0;0; 0;0;0];%初始状态变量
X_out(:,1)=X_init;%仿真输出状态变量 相当于观测器 这个最后记录了整个过程所有的状态
P_foot=reshape(Body.hipPos,[],1);
pFoot=reshape(Body.hipPos,3,4);
pFoot_next=reshape(Body.hipPos,3,4);
swingTimeCount=0;
first_run=1;
pic_num = 1;%保存gif用
time=['NLP','_',datestr(datetime('now'),'yyyy-mm-dd-HH-MM'),'_trot.gif'];

%%  参考轨迹生成
X_ref=zeros(12*(N+1),1);
U_ref=zeros(12*(N),1);
r_ref=zeros(12*(N),1);


%%  主循环
for iter = 0:runtime/delta_t+1
    
%       if iter <500
%         v_des_world=[0;0;0];
%         yaw_turn_rate=0;
%       else
%         v_des_world=[0.6;0;0];
%         yaw_turn_rate=0.0;
%       end
    
    %% 速度转换
    Rbody=rotsb(X_out(1:3,iter+1));
    if onmiMode
        des_vel=Rbody'*v_des_world;
    else
        v_des_world=Rbody*des_vel;
    end
    
    %%  首次设置起始点
    if first_run
        for i=1:4
            footSwingTrajectories(i).startPoint=pFoot(:,i);
            footSwingTrajectories(i).endPoint=pFoot(:,i);
        end
%         first_run=0;
    else
        %         %%  仿真足端位置
        %         for i=1:4
        %             pFoot(i)=pFoot(i);
        %         end
    end
    %% 步态计算
    [Gait.gait_iteration,Gait.phsae]=setGaitIter(Gait,iter+1);%计算步态相位
    contactstate=getConectState(Gait);%计算接触相位
    swingstate=getSwingState(Gait);%计算摆动相位
    %contact_list=[contact_list;contactstate];%保存接触相位
    swing_list=[swing_list;swingstate];%保存摆动相位
    %% 剩余摆动时间计算
    for leg=1:4
        xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
        if firstSwing(leg)
            %此时在接触
            swingTimeRemaining(leg)=Gait.swingTime;
        else
            %此时在摆动
            swingTimeRemaining(leg)=swingTimeRemaining(leg)-delta_t;
        end
    end
    %%  控制计算 abs(Gait.phsae-0)<0.0126 ||abs(Gait.phsae-0.5)<0.0126
    if rem(iter,IterMPCRate)==0
        %%	落足点计算 参考mit代码
        for leg=1:4
            if swingstate(leg)>0 && pre_contact(leg)~=1%只计算摆动腿
                xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
                offest_hip=[0.0*side_signx(leg);0.0*side_signy(leg);0];%落足点偏置 可以实现腿内收外扩
                pRobotFrame=Body.hipPos(xyz_idx)-offest_hip;
                pYawCorrected=rotz(-yaw_turn_rate*Gait.standTime/2)*pRobotFrame;%负值为了促进旋转
                Pf=X_out(4:6,1)+Rbody*(pYawCorrected+des_vel*swingTimeRemaining(leg));
                
                pfx_rel = X_out(10,1) * (0.5 + 0) * Gait.standTime +0.03*(X_out(10,1)-v_des_world(1)) +(0.5* X_out(6,1)/9.81) * (X_out(11,1)*yaw_turn_rate);
                pfy_rel = X_out(11,1) * 0.5 * Gait.standTime +0.03*(X_out(11,1)-v_des_world(2)) +(0.5*X_out(6,1)/9.81) * (-X_out(10,1)*yaw_turn_rate);
                pfx_rel = min(max(pfx_rel, -0.3), 0.3);
                pfy_rel = min(max(pfy_rel, -0.3), 0.3);
                
                Pf(1)=Pf(1)+pfx_rel;
                Pf(2)=Pf(2)+pfy_rel;
                Pf(3)=0;
                %更新摆动终点
                %                 footSwingTrajectories(leg).endPoint=Pf;
                footRef(xyz_idx)=Pf;
            end
        end
        %% 步态时间表生成
        table=getMPCTable(Gait,N,zeros(4,1));
        %% 简单计算参考轨迹
        yaw_des=yaw_des+yaw_turn_rate*Gait.dt;
        X_ref=X_out(:,iter+1);
        X_to_sum = [ 0;0;yaw_turn_rate*Gait.dt;v_des_world(1)*Gait.dt;v_des_world(2)*Gait.dt;0.0; 0.0;0;0; 0;0;0];%参考值累计用
        %%如果跟踪一个给定轨迹 在这里设置 X_out(4,iter) X_out(5,iter)
        X_init_state=[0;0;yaw_des;X_out(4,iter+1);X_out(5,iter+1);0.28;0;0;yaw_turn_rate;v_des_world(1);v_des_world(2);0];%参考值初值
        X_ref(13:24)=X_init_state;
        %速度累计计算参考轨迹
        for i =3:N+1
            X_idx=(i-1)*12+1:(i-1)*12+12;
            last_X_idx=(i-2)*12+1:(i-2)*12+12;
            X_ref(X_idx)=X_ref(last_X_idx)+X_to_sum;
        end
        if first_run
            U_ref(1:12)=zeros(12,1);
            r_ref(1:12)=zeros(12,1);
        first_run=0;
        else
            U_ref(1:12)=F_use;
            r_ref(1:12)=rin;
        end
        for i =2:N
            X_idx=(i-1)*12+1:(i-1)*12+12;
            Csk=table(i,:);
            Xk=X_ref(X_idx);
            r_k=zeros(12,1); 
            U_k=zeros(12,1);
            for leg=1:4
                xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
                r_k(xyz_idx)=Csk(leg)*(footRef(xyz_idx)-Xk(4:6));
                U_k(xyz_idx)=((Body.m*Gait.Tp)/(Gait.standTime*sum(Csk)))*...
                    [cos(Xk(3))*Xk(9),0,0;
                    0,sin(Xk(3))*Xk(9),0;
                    0,0,world.g]*abs(v_des_world);
            end
            U_ref(X_idx)=U_k;
            r_ref(X_idx)=r_k;
        end
        %%  求解 准备参考 初值
        args.p=[X_ref;U_ref;r_ref;reshape(table',4*N,1)];
%         if first_run
%             args.x0=[X_ref;U_ref;r_ref];
%             first_run=0;
%         else
%             args.x0=[repmat(X_out(:,iter+1),N+1,1);U_ref;r_ref];
%         end
        args.x0=[X_ref;U_ref;r_ref];
        %%	决策变量上下界面 args
        %%  状态上边界
        tempub=[Body.m*world.g*world.mu*6; Body.m*world.g*world.mu*6 ;1000];
        args.ubx=[];
        UBx=[pi*3*ones(3,1);X_out(4:5,iter+1)+[5;5];1;...
            pi*3*ones(3,1);40*ones(3,1)];
        UBx=repmat(UBx,N+1,1);
        UBf=[tempub;tempub;tempub;tempub];
        UBf=repmat(UBf,N,1);
        UBp=repmat([0.4;0.4;2],4,1);
        UBp=repmat(UBp,N,1);
        args.ubx=[args.ubx;UBx;UBf;UBp];
        %%  状态下边界
        templb=[-Body.m*world.g*world.mu*6; -Body.m*world.g*world.mu*6 ;0];
        args.lbx=[];
        LBx=[-pi*3*ones(3,1);X_out(4:5,iter+1)-[5;5];0;...
            -pi*3*ones(3,1);-40*ones(3,1)];
        LBx=repmat(LBx,N+1,1);
        LBf=[templb;templb;templb;templb];
        LBf=repmat(LBf,N,1);
        LBp=repmat([-0.4;-0.4;-2],4,1);
        LBp=repmat(LBp,N,1);
        args.lbx=[args.lbx;LBx;LBf;LBp];
        %%  求解
        sol=solver('x0',args.x0,'lbx', args.lbx,'ubx', args.ubx,'lbg', args.lbg,'ubg', args.ubg,'p',args.p);
        %%  获得结果
        x_li=sol.x(1:n_state*(N+1));
        X_li=reshape(full(x_li),n_state,(N+1));
        
        f_sol=sol.x(n_state*(N+1)+1:n_state*(N+1)+n_F*N);
        f_li=reshape(full(f_sol),n_F,N);
        
        % 预测N个步态分段时长 步态分段时长按iter计数为Gait.iterLength 起时为iter
        % 结束为iter+Gait.iterLength;此用来插值使用;
        %%
        T_t=iter+1:Gait.iterationsBetween:iter+1+Gait.iterLength-Gait.iterationsBetween;
        Tt=iter+1:1:iter+1+Gait.iterLength-Gait.iterationsBetween;
        %%
%         for leg=1:4
%             tableleg=table(:,leg);
%             tableleg.*T_t
%             interp1(T_t,f_li(k,:),Tt)
%         end
%         F_li=[];        
%         for k=1:12
%             F_li=[F_li;interp1(T_t,f_li(k,:),Tt,'pchip')];
%         end
%         for k=1:size(F_li,2)
%             for leg =1:4
%                 xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
%                 p_foot=Body.hipPos(xyz_idx);
%                 feetforce_used=F_li(xyz_idx,k);
%                 plot3([p_foot(1),p_foot(1)+0.01*feetforce_used(1)],...
%                     [p_foot(2),p_foot(2)+0.01*feetforce_used(2)],...
%                     [p_foot(3),p_foot(3)+0.01*feetforce_used(3)]);
%             end
%             hold on;grid on; axis equal;
%         end
        %%
        r_sol=sol.x(n_state*(N+1)+n_F*N+1:n_state*(N+1)+n_F*N+n_r*N);
        r_li=reshape(full(r_sol),n_F,N);
        
        p_li=r_li+repmat(X_li(4:6,1:end-1),4,1);
        %更新落点
        for leg=1:4
            xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
            if table(Gait.nIterations-2,leg)==1
                footSwingTrajectories(leg).endPoint=p_li(xyz_idx,Gait.nIterations-2);
            end
        end
        F_use=f_li(:,2);
        
        % 绘制规划结果
%         for i=1:N
%             cube_animate(X_li(:,i),i,p_li(:,i),~table(i,:),[0;0;0;0],...
%                 f_li(:,i),3,[],[],[],[],[],[-30,30],dt_val,[]);
%             pause(0.3);
%         end
        
        %%
    end
    %% 测试
        F_li=[];
        for k=1:12
            F_li=[F_li;interp1(T_t,f_li(k,:),iter+1,'previous')];
        end
        F_use= F_li;
    %%  腿部控制
    for leg=1:4
        if swingstate(leg)>0 && pre_contact(leg)~=1%%摆动腿控制
            if firstSwing(leg) == 1
                firstSwing(leg) = 0;
                %即刚从支撑切换到摆动 起始点为原点
                footSwingTrajectories(leg).startPoint=pFoot(:,leg);
            end
            F_use((leg-1)*3+1:(leg-1)*3+3)=[0;0;0]; %摆动腿力设零
            %% 摆动腿控制
            % 摆动结束时认为腿到达终点
            if swingstate(leg)==1
                pFoot(:,leg)=footSwingTrajectories(leg).endPoint;
            end
        else%支撑腿控制
            firstSwing(leg) = true;
            % 将计算的对应腿的反力使用
        end
    end
    
    %%  仿真计算
%     F_use=repmat([0,0,9*9.8/4]',4,1);
%     rin=Body.hipPos-repmat([0,0,0.28]',4,1);


    rin=reshape(pFoot,[],1)-repmat(X_out(4:6,iter+1),4,1);
%     k1 = f(X_out(:,iter+1),F_use,rin);   % new
%     X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t*k1);% 欧拉
    
    k1 = fa(X_out(:,iter+1),F_use,rin);   % new
    k2 = fa(X_out(:,iter+1) + delta_t/2*k1,F_use,rin); % new
    k3 = fa(X_out(:,iter+1) + delta_t/2*k2,F_use,rin); % new
    k4 = fa(X_out(:,iter+1) + delta_t*k3,F_use,rin); % new
%     X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t*k1);% 欧拉
    X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t/6*(k1+2*k2+2*k3+k4)); % new runge kutta method
    


    %%  动画
    if rem(iter,IterMPCRate*2)==0
        cube_animate(X_out(:,iter+1),iter+1,pFoot,swingstate,[0;0;0;0],...
            F_use,2,[],[],x_li,[],[],[-30,40],delta_t,footSwingTrajectories,1);%footSwingTrajectories
%%保存gif
                frame = getframe(figure(2));
    [A,map]=rgb2ind(frame2im(frame),256);
    if pic_num==1
        imwrite(A,map,time,'gif','LoopCount',Inf,'DelayTime',delta_t*IterMPCRate*2);
    else
        imwrite(A,map,time,'gif','WriteMode','append','DelayTime',delta_t*IterMPCRate*2);
    end
    pic_num=pic_num+1;
    end
end








%% 步态函数
%% 步态相位
function [s,p]=setGaitIter(Gait,currentIter)
s=rem(fix(currentIter/Gait.iterationsBetween),Gait.nIterations);
p=rem(currentIter,Gait.iterationsBetween*Gait.nIterations)/(Gait.iterationsBetween*Gait.nIterations);
end
%% 接触状态
function p=getConectState(Gait)

p=Gait.phsae-Gait.offsetFloat;
for leg=1:4
    if p(leg)<0
        p(leg)=p(leg)+1;
    end
    if p(leg)>Gait.durationsFloat(leg)
        p(leg)=0;
    else
        p(leg)=p(leg)/Gait.durationsFloat(leg);
    end
    
end
% p(1)=0;%瘸腿
end
%% 摆动状态
function p=getSwingState(Gait)
swing_offset=Gait.offsetFloat+Gait.durationsFloat;
for leg=1:4
    if swing_offset(leg)>1
        swing_offset(leg)=swing_offset(leg)-1;
    end
end
swing_duration=1-Gait.durationsFloat;
p=Gait.phsae-swing_offset;
for leg=1:4
    if p(leg)<0
        p(leg)=p(leg)+1;
    end
    if p(leg)>swing_duration(leg)
        p(leg)=0;
    else
        p(leg)=p(leg)/swing_duration(leg);
    end
    
end
% p(1)=0.5;%瘸腿
end
%% 生成接触状态列表将接触检测加入
%将接触检测加入
function [table,ContactPhase]=getMPCTable(Gait,Iterlen,pre_contact)
%生成从当前时刻向后一段时间的接触状态
%输入：
%     durations:支撑时间
%     offsets:相位差
%     currentIter:当前步态所处片段
%     Iterlen：要产生多少个片段
%     gaitIterlen:步态片段数
%     pre_contact:接触状态
%输出：
%     table:mpc接触预测
table=[];
ContactPhase=[];
CP=zeros(4,1);
switchSwing=zeros(1,4);%提前接触标志
for i=0:Iterlen-1
    iter_in=rem((i+Gait.gait_iteration+1),Gait.nIterations);
    p=iter_in-Gait.offsets;%每条腿在此刻向后预测时的相位 如果此刻接触产生
    for leg=1:4
        if p(leg)<0
            p(leg)=p(leg)+Gait.nIterations;
        end
        if p(leg)>Gait.durations(leg)
            if pre_contact(leg)==1
                p(leg)=1;
                if switchSwing(leg)==0
                    switchSwing(leg)=1;
                end
            else
                p(leg)=0;
            end
            if  switchSwing(leg)==2
                p(leg)=0;
            end
        else
            CP(1)=p(leg)/Gait.durations(leg);
            p(leg)=1;
            if switchSwing(leg)==1
                switchSwing(leg)=2;
            end
        end
    end
%     p(1)=0;%瘸腿
    table=[table;p];
    ContactPhase=[ContactPhase;CP];
end
end

%% 工具函数
function rotxm=rotx(theta)
s=sin(theta);
c=cos(theta);
rotxm=[1,0,0;
    0,c,-s
    0,s c];
end
function rotym=roty(theta)
s=sin(theta);
c=cos(theta);
rotym =[c,0,s;
    0,1,0;
    -s,0,c];
end
function rotzm=rotz(theta)
s=sin(theta);
c=cos(theta);
rotzm=[c,-s,0;
    s,c,0;
    0,0,1];
end
function R=rotsb(theta)
% R=rotx(theta(1))*roty(theta(2))*rotz(theta(3));
R=rotz(theta(3))*roty(theta(2))*rotx(theta(1));
end
function s=Skew(in)
s = [0 -in(3) in(2);
    in(3) 0 -in(1);
    -in(2) in(1) 0];
end

