%% �������Է�����ģ��Ԥ�����
clc;
clear;
warning off
addpath('E:\Matlab\casadi-matlab')
import casadi.*
%%
mpcHz=60;%mpcƵ��
delta_t = 0.002;%����ʱ����
IterMPCRate=fix(1/mpcHz/delta_t);%mpcƵ�ʿ���
N=8;
runtime=14;
%% %�������

Body.m = 9;%����������
%�������
Body.Ib = diag([0.07,0.26,0.242]);%ת����������
Body.length_body=0.38;
Body.width_body=0.22;
Body.hipPos=diag([1 1 1, 1 -1 1, -1 1 1, -1 -1 1])*repmat([0.2 0.1 0],1,4)';
world.g = 9.8;%�������ٶ�
world.mu=0.4;%Ħ��ϵ��

%%	��̬����
%һ����̬Ƭ��ʱ��
Gait.dt=0.05;
%һ����̬Ƭ������ѭ���е�iter��
Gait.iterationsBetween=fix(Gait.dt/delta_t);
dt_val = repmat(Gait.dt,1,N)';
%��̬��λ
Gait.phsae=0;
%��̬Ƭ����
Gait.nIterations=8;
Gait.Tp=Gait.nIterations*Gait.dt;
Gait.iterLength=Gait.iterationsBetween*Gait.nIterations;
%��̬��λ��Ƭ�μ�����ʽ
% offsets=[3,0,8,5]*5;
Gait.offsets=[0,4,4,0];
%��̬��λ��С����ʽ
Gait.offsetFloat=Gait.offsets/Gait.nIterations;
%ÿ����֧�ų���ʱ��ʱ��Ƭ�μ���
Gait.durations=ones(1,4)*Gait.nIterations*0.5;
%ÿ����֧�ų���ʱ����ѭ��ʱ��Ƭ�μ���
Gait.durationsFloat=Gait.durations/Gait.nIterations;
%��̬Ƭ�μ���ʱ��Ƭ�μ���
Gait.gait_iteration=0;
%֧��ʱ��Ƭ�μ���
Gait.standtime=Gait.durations(1);
%�ڶ�ʱ��Ƭ�μ���
Gait.swingtime=Gait.nIterations-Gait.standtime;
%֧��ʱ�� s
Gait.standTime=Gait.standtime*Gait.dt;%��̬�ֶ�ʱ��
%�ڶ�ʱ�� s
Gait.swingTime=Gait.swingtime*Gait.dt;
%%	����������Ҫ

firstSwing=ones(1,4);%��ʼ�ڶ���־
swingTimeRemaining=ones(4,1)*Gait.swingTime;%ʣ��ڶ�ʱ��
side_signy=[1,-1,1,-1];%��ŵ�ƫִ����
side_signx=[1,1,-1,-1];%��ŵ�ƫִ����
pre_contact=[0,0,0,0];%��ǰ�Ӵ���־
footRef=Body.hipPos;
% �ڶ��켣��ʼ��
footSwingTrajectories(1:4)=struct;
for i=1:4
    footSwingTrajectories(i).startPoint=zeros(3,1);
    footSwingTrajectories(i).endPoint=zeros(3,1);
end

%% %״̬��Ȩ��

weight.QX = [5000 1000 2500, 0 0 5000, 2.5 2.5 300, 2500 2500 2500 ]';
weight.QN = [5000 1000 2500, 0 0 5000, 2.5 2.5 300, 2500 2500 2500 ]';
weight.Qc = [350 350 350]';
weight.Qf = [0.25 0.25 0.0001]';
%% ����΢�ַ���
X_=SX.sym('X_', 12, 1);
n_state=size(X_,1);
F_=SX.sym('F_', 12, 1);
n_F=size(F_,1);
R_=SX.sym('R_', 12, 1);
n_r=size(R_,1);
%% ����΢�ַ���
I3=eye(3);
Rbody=rotsb(X_(1:3));

cy = cos(X_(3));
sy = sin(X_(3));
cp = cos(X_(2));
sp = sin(X_(2));

R_yaw =[cy sy 0;
    -sy cy 0;
    0 0 1];%���絽����
R_w=[cy/cp,sy/cp,0;
    -sy,cy,0;
    cy*sp/cp,sy*sp/cp,1];
Ig = Rbody*Body.Ib*Rbody';
Ig_inv=Ig\I3;
%hip������ϵ

A = [zeros(3) zeros(3) R_yaw zeros(3)  ;
    zeros(3) zeros(3) zeros(3) I3 ;
    zeros(3) zeros(3) zeros(3) zeros(3);
    zeros(3) zeros(3) zeros(3) zeros(3) ;
    ];%״̬����
AA=A;
AA(1:3,7:9)=R_w;
B=[zeros(3)           zeros(3)           zeros(3)            zeros(3);
    zeros(3)           zeros(3)           zeros(3)            zeros(3);
    Ig_inv*Skew(R_(1:3)) Ig_inv*Skew(R_(4:6)) Ig_inv*Skew(R_(7:9))  Ig_inv*Skew(R_(10:12));
    I3/Body.m   I3/Body.m   I3/Body.m    I3/Body.m;];%���ƾ���
grav=zeros(12,1);
grav(12)=-world.g;

dotX=A*X_+B*F_+grav;
Dotx=AA*X_+B*F_+grav;%
%% ����΢�ַ��̺���
% ��
f=Function('f',{X_;F_;R_},{dotX},{'input_states','control_inputs','foot_input'},{'dotX'});
% ��ȷ
fa=Function('fa',{X_;F_;R_},{Dotx},{'input_states','control_inputs','foot_input'},{'dotX'});

%%  ������ۺ�Լ�� ��������
X = SX.sym('X', n_state, N+1); % N+1��״̬
F = SX.sym('F', n_F, N); % N���ڵĿ���
r = SX.sym('r', n_r, N); % N���ڵĿ���

RefX = SX.sym('RefX', n_state, N+1); % N���ڵĿ������
RefF = SX.sym('RefF', n_F, N); % N���ڵĿ������
Refr = SX.sym('Refr', n_r, N); % N���ڵĿ������

ContactState=SX.sym('ConState', 4, N);
obj=0;
% g=[];
%%  ������ۺ�Լ�� Ħ�����˶�ѧ����
mu_inv = 1.0/world.mu;
%Ħ��Լ��
f_block =[ mu_inv, 0,  -1.0;
    -mu_inv, 0,  -1.0;
    0,  mu_inv, -1.0;
    0, -mu_inv, -1.0;];
%�˶�ѧԼ��
kin_box_x = 0.05;
kin_box_y = 0.05;
kin_box_z = 0.3;

Kin_block =[ 1, 0,  0,-kin_box_x;
    -1, 0,  0,-kin_box_x;
    0, 1, 0,-kin_box_y;
    0, -1, 0,-kin_box_y;
    0, 0, 1,0.06;
    0, 0, -1,-kin_box_z];
%��˻�����
Phip_swing=Body.hipPos;
Phip_swing([3,6,9,12])=[-0.2,-0.2,-0.2,-0.2];
%%  Լ���ݴ�������� %��̬Լ��

defect_init=X(:,1)-RefX(:,1);%12*1 eq

defect_state=SX.zeros(12*(N+1),1);%12(N+1)*1 eq
defect_FootOnGround=SX.zeros(4*(N),1);%4(N)*1 eq
defect_footStance=SX.zeros(12*(N),1);%(3*4)(N)*1 eq
n_equa_c=12+12*(N+1)+4*(N)+12*(N);%12+
%��
defect_legLimits=SX.zeros(24*(N),1);%(4*6)(N)*1
defect_footforce=SX.zeros(16*(N),1);%(4*4)(N)*1 xyĦ��Լ��4��
defect_ForceNormal=SX.zeros(N,1);% (N)*1
defect_footswing=SX.zeros(4*(N),1);%4(N)*1
n_inequa_c=24*(N)+16*(N)+N+4*(N);
%%	Լ���ʹ��ۼ���
for k = 1:N
    %%	���ۼ���
    Xk=X(:,k);
    Fk=F(:,k);
    rk=r(:,k);
    Pk=repmat(Xk(4:6),4,1)+rk;
    ContactStatek=ContactState(:,k);
    dtk=dt_val(k);
    X_err = Xk - RefX(:,k);                                         % ����״̬���
    pf_err = repmat(Xk(4:6),4,1) + Phip_swing - Pk;                      %  ����ʱԼ��footλ��
    r_err=rk-Refr(:,k);
    U_err = Fk - RefF(:,k);                                         % GRF ���
    obj = obj + (X_err'*diag(weight.QX)*X_err+...                     % ������
        r_err'*diag(repmat(weight.Qc,4,1))*r_err+...
        U_err'*diag(repmat(weight.Qf,4,1))*U_err)*dtk;
    %%	Լ������
    %״̬Լ��
    %% runge kutta method
    %     k1 = f(Xk,Fk,Pk);   % new
    %     k2 = f(Xk + dtk/2*k1,Fk,rk); % new
    %     k3 = f(Xk + dtk/2*k2,Fk,rk); % new
    %     k4 = f(Xk + dtk*k3,Fk,rk); % new
    %     st_next_RK4=Xk +dtk/6*(k1+2*k2+2*k3+k4); % new
    %     defect_state((k-1)*12+1:(k-1)*12+12)=X(:,k+1)-(st_next_RK4);
    %%
    defect_state((k-1)*12+1:(k-1)*12+12)=X(:,k+1)-(Xk+f(Xk,Fk,rk)*dtk);
    %����������0 ����ʽ
    defect_ForceNormal(k)=-dot(Fk,repmat([0;0;1],4,1));
    %��Ϸ���������0��Ħ��Լ����Լ���ڶ�����Ϊ0 ������� ����ʽ
    defect_footswing((k-1)*4+1:(k-1)*4+4)=Fk([3,6,9,12])-ContactStatek.*repmat(1000,4,1);
    for leg=1:4
        xyz_idx = 3*(leg-1)+1:3*(leg-1)+3;
        %���ڵ���Լ�� 0�Ǵ�ʱ����߶ȵ�ʽ
        defect_FootOnGround((k-1)*4+leg)=ContactStatek(leg)*Pk(3*(leg-1)+3);
        %�����ȳ� ���Ʒ�Χ����ʽ
        Rbody=rotsb(Xk(1:3));
        Phip=Rbody*Body.hipPos(xyz_idx)+Xk(4:6);
        p_rel = (Pk(xyz_idx) - Phip);%hip->���
        defect_legLimits((k-1)*24+(leg-1)*6+1:(k-1)*24+(leg-1)*6+6)= Kin_block*[p_rel;1];
        %�Ӵ��нŲ�����
        if (k < N)
            Pk1=repmat(X(4:6,k+1),4,1)+r(:,k+1);
            defect_footStance((k-1)*12+(leg-1)*3+1:(k-1)*12+(leg-1)*3+3)=ContactStatek(leg)*(Pk1(xyz_idx)-Pk(xyz_idx));%ContactState(leg,k+1)*
        end
        %Ħ��Լ�� ����ʽ
        defect_footforce((k-1)*16+(leg-1)*4+1:(k-1)*16+(leg-1)*4+4)=f_block*Fk(xyz_idx);
    end
end
%%	Լ������ defect_init;
g=[defect_init;defect_state;defect_FootOnGround;defect_footStance;...
    defect_legLimits;defect_footforce;defect_ForceNormal;defect_footswing];
display_str=['��ʽԼ������',num2str(n_equa_c),'   ����ʽԼ������',num2str(n_inequa_c)];
disp(display_str);
%%	�ն� cost
X_err = X(:,end)-RefX(:,end);
obj = obj + X_err'*diag(weight.QN)*X_err;
%%	����������������
OPT_variables = [reshape(X,n_state*(N+1),1);reshape(F,n_F*N,1);reshape(r,n_r*N,1)];
OPT_Param = [reshape(RefX,n_state*(N+1),1);reshape(RefF,n_F*N,1);reshape(Refr,n_r*N,1);reshape(ContactState,4*N,1)];
nlp_prob =struct('f', obj, 'x',OPT_variables,'p',OPT_Param, 'g',g);
%%  �Ż�����
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

%% ���������
solver = nlpsol('solver', 'ipopt', nlp_prob,opts_setting);
%%	Լ�����½��� args
args.lbg(1:n_equa_c) = 0;  % -1e-20  % Equality constraints
args.ubg(1:n_equa_c) = 0;  % 1e-20   % Equality constraints

args.lbg(n_equa_c+1 : n_equa_c+ n_inequa_c) = -inf; % inequality constraints
args.ubg(n_equa_c+1 : n_equa_c+ n_inequa_c) = 0; % inequality constraints

%%  ��ͼ�þ���
contact_list=[];%���ƽӴ���λͼ
swing_list=[];%���ưڶ���λͼ

%%  �����ٶ�
des_vel=[0.0,0,0]';%�������ٶ�
v_des_world=[0.6,0,1]';%�������ٶ� ���һ������Ϊ1
yaw_turn_rate=0.0;%yaw�ٶ�����
yaw_des=0;%yaw����
onmiMode=true;%�Ƿ���ͷģʽ��1 ʹ�������ٶȿ��ƻ����ˣ�0 ʹ�ñ����ٶȿ���

%%  �������
X_init = [0;0.0;0; 0.0;0.0;0.28 ;0;0;0; 0;0;0];%��ʼ״̬����
X_out(:,1)=X_init;%�������״̬���� �൱�ڹ۲��� �������¼�������������е�״̬
P_foot=reshape(Body.hipPos,[],1);
pFoot=reshape(Body.hipPos,3,4);
pFoot_next=reshape(Body.hipPos,3,4);
swingTimeCount=0;
first_run=1;
pic_num = 1;%����gif��
time=['NLP','_',datestr(datetime('now'),'yyyy-mm-dd-HH-MM'),'_trot.gif'];

%%  �ο��켣����
X_ref=zeros(12*(N+1),1);
U_ref=zeros(12*(N),1);
r_ref=zeros(12*(N),1);


%%  ��ѭ��
for iter = 0:runtime/delta_t+1
    
%       if iter <500
%         v_des_world=[0;0;0];
%         yaw_turn_rate=0;
%       else
%         v_des_world=[0.6;0;0];
%         yaw_turn_rate=0.0;
%       end
    
    %% �ٶ�ת��
    Rbody=rotsb(X_out(1:3,iter+1));
    if onmiMode
        des_vel=Rbody'*v_des_world;
    else
        v_des_world=Rbody*des_vel;
    end
    
    %%  �״�������ʼ��
    if first_run
        for i=1:4
            footSwingTrajectories(i).startPoint=pFoot(:,i);
            footSwingTrajectories(i).endPoint=pFoot(:,i);
        end
%         first_run=0;
    else
        %         %%  �������λ��
        %         for i=1:4
        %             pFoot(i)=pFoot(i);
        %         end
    end
    %% ��̬����
    [Gait.gait_iteration,Gait.phsae]=setGaitIter(Gait,iter+1);%���㲽̬��λ
    contactstate=getConectState(Gait);%����Ӵ���λ
    swingstate=getSwingState(Gait);%����ڶ���λ
    %contact_list=[contact_list;contactstate];%����Ӵ���λ
    swing_list=[swing_list;swingstate];%����ڶ���λ
    %% ʣ��ڶ�ʱ�����
    for leg=1:4
        xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
        if firstSwing(leg)
            %��ʱ�ڽӴ�
            swingTimeRemaining(leg)=Gait.swingTime;
        else
            %��ʱ�ڰڶ�
            swingTimeRemaining(leg)=swingTimeRemaining(leg)-delta_t;
        end
    end
    %%  ���Ƽ��� abs(Gait.phsae-0)<0.0126 ||abs(Gait.phsae-0.5)<0.0126
    if rem(iter,IterMPCRate)==0
        %%	�������� �ο�mit����
        for leg=1:4
            if swingstate(leg)>0 && pre_contact(leg)~=1%ֻ����ڶ���
                xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
                offest_hip=[0.0*side_signx(leg);0.0*side_signy(leg);0];%�����ƫ�� ����ʵ������������
                pRobotFrame=Body.hipPos(xyz_idx)-offest_hip;
                pYawCorrected=rotz(-yaw_turn_rate*Gait.standTime/2)*pRobotFrame;%��ֵΪ�˴ٽ���ת
                Pf=X_out(4:6,1)+Rbody*(pYawCorrected+des_vel*swingTimeRemaining(leg));
                
                pfx_rel = X_out(10,1) * (0.5 + 0) * Gait.standTime +0.03*(X_out(10,1)-v_des_world(1)) +(0.5* X_out(6,1)/9.81) * (X_out(11,1)*yaw_turn_rate);
                pfy_rel = X_out(11,1) * 0.5 * Gait.standTime +0.03*(X_out(11,1)-v_des_world(2)) +(0.5*X_out(6,1)/9.81) * (-X_out(10,1)*yaw_turn_rate);
                pfx_rel = min(max(pfx_rel, -0.3), 0.3);
                pfy_rel = min(max(pfy_rel, -0.3), 0.3);
                
                Pf(1)=Pf(1)+pfx_rel;
                Pf(2)=Pf(2)+pfy_rel;
                Pf(3)=0;
                %���°ڶ��յ�
                %                 footSwingTrajectories(leg).endPoint=Pf;
                footRef(xyz_idx)=Pf;
            end
        end
        %% ��̬ʱ�������
        table=getMPCTable(Gait,N,zeros(4,1));
        %% �򵥼���ο��켣
        yaw_des=yaw_des+yaw_turn_rate*Gait.dt;
        X_ref=X_out(:,iter+1);
        X_to_sum = [ 0;0;yaw_turn_rate*Gait.dt;v_des_world(1)*Gait.dt;v_des_world(2)*Gait.dt;0.0; 0.0;0;0; 0;0;0];%�ο�ֵ�ۼ���
        %%�������һ�������켣 ���������� X_out(4,iter) X_out(5,iter)
        X_init_state=[0;0;yaw_des;X_out(4,iter+1);X_out(5,iter+1);0.28;0;0;yaw_turn_rate;v_des_world(1);v_des_world(2);0];%�ο�ֵ��ֵ
        X_ref(13:24)=X_init_state;
        %�ٶ��ۼƼ���ο��켣
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
        %%  ��� ׼���ο� ��ֵ
        args.p=[X_ref;U_ref;r_ref;reshape(table',4*N,1)];
%         if first_run
%             args.x0=[X_ref;U_ref;r_ref];
%             first_run=0;
%         else
%             args.x0=[repmat(X_out(:,iter+1),N+1,1);U_ref;r_ref];
%         end
        args.x0=[X_ref;U_ref;r_ref];
        %%	���߱������½��� args
        %%  ״̬�ϱ߽�
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
        %%  ״̬�±߽�
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
        %%  ���
        sol=solver('x0',args.x0,'lbx', args.lbx,'ubx', args.ubx,'lbg', args.lbg,'ubg', args.ubg,'p',args.p);
        %%  ��ý��
        x_li=sol.x(1:n_state*(N+1));
        X_li=reshape(full(x_li),n_state,(N+1));
        
        f_sol=sol.x(n_state*(N+1)+1:n_state*(N+1)+n_F*N);
        f_li=reshape(full(f_sol),n_F,N);
        
        % Ԥ��N����̬�ֶ�ʱ�� ��̬�ֶ�ʱ����iter����ΪGait.iterLength ��ʱΪiter
        % ����Ϊiter+Gait.iterLength;��������ֵʹ��;
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
        %�������
        for leg=1:4
            xyz_idx=3*(leg-1)+1:3*(leg-1)+3;
            if table(Gait.nIterations-2,leg)==1
                footSwingTrajectories(leg).endPoint=p_li(xyz_idx,Gait.nIterations-2);
            end
        end
        F_use=f_li(:,2);
        
        % ���ƹ滮���
%         for i=1:N
%             cube_animate(X_li(:,i),i,p_li(:,i),~table(i,:),[0;0;0;0],...
%                 f_li(:,i),3,[],[],[],[],[],[-30,30],dt_val,[]);
%             pause(0.3);
%         end
        
        %%
    end
    %% ����
        F_li=[];
        for k=1:12
            F_li=[F_li;interp1(T_t,f_li(k,:),iter+1,'previous')];
        end
        F_use= F_li;
    %%  �Ȳ�����
    for leg=1:4
        if swingstate(leg)>0 && pre_contact(leg)~=1%%�ڶ��ȿ���
            if firstSwing(leg) == 1
                firstSwing(leg) = 0;
                %���մ�֧���л����ڶ� ��ʼ��Ϊԭ��
                footSwingTrajectories(leg).startPoint=pFoot(:,leg);
            end
            F_use((leg-1)*3+1:(leg-1)*3+3)=[0;0;0]; %�ڶ���������
            %% �ڶ��ȿ���
            % �ڶ�����ʱ��Ϊ�ȵ����յ�
            if swingstate(leg)==1
                pFoot(:,leg)=footSwingTrajectories(leg).endPoint;
            end
        else%֧���ȿ���
            firstSwing(leg) = true;
            % ������Ķ�Ӧ�ȵķ���ʹ��
        end
    end
    
    %%  �������
%     F_use=repmat([0,0,9*9.8/4]',4,1);
%     rin=Body.hipPos-repmat([0,0,0.28]',4,1);


    rin=reshape(pFoot,[],1)-repmat(X_out(4:6,iter+1),4,1);
%     k1 = f(X_out(:,iter+1),F_use,rin);   % new
%     X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t*k1);% ŷ��
    
    k1 = fa(X_out(:,iter+1),F_use,rin);   % new
    k2 = fa(X_out(:,iter+1) + delta_t/2*k1,F_use,rin); % new
    k3 = fa(X_out(:,iter+1) + delta_t/2*k2,F_use,rin); % new
    k4 = fa(X_out(:,iter+1) + delta_t*k3,F_use,rin); % new
%     X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t*k1);% ŷ��
    X_out(:,iter+2)=X_out(:,iter+1) +full(delta_t/6*(k1+2*k2+2*k3+k4)); % new runge kutta method
    


    %%  ����
    if rem(iter,IterMPCRate*2)==0
        cube_animate(X_out(:,iter+1),iter+1,pFoot,swingstate,[0;0;0;0],...
            F_use,2,[],[],x_li,[],[],[-30,40],delta_t,footSwingTrajectories,1);%footSwingTrajectories
%%����gif
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








%% ��̬����
%% ��̬��λ
function [s,p]=setGaitIter(Gait,currentIter)
s=rem(fix(currentIter/Gait.iterationsBetween),Gait.nIterations);
p=rem(currentIter,Gait.iterationsBetween*Gait.nIterations)/(Gait.iterationsBetween*Gait.nIterations);
end
%% �Ӵ�״̬
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
% p(1)=0;%ȳ��
end
%% �ڶ�״̬
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
% p(1)=0.5;%ȳ��
end
%% ���ɽӴ�״̬�б��Ӵ�������
%���Ӵ�������
function [table,ContactPhase]=getMPCTable(Gait,Iterlen,pre_contact)
%���ɴӵ�ǰʱ�����һ��ʱ��ĽӴ�״̬
%���룺
%     durations:֧��ʱ��
%     offsets:��λ��
%     currentIter:��ǰ��̬����Ƭ��
%     Iterlen��Ҫ�������ٸ�Ƭ��
%     gaitIterlen:��̬Ƭ����
%     pre_contact:�Ӵ�״̬
%�����
%     table:mpc�Ӵ�Ԥ��
table=[];
ContactPhase=[];
CP=zeros(4,1);
switchSwing=zeros(1,4);%��ǰ�Ӵ���־
for i=0:Iterlen-1
    iter_in=rem((i+Gait.gait_iteration+1),Gait.nIterations);
    p=iter_in-Gait.offsets;%ÿ�����ڴ˿����Ԥ��ʱ����λ ����˿̽Ӵ�����
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
%     p(1)=0;%ȳ��
    table=[table;p];
    ContactPhase=[ContactPhase;CP];
end
end

%% ���ߺ���
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

