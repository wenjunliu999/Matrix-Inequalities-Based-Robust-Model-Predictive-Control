%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Wenjun Liu
% Date: 03/12/2020
%
% Matrix Inequalities Based Robust Model Predictive Control for Vehicle Considering
% Model Uncertainties, External Disturbances
%
% Installation package to be installed---yalmip,penlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

%Vehicle Parameters
vx  =10;  % m/s    [Longitudinal Velocity]
cf =3000; % N/rad  [Front wheel coefficient]
cr =3000; % N/rad  [Rear wheel coefficient]
a1 =1.0;  % m      [Front to CG distance]
a2 =1.6;  % m      [Rear to CG distance]
L  =2.6;  % m      [Wheel Base]
Iz =1650; % Kg.m^2 [Moment of Interia]
m  =1000; % Kg     [Mass]

umax = 30*pi/180; % maximum steering angle
umin =-30*pi/180; % minimum steering angle

%Lateral Control Model: time invariant model fixed longitudinal velocity
Ac =[-(cf+cr)/(m*vx),(-a1*cf+a2*cr)/(m*vx*vx)-1;(-a1*cf+a2*cr)/Iz,-(a1*a1*cf+a2*a2*cr)/(Iz*vx)];
Bc =[cf/(m*vx);a1*cf/Iz];
Cc =[0,1];
Dc = 0;

dt =0.01;% sec
%discretize model
[A,B,C,~]=c2dm(Ac,Bc,Cc,Dc,dt);

% model data

E  = [0.01;0.1];
NA = 0.02*A;
NB = 0.02*B;
nx = 2; % Number of states
nu = 1; % Number of inputs
%ndyn = 2; % Number of polytopic systems

% MPC data
Q   = 5*eye(2);
R   = 1;
%tao = 1;
N   = 499;% iteration times

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
x_d = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
x_r = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

%d = binvar(repmat(2,1,N),repmat(1,1,N));
load('phi_ref.mat');
load('Beta.mat');
load('Beta_mpc.mat');
load('phi_mpc.mat');

% Initial state
x{1} = [0; 0];
x_r{1} = [Beta(1);phi_ref(1)];
x_d{1} = [Beta(1);phi_ref(1)];
% y_r{1} = [phi_ref(1)];
% temp = x_r{1};
% y_d{1} = [phi_ref(1)];

%A  = sdpvar(nx,nx);
%B  = sdpvar(nx,1);

% Controller
LMI1 = cell(9,9);
LMI2 = cell(3,3);
LMI3 = cell(13,13);
LMI4 = cell(3,3);
X = sdpvar(2,2);
G = sdpvar(2,2);
%G    = sdpvar(repmat(2,2,2),repmat(2,1,2));
Y = sdpvar(1,2);
%Y    = sdpvar(repmat(2,2,2),repmat(1,1,2));
Z    = sdpvar;
M    = eye(2);
lambda = sdpvar(1);
xi     = sdpvar(1);
eta     = sdpvar(1);
eta1     = sdpvar(1);
epsilon = sdpvar(1);
gamma  = 0.000001;
tao = [1 0.6 0.2];


 LMI1 = [(-1+lambda)*X  zeros(2,1)     (A*X+B*Y)' (NA*X+NB*Y)' zeros(2,2);
              zeros(1,2)  -lambda/(gamma^2)*eye(1)   E'  zeros(1,2)  zeros(1,2);
              A*X+B*Y         E           -X          zeros(2,2)    (eta1*M')';
              NA*X+NB*Y    zeros(2,1)       zeros(2,2) -eta1*eye(2) zeros(2,2);
              zeros(2,2)   zeros(2,1)      eta1*M'      zeros(2,2)  -eta1*eye(2)];              
          
F =[];          
F = X > 0;
F =[F, 0<lambda<1];
F =[F,Z <= umax*umax];
F = F + [LMI1 <= 0];

[model,recoverdata,diagnostic,interfacedata] = export(F,lambda,sdpsettings('solver','penbmi'),[],[],1);
bmi=yalmip2bmi(model);
penm = bmi_define(bmi);
prob = penlab(penm);
solve(prob);
lambda = prob.x(7)
 

 %%
 LMI1 = [(-1+lambda)*X  zeros(2,1)     (A*X+B*Y)' (NA*X+NB*Y)' zeros(2,2);
              zeros(1,2)  -lambda/(gamma^2)*eye(1)   E'  zeros(1,2)  zeros(1,2);
              A*X+B*Y         E           -X          zeros(2,2)    (eta1*M')';
              NA*X+NB*Y    zeros(2,1)       zeros(2,2) -eta1*eye(2) zeros(2,2);
              zeros(2,2)   zeros(2,1)      eta1*M'      zeros(2,2)  -eta1*eye(2)];          

 LMI2 = [Z   Y;
         Y'  X];                                        %(3.61)
 LMI3 = [-X         zeros(2,2)      (A*X+B*Y)'    (Q*X)'       (R*Y)'      (NA*X+NB*Y)'  zeros(2,2);
        zeros(2,2)  -1*tao(1)*xi*eye(2) (xi*eye(2))'  zeros(2,1)  zeros(2,2) zeros(2,2)  zeros(2,2)
         A*X+B*Y    xi*eye(2)       -X     zeros(2,2)    zeros(2,1)    zeros(2,2)   (epsilon*M')';
          Q*X      zeros(2,2)       zeros(2,2)  -xi*Q  zeros(2,1)    zeros(2,2)    zeros(2,2);
          R*Y     zeros(1,2)        zeros(1,2)   zeros(1,2)  -xi*R   zeros(1,2)    zeros(1,2);
          NA*X+NB*Y  zeros(2,2)     zeros(2,2) zeros(2,2) zeros(2,1)   -epsilon*eye(2) zeros(2,2);
          zeros(2,2)   zeros(2,2)   epsilon*M' zeros(2,2) zeros(2,1)   zeros(2,2)    -epsilon*eye(2) ];

f = [];
f = f + [LMI1 <= 0] + [LMI2 >= 0]+ [LMI3 < 0];
    
for k = 1:N   
    d  = 0.0000001*sin(k);
    H = sin(k);
    A = A + M*H*NA;
    B = B + M*H*NB;
    LMI4 = [ones(1)  x{k}';   
               x{k}       X];
    f = f + [LMI4 >= 0]; 
    f = [f,umin <= u{k} <= umax];%+ [ H'*H<= eye];
    ops = sdpsettings('warning',1,'verbose',1,'solver','sedumi','cachesolvers',1);
    obj = xi;
    %optimize(F,[],ops);
    sol = optimize(f,obj);
    if sol.problem == 0
        disp('Solver thinks it is feasible');
    end
    
    %X = cell2mat(X);
     X_ = value(X);
     Y_ = value(Y);
     u{k} = Y_*inv(X_)*x{k};
     x{k+1} = A*x{k}+ B*u{k} + E*d; 
     %x{k} = x_d{k} - x_r{k};
     x_d{k+1} = [Beta(k); phi_ref(k)];
     x_r{k+1} = x_d{k+1} + x{k+1}
end

%%
x_ = cell2mat(x_r(1:end));
u_ = cell2mat(u);
x = cell2mat(x);
x_d = cell2mat(x_d(1:end));

figure
plot(x_d(1,:),'linewidth',2);
hold on
plot(x_(1,:),'r--','linewidth',2);
hold on
plot(Beta_mpc(1:100),'k','linewidth',2);
legend( 'Desired','Robust MPC','MPC');
figure
plot(x_d(2,:),'linewidth',2);
hold on
plot(x_(2,:),'r--','linewidth',2);
hold on
plot(phi_mpc(1:100),'k','linewidth',2);
legend( 'Desired','Robust MPC','MPC');

length_evaluate = length(x_d(1,:));
Err_beta = sqrt(norm(x_(1,1:length_evaluate)-x_d(1,1:length_evaluate))^2 / length_evaluate);
disp(['beta tracking errors  = ', num2str(Err_beta)])
Err_r = sqrt(norm(x_(2,1:length_evaluate)-x_d(2,1:length_evaluate))^2 / length_evaluate);
disp(['beta tracking  = ', num2str(Err_r)])

