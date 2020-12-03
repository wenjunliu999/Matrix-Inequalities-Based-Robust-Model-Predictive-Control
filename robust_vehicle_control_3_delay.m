%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Wenjun Liu
% Date: 03/12/2020
%
% Matrix Inequalities Based Robust Model Predictive Control for Vehicle Considering
% Model Uncertainties, External Disturbances and Time-varying Delay (bounded in 3s)
%
% Installation package to be installed---yalmip,penlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
yalmip('clear');
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
[A_temp,B,C,~]=c2dm(Ac,Bc,Cc,Dc,dt);
% model data (time delay considered)

alpha = 0.8; % The limits 1 and 0 correspond to no delay term and to a completed delay term,respectively.

A = alpha*A_temp;
A_d = (1-alpha)*A_temp;
% model data

E  = [0.01;0.1];
NA = 0.05*A;
NB = 0.05*B;
NAd = 0.05*A_d;
nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q   = 5*eye(2);
R   = 1;
%tao = 1;
delay_bound = 3;
N   = 499+delay_bound;% iteration times

% PRI parameter
gamma = 0.8;
gamma_d = 1-gamma;
%
d_s = delay_bound-1;

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
x_d = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
x_r = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

load('phi_ref.mat');
load('Beta.mat');


% Initial state
for i = 1:1:delay_bound
    x{i} = [0; 0];
    x_d{i} = [Beta(i+500);phi_ref(i+500)];
    x_r{i} = x_d{i} + x{i};
end

% Controller
LMI1 = cell(9,9);
LMI2 = cell(3,3);
LMI3 = cell(13,13);
LMI4 = cell(9,9);

lambda = sdpvar(1);
X = sdpvar(2,2);
X_d = sdpvar(2,2);
Y = sdpvar(1,2);
Z    = sdpvar;
M    = eye(2);

xi     = sdpvar(1);
eta     = sdpvar(1);
eta1     = sdpvar(1);
epsilon = sdpvar(1);
sigma = sdpvar(1);
rou  = 0.001;
tao = [1 0.6 0.2];


 % RPI
LMI1 = [gamma*(-1+lambda)*X    zeros(2,2)             zeros(2,1)        (A*X+B*Y)'  (NA*X+NB*Y)'     sigma*M;
        zeros(2,2)             gamma_d*(-1+lambda)*X  zeros(2,1)        (A_d*X)'    (NAd*X)'         zeros(2,2);
        zeros(1,2)             zeros(1,2)             -lambda/(rou^2)    E'          zeros(1,2)      zeros(1,2);
        A*X+B*Y                A_d*X                  E                  -X          zeros(2,2)      zeros(2,2)
        NA*X+NB*Y              NAd*X                  zeros(2,1)         zeros(2,2)  -sigma*eye(2)   zeros(2,2);
        sigma*M'               zeros(2,2)             zeros(2,1)         zeros(2,2)  zeros(2,2)      -sigma*eye(2);];             
          
F =[];          
F = X > 0;
F =[F, 0<lambda<1];
F = F + [LMI1 <= 0];

[model,recoverdata,diagnostic,interfacedata] = export(F,lambda,sdpsettings('solver','penbmi'),[],[],1);
bmi=yalmip2bmi(model);
penm = bmi_define(bmi);
prob = penlab(penm);
solve(prob);
lambda = prob.x(1)
 

 %%
 % RPI
LMI1 = [gamma*(-1+lambda)*X    zeros(2,2)             zeros(2,1)        (A*X+B*Y)'  (NA*X+NB*Y)'     sigma*M;
        zeros(2,2)             gamma_d*(-1+lambda)*X  zeros(2,1)        (A_d*X)'    (NAd*X)'         zeros(2,2);
        zeros(1,2)             zeros(1,2)             -lambda/(rou^2)    E'          zeros(1,2)      zeros(1,2);
        A*X+B*Y                A_d*X                  E                  -X          zeros(2,2)      zeros(2,2)
        NA*X+NB*Y              NAd*X                  zeros(2,1)         zeros(2,2)  -sigma*eye(2)   zeros(2,2);
        sigma*M'               zeros(2,2)             zeros(2,1)         zeros(2,2)  zeros(2,2)      -sigma*eye(2);];              

 LMI2 = [Z   Y;
         Y'  X];                                        
 % guaranted upper bound
 LMI3 = [-X          zeros(2,2)   zeros(2,1)   (A*X+B*Y)'  X                   (Q*X)'       (R*Y)'      (NA*X+NB*Y)'     zeros(2,2);
         zeros(2,2)  -X_d         zeros(2,1)   (A_d*X_d)'  zeros(2,2)          zeros(2,2)   zeros(2,1)  (NAd*X_d)'       zeros(2,2);
         zeros(1,2)  zeros(1,2)   -1*tao(1)*xi (xi*E)'     zeros(1,2)          zeros(1,2)   0            zeros(1,2)      zeros(1,2);
         A*X+B*Y     A_d*X_d      xi*E          -X         zeros(2,2)          zeros(2,2)   zeros(2,1)   zeros(2,2)      eta*M;
         X           zeros(2,2)   zeros(2,1)    zeros(2,2) -1*(d_s+1)^(-1)*X_d zeros(2,2)   zeros(2,1)   zeros(2,2)      zeros(2,2);
         Q*X         zeros(2,2)   zeros(2,1)    zeros(2,2) zeros(2,2)          -xi*Q        zeros(2,1)   zeros(2,2)      zeros(2,2);
         R*Y         zeros(1,2)   zeros(1,1)    zeros(1,2) zeros(1,2)          zeros(1,2)   -xi*R        zeros(1,2)      zeros(1,2);
         NA*X+NB*Y   NAd*X_d      zeros(2,1)    zeros(2,2) zeros(2,2)          zeros(2,2)   zeros(2,1)   -eta*eye(2)     zeros(2,2);
         zeros(2,2)  zeros(2,2)   zeros(2,1)    eta*M'     zeros(2,2)          zeros(2,2)   zeros(2,1)   zeros(2,2)      -eta*eye(2) ];
f = [];
f = f + [LMI1 <= 0] + [LMI2 >= 0]+ [LMI3 < 0];
% f = f + [LMI2 >= 0]+ [LMI3 < 0]; %without RPI
 

%data save
u_record = [];
x_r_record = [];
x_d_record = [];
d_record = [];
x_d_record = [x_d_record value(x_d{delay_bound})];
x_r_record = [x_r_record value(x{delay_bound})+value(x_d{delay_bound})];
for k = delay_bound+1:N   
    d  = unidrnd(delay_bound);
    d_record = [d_record d];
    p  = 0.001*sin(k);
    H = sin(k);
    A = A + M*H*NA;
    B = B + M*H*NB;
    A_d = A_d + M*H*NAd;
    xi_3 = [x{k}',x{k-1}',x{k-2}',x{k-3}']';
    big_gamma = [X          zeros(2,2) zeros(2,2)  zeros(2,2);
                 zeros(2,2) X_d/3      zeros(2,2)  zeros(2,2);
                 zeros(2,2) zeros(2,2) X_d/2       zeros(2,2);
                 zeros(2,2) zeros(2,2) zeros(2,2)  X_d];
%     
    LMI4 = [ones(1)  xi_3';
            xi_3   big_gamma];

    f = f + [LMI4 >= 0]; 
    f = [f,umin <= u{k} <= umax];
    ops = sdpsettings('warning',1,'verbose',1,'solver','sedumi','cachesolvers',1);
    obj = xi;
    
    sol = optimize(f,obj);
    if sol.problem == 0
        disp('Solver thinks it is feasible');
    end
    
     X_ = value(X);
     Y_ = value(Y);
     u{k} = Y_*inv(X_)*x{k};
     x{k+1} = A*x{k}+ A_d*x{k-d}+ B*u{k} + E*p; 
     
     x_d{k+1} = [Beta(k+500); phi_ref(k+500)];
     x_r{k+1} = x_d{k+1} + x{k+1};
     
     u_record = [u_record value(u{k})];
     x_d_record = [x_d_record value(x_d{k+1})];
     x_r_record = [x_r_record value(x_r{k+1})];
  
end

%%


figure
plot(x_d_record(1,:),'linewidth',3);
hold on
plot(x_r_record(1,:),'r--','linewidth',3);
legend( 'Desired','Robust MPC');

figure
plot(x_d_record(2,:),'linewidth',2);
hold on
plot(x_r_record(2,:),'r--','linewidth',2);
legend( 'Desired','Robust MPC');

length_evaluate = length(x_d_record(1,:));
Err_beta = sqrt(norm(x_r_record(1,1:length_evaluate)-x_d_record(1,1:length_evaluate))^2 / length_evaluate);
disp(['beta tracking errors  = ', num2str(Err_beta)])
Err_r = sqrt(norm(x_r_record(2,1:length_evaluate)-x_d_record(2,1:length_evaluate))^2 / length_evaluate);
disp(['r tracking  = ', num2str(Err_r)])
