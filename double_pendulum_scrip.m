syms t l1 l2 theta1 theta2 omega1 omega2 tau1 tau2 m1 m2 M1 M2 alpha1 alpha2 T2x T2y T1x T1y
%% 
% Centre of mass:

%centre of mass of the 2 pendulum-rod systems
L1 = (m1*l1 + M1*l1/2)/(m1 + M1);
L2 = (m2*l2 + M2*l2/2)/(m2 + M2);
I1 = m1*l1^2 + (M1*l1^2)/3;
I2 = m2*((l2/4)^2) + M2*(((l2)^2/12+(l2/4)^2));
%% 
% Defining positions:

x1 = (L1)*sin(theta1);
y1 = -(L1)*cos(theta1);
x2 = (L1)*sin(theta1) + (L2)*sin(theta2);
y2 = (L1)*cos(theta1) - (L2)*cos(theta2);
%% 
% Defining velocities:
% 
% $\omega_1 =\overset{\ldotp }{\theta_1 }$,  $\omega_2 =\overset{\ldotp }{\theta_2 
% }$

%velocities
vx1 = omega1*L1*cos(theta1);
vy1 = omega1*L1*sin(theta1);
vx2 = omega1*L1*cos(theta1)+omega2*L2*cos(theta2);
vy2 = omega1*L1*sin(theta1)+omega2*L2*sin(theta2);
%% 
% Defining accelerations:
% 
% $\alpha_1 =\overset{\ldotp \ldotp }{\theta_1 }$,  $\alpha_2 =\overset{\ldotp 
% \ldotp }{\theta_2 }$

%acceleration
ax1 = L1*(alpha1*cos(theta1) - omega1^2*sin(theta1));
ay1 = L1*(alpha1*sin(theta1) - omega1^2*cos(theta1));
ax2 = L1*(alpha1*cos(theta1) - omega1^2*sin(theta1)) + ...
    L2*(alpha2*cos(theta2) - omega1^2*sin(theta2));
ay2 = L1*(alpha1*sin(theta1) + omega1^2*cos(theta1))+...
    L2*(alpha2*sin(theta2) + omega1^2*cos(theta2));
%% 
% Force and torque equations:

eqa = T2x == (m2+M2)*ax2;
eqb = T2y == (m2+M2)*ay2 + (m2+M2)*9.8;
eqc = T1x - T2x == (m2+M2)*ax1;
eqd = T1y - T2y == (m2+M2)*ay1 + (M1+m1)*9.8;
eqe = I1*alpha1 == tau1 - l1*9.8*sin(theta1)*(m1+M1/2) - T2x*l1*cos(theta1) - T2y*l1*sin(theta1);
eqf = I2*alpha2 == tau2 - T2x*L2*cos(theta2) - T2y*L2*sin(theta2);
%% 
% Using the above equations to get the expressions of $\overset{\ldotp \ldotp 
% }{\theta_1 } \;\textrm{and}\;\overset{\ldotp \ldotp }{\theta_2 }$:

eqn1 = alpha2 == (1/(I2+L2^2*(M2+m2)))*(tau2 - L2*(M2+m2)*9.8*sin(theta2)-...
    L1*L2*(M2+m2)*(alpha1*cos(theta2-theta1)+omega1^2*sin(theta2-theta1)));
eqn2 = alpha1 == (1/(I1+L1*l1*(M2+m2)))*(tau1 - l1*9.8*sin(theta1)*(m1+(M1/2)+m2+M2)-...
    l1*(M2+m2)*(L2*alpha2*cos(theta2-theta1)+L2*omega2^2*sin(theta1-theta2)));
eqns = [eqn1,eqn2];
%% 
% Using $X=A^{-1} B:$

Mat1 = [L1*L2*(M2+m2)*cos(theta2-theta1) I2+L2^2*(M2+m2);I1+L1*l1*(M2+m2) l1*(M2+m2)*L2*cos(theta2-theta1)];
res = [tau2-L2*(M2+m2)*9.8*sin(theta2)-L1*L2*(M2+m2)*(omega1^2*sin(theta2-theta1));...
    tau1 - l1*9.8*sin(theta1)*(m1+(M1/2)+m2+M2)-l1*(M2+m2)*L2*omega2^2*sin(theta1-theta2)];
alphas = inv(Mat1)*res;
% substitute some constants 
out = subs(alphas,[l1 l2 m1 m2 M1 M2],[1 1 1 1 1 1]);
out2 = vpa(simplify(out));


% linearization
A = [0 0 1 0;...
    0 0 0 1;...
    diff(out(1),theta1) diff(out(1),(theta2)) diff(out(1),omega1) diff(out(1),(omega2));...
    diff(out(2),theta1) diff(out(2),(theta2)) diff(out(2),omega1) diff(out(2),(omega2))];
B = [0 0;...
    0 0;...
    diff(out(1),tau1) diff(out(1),tau2);...
    diff(out(2),tau2) diff(out(2),tau2)];

C_new = [1 0 0 0;0 1 0 0];
D_new = zeros(2,2);

xx_s = [theta1; theta2; omega1; omega2];
xxDot_s = [omega1, omega2, alpha1, alpha2];
uu_s = [tau1;tau2];
%linear state-space model
xxDotExpr_s = A*xx_s + B*uu_s;

%simulation
%xxFunc_s = matlabFunction(xxDotExpr_s, 'vars', {t, xx_s, uu_s});
%ode45(@(t,x) xxFunc_s(t, x, [0;0]), 0:1e-2:10, [15;-10; zeros(2,1)]/180*pi);
%legend ('theta1','theta2','omega1','omega2')
%xlabel('time(0-10s)') 
%ylabel('angular positions and velocities')

A_new = double(subs(A,[theta1; theta2; omega1; omega2; tau1;tau2],...
    [0;0;0;0;0;0]));
B_new = double(subs(B,[theta1; theta2; omega1; omega2; tau1;tau2],...
    [0;0;0;0;0;0]));

%% plotting state space with initial condition

ssobj = ss(A_new,B_new,C_new,D_new);
tf_obj = tf(ssobj);

%% Define Parameters
h=1/2;            %time step
n=500;           %no of steps
ySP = [pi/3;pi/2];    %set-point
m=4;            % Control horizon
p=10;           % Prediction horizon
Q=[1 0;0 1];    % Output weight
R=[1 0;0 1];         % Input weight
duMax=[0.2; 0.2];     % Input rate constraint
duMin=-duMax;
uMin=[-50;-50];      % Input constraint
uMax=[50;50];

dmc_MIMO(tf_obj,n,h,ySP,m,p,Q,R,duMin,duMax,uMin,uMax)
%%
%unto simulink
set_param('simulink_model/State-Space','A','A_new','B','B_new','C','C_new','D','D_new') 
%%

%decouple
G = tf_obj;
G.InputName = {'Tau1','Tau2'};
G.OutputName = 'y';
DC = tunableGain('Decoupler',eye(2));
DC.InputName = 'e';
DC.OutputName = {'ta1','ta2'};

PID_1 = tunablePID('PID_1','pid');
PID_1.InputName = 'ta1';
PID_1.OutputName = 'Tau1';
  
PID_2 = tunablePID('PID_2','pid'); 
PID_2.InputName = 'ta2';
PID_2.OutputName = 'Tau2'; 

sum1 = sumblk('e = r - y',2);
C0 = connect(PID_1,PID_2,DC,sum1,{'r','y'},{'Tau1','Tau2'});
wc = [0.1,100];
[G,C,gam,Info] = looptune(G,C0,wc);
showTunable(C)
T = connect(G,C,'r','y');
step(T)
figure('Position',[100,100,520,1000])
loopview(G,C,Info)