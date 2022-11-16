%% Newtonian Method
syms t l1 l2 theta1 theta2 omega1 omega2 alpha1 alpha2 m1 m2 tau1 tau2 g
%position
x1 = (l1)*sin(theta1);
y1 = -(l1)*cos(theta1);
x2 = (l1)*sin(theta1) + (l2)*sin(theta2);
y2 = -(l1)*cos(theta1) - (l2)*cos(theta2);
%vel
vx1 = omega1*l1*cos(theta1);
vy1 = omega1*l1*sin(theta1);
vx2 = omega1*l1*cos(theta1)+omega2*l2*cos(theta2);
vy2 = omega1*l1*sin(theta1)+omega2*l2*sin(theta2);
%acc
ax1 = l1*(alpha1*cos(theta1) - omega1^2*sin(theta1));
ay1 = l1*(alpha1*sin(theta1) + omega1^2*cos(theta1));
ax2 = l1*(alpha1*cos(theta1) - omega1^2*sin(theta1)) + ...
    l2*(alpha2*cos(theta2) - omega1^2*sin(theta2));
ay2 = l1*(alpha1*sin(theta1) + omega1^2*cos(theta1))+...
    l2*(alpha2*sin(theta2) + omega1^2*cos(theta2));
%alpha
Mat1 = [l1^2*(m1+m2) m2*l1*l2*cos(theta2-theta1);m2*l1*l2*cos(theta2-theta1) m2*l2^2];
res = [tau1-(m1+m2)*l1*g*sin(theta1)-m2*l1*l2*omega2^2*sin(theta1-theta2);tau2-m2*g*l2*sin(theta2)];
alphas = inv(Mat1)*res;

out = subs(alphas,[l1 l2 m1 m2 g],[1 1 1 1 9.8]);
xx_s = [theta1; theta2; omega1; omega2];
xxDot_s = [omega1; omega2; alpha1; alpha2];
uu_s = [tau1;tau2];
xxDotExpr_nn = [xxDot_s(1:2,1); out]; % + [zeros(2);eye(2)]*uu_s;
%% simulation
xxFunc_s = matlabFunction(xxDotExpr_nn, 'vars', {t, xx_s, uu_s});
[tt2,yy2] = ode45(@(t,x) xxFunc_s(t, x, [0;0]), 0:1e-2:10, [15;5; zeros(2,1)]/180*pi);

%% Euler-Lagrange Equation

syms t g real;
qq = sym("q%d__dt_0_", [2,1]);
uu = sym("tau%d", [2,1]);
mm = sym("m%d", [2,1]);
ll = sym("l%d", [2,1]);

II = diag(mm) * ll.^2;
angPos = qq + [0; qq(1)];
angVel = pDiff(angPos, t);
angAcc = pDiff(angVel, t);

xxCom = -[sin(angPos(1))*ll(1); sin(angPos(2))*ll(2)];
xxCom(2) = sum(xxCom);
xxcom_vel = pDiff(xxCom,t);

yyCom = -[cos(angPos(1))*ll(1); cos(angPos(2))*ll(2)];
yyCom(2) = sum(yyCom);
yycom_vel = pDiff(yyCom,t);

kEnergy = simplify(expand(0.5 * yycom_vel.' * diag(mm) * yycom_vel)+...
    expand(0.5 * xxcom_vel.' * diag(mm) * xxcom_vel));

kPotential = yyCom.' * mm * g;

lagran = kEnergy - kPotential;

eqRight = jacobian(lagran, qq).'; 
eqLeft = pDiff(jacobian(lagran, pDiff(qq,t)),t).';

xx = [qq; pDiff(qq,t)];
xxDot = pDiff(xx,t);

B_right = [-g*(mm(2)*(ll(2)*sin(xx(1)+xx(2))+ll(1)*sin(xx(1)))+mm(1)*ll(1)*sin(xx(1)))+xx(4)*ll(1)*ll(2)*mm(2)*sin(xx(2))*(2*xx(3)+xx(4));...
    -mm(2)*ll(2)*(ll(1)*sin(xx(2))*xx(3)*xx(3)+g*sin(xx(1)+xx(2)))];
A_left = [mm(1)*ll(1)^2+mm(2)*(ll(1)^2+ll(2)^2)+2*ll(1)*ll(2)*mm(2)*cos(xx(2)),...
    ll(2)*(ll(2)*mm(2)+ll(1)*mm(1)*cos(xx(2)));...
    mm(2)*ll(2)*(ll(2)+ll(1)*cos(xx(2))),mm(2)*ll(2)^2];
alphas2 = inv(A_left)*B_right;

xxDotExpr2 = [xxDot(1:2,1); alphas2];
xxFunc = matlabFunction(xxDotExpr2, 'vars', {t, xx, uu, [g; mm; ll]});
[tt,y]=ode45(@(t,x) xxFunc(t, x, zeros(2,1), [9.8; 1;1; 1;1]), 0:1e-2:10, [15;-10; zeros(2,1)]/180*pi);

%% Plots
figure(1);
plot(tt,y(:,1),tt2,yy2(:,1));
xlabel("time(s)")
ylabel("Pendulum 1")
legend('lagrangian','newtonian')
figure(2);
plot(tt,y(:,1)+y(:,2),tt2,yy2(:,2));
xlabel("time(s)")
ylabel("Pendulum 2")
legend('lagrangian','newtonian')
