function mimodmc = dmc_MIMO(G_sys,n,h,ySP,m,p,Q,R,duMin,duMax,uMin,uMax)
%% Error Checks
if size(ySP,1)~=size(G_sys,1)
    error('Error. Set points not equal to outputs')
end
if size(uMin,1)~=size(G_sys,2) || size(uMax,1)~=size(G_sys,2)||...
    size(duMin,1)~=size(G_sys,2) || size(duMax,1)~=size(G_sys,2)
    error('Error. Input constraints not same size as input')
end

%% Step Response Model

ny = size(G_sys,1); nu = size(G_sys,2);

%step-response data
[y,t]=step(G_sys,0:h:n*h);
S = zeros(n*ny,nu);
for i = 1:(length(t)-1)
    for j = 0:1:ny-1
        for k = 1:1:nu
            S(ny*i-j,k) = y(i+1,ny-j,k);
        end
    end         
end

%% Initialization
uPrev=zeros(nu,1);
Yk0=zeros(n*ny,1);
maxTime=100;

Y_SAVE=zeros((maxTime+1)*ny,1);
T_SAVE=zeros(maxTime+1,1);
U_SAVE=zeros(maxTime*nu,1);

%clear G y t

%% Pre-compute Matrices
%Su_col=S(1:p,:,:);
%Im_col=ones(m,1);
bigSu=zeros(p*ny,m*nu);
bigIm=zeros(m*nu,m*nu);
for i = 1:1:p
    for j = 1:1:m
        for k = 0:1:ny-1
            for l = 0:1:nu-1
                if j<=i
                    bigSu(ny*i-k,nu*j-l) = S(ny*(i-j+1)-k,nu-l);
                end
            end
        end
    end
end


for i =1:1:m
    for j = 1:1:m
        for l =0:1:nu-1
            if j<=i
                bigIm(nu*i-l,nu*j-l) = 1;
            end
        end
    end
end

% For Objective
bigR=repmat(ySP,p,1);
GammaY=zeros(ny*p,ny*p);
for i = 1:1:ny
    for j = 0:1:p-1
        GammaY(ny*j+i,ny*j+i)=Q(i,i);
    end
end
GammaU=zeros(nu*m,nu*m);
for i = 1:1:nu
    for j = 0:1:m-1
        GammaY(nu*j+i,nu*j+i)=R(i,i);
    end
end

% For constraints
Im=eye(m*nu);

% Pre-compute Hessian
Hess=bigSu'*GammaY*bigSu + GammaU;
% Pre-compute LHS of constraints
C_LHS=[Im;-Im;bigIm;-bigIm];

%% Implementation
bu = [];
for k=1:maxTime+1
    time=(k-1)*h;   % Current time 
        
    % Calculate error
    if (k==1)
        YHAT=Yk0;
        %err=zeros(ny,1);
    else
        %err=yk-YHAT(1:ny);    
    end
    
    
    % Calculate gradient
    predErr=YHAT(ny+1:ny*(p+1))-bigR;
    grad=bigSu'*GammaY*predErr;
    % Calculate RHS of constraint
    cRHS=[repmat(duMax,m,1); repmat(-duMin,m,1)];
    cRHS=[cRHS; repmat( (uMax-uPrev),m,1)];
    cRHS=[cRHS; repmat(-(uMin-uPrev),m,1)];
            
    % Calculating current step
    % big_dU=-inv(Hess)*grad;
    big_dU=quadprog(Hess,grad,C_LHS,cRHS);
    du=big_dU(1:nu);
    uk=uPrev+du;
    
    % Store results
    for i=1:1:nu
        U_SAVE(nu*(k-1)+i)=uk(i);
    end
    uPrev=uk;

    % Implementing current step
    YHAT=[YHAT(ny+1:end);YHAT(end-ny+1:end)] + S*du ;
    yk=YHAT(1:ny);
    
    % Storing results
    for i=1:1:ny
        Y_SAVE(ny*(k-1)+i)=yk(i);
    end
    T_SAVE(k+1)=k*h;
    bu = [bu big_dU];
end

mimodmc = U_SAVE;

%% Plotting results
figure
for sp = 1:1:nu
    subplot(nu,1,sp);
    stairs(T_SAVE(1:maxTime),U_SAVE(sp:nu:nu*maxTime));
    ylabel(sprintf('Input_{%d}',sp)); xlabel('time, t');
    xlim([0,h*maxTime]);
end
figure
for sp = 1:1:ny
    subplot(ny,1,sp);
    plot(T_SAVE(1:maxTime),Y_SAVE(sp:ny:ny*maxTime));
    ylabel(sprintf('Output_{%d}',sp)); xlabel('time, t');
    xlim([0,h*maxTime]);
end


%clear S n Y0 k du


