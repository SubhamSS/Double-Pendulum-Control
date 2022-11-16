PathSetupSub; % Include folders to path directory and reset workspace

%%
Sim=sAxes('Double Pendulum Simulation',3,'Path_DoublePendulum.mat',@numTF_DoublePendulum_mex);
Sim.setAxesProp('WORLD',1*[-1 -1 -1; 1 1 1],[0 0]).setPresetTF('WORLD');
xlabel('x');

[Links, Link1Traj, Link2Traj] = Sim.genPlot({'Links'; 'Link1Traj'; 'Link2Traj'});
Links.setLineSpec('-','*','r',2).setPlotPoint({'WORLD';'Frame1COM';'Frame2COM'});
Link1Traj.setLineSpec(':','none','b',1).setTimedPlotPoint('Frame1COM');
Link2Traj.setLineSpec(':','none','k',1).setTimedPlotPoint('Frame2COM');

%% Load Model Info

load('ModelInfo_DoublePendulum.mat');

pp = [9.81; 0.05; 1; 1; 0.3; 0.3];
MM = ModelInfo.Dynamics.InertialMatrix;




%%
Link1Traj.clearTimedPlot();
Link2Traj.clearTimedPlot();

hh = 1e-2;
tSpan = 0:hh:120; %Simulation Time Span
xx = [1; 2; 0; 0];
uu = [0;0];
nh = 0;
cF = 0;
flow = 1;

kp = 5;
kd = 1;
ki = 0.5;

ref = [pi; 0; 0; 0;];
ee = xx - ref;
eeInt = zeros(4,1);

Sim.drawNow(flow,xx,pp,uu,nh); % requires 3D environment

ii = 0;
for tt = tSpan
    ii = ii+1;
    
    ref = [pi;0;0;0];
    ee = xx - ref;
    eeInt = eeInt + ee*hh;
    
    [~,~,~,~,~,~,~,~,GFF] = System_DoublePendulum_mex(tt,xx,pp,uu,nh);
    
    uu = - kp * ee(1:2) - kd * ee(3:4) - (sum(GFF,2));
    
    [xx,nh,cF]=rk4Hybrid(@Flow_DoublePendulum_mex,hh,flow,tt,xx,pp,uu,nh,cF);
    
    if rem(ii,floor(100/30)) == 1
        Sim.drawNow(flow,xx,pp,uu,nh); % requires 3D environment
        pause(floor(100/30)*hh*0.125);
    end
end

