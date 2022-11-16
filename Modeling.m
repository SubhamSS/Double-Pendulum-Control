PathSetupSub;


%% 
t=sym('t','real');
System=kSystem('DoublePendulum',t);
WORLD=System.Space.RootFrame;
System.Space.setPrecompile(true);
System.Model.setPrecompile(true);


%% Declare Constant Parameters
C = System.genParam({ %Constant Parameter
             'empty' 'empty' 0; 
             });
         
[q1, q2]...
=System.genCont({
            'q1' 'Angle 1';
            'q2' 'Angle 2';
            });
        
P_Sys=System.genDisc({
             'g' 'Gravity Acceleration';
             'b' 'Global Damping Factor';
             });
        
P_Mech=System.genDisc({
             'm1' 'Link 1 End Point Mass';
             'm2' 'Link 2 End Point Mass';
             'l1' 'Link 1 Length';
             'l2' 'Link 2 Length';
             });
U_Mech=System.genInput({
             'u1' 'Link 1 Control Torque';
             'u2' 'Link 2 Control Torque';
             });
         
%% Declare Kinematics
[Frame0, Frame1COM, Frame2, Frame2COM]...
=System.Space.genNode({
                     'Frame0';
                     'Frame1COM';
                     'Frame1';
                     'Frame2COM';
                     });
                 
World2Frame0=System.Space.genEdge({'World2Frame1' 'WORLD' 'Frame0'});
World2Frame0.setAxang([q1.dot(0); 0; -1; 0]).genProp;
                 
Frame0ToFrame1COM=System.Space.genEdge({'Frame0ToFrame1COM' 'Frame0' 'Frame1COM'});
Frame0ToFrame1COM.setTransDis([0;0;-P_Mech.l1]).genProp;

Frame1COMToFrame1=System.Space.genEdge({'Frame1COMToFrame1' 'Frame1COM' 'Frame1'});
Frame1COMToFrame1.setAxang([q2.dot(0); 0; -1; 0]).genProp;

Frame1ToFrame2COM=System.Space.genEdge({'Frame1ToFrame2COM' 'Frame1' 'Frame2COM'});
Frame1ToFrame2COM.setTransDis([0;0;-P_Mech.l2]).genProp;

System.Space.plotGraph(1);
System.Space.makeNumKinematics();


%% Declare Dynamics
qqdot=System.getContVector(1); % 

%% Declare Link Bodies
[Link1, Link2]...
=System.genBody({
                'Link1' 'Frame1COM';
                'Link2' 'Frame2COM';
                });
            
Link1.setProp(P_Mech.m1, sym(zeros(3))).genForce({'gLink1' Frame1COM WORLD}).setProp([0;0;-1 * P_Mech.m1 * P_Sys.g]);
Link2.setProp(P_Mech.m2, sym(zeros(3))).genForce({'gLink2' Frame2COM WORLD}).setProp([0;0;-1 * P_Mech.m2 * P_Sys.g]);

%% Declare Torque Inputs
Link1.genTorque({'uLink1' Frame1COM Frame1COM}).setProp([0;-U_Mech.u1;0]);
Link2.genTorque({'uLink1' Frame1COM Frame1COM}).setProp([0;U_Mech.u2;0]);
Link2.genTorque({'uLink2' Frame2COM Frame2COM}).setProp([0;-U_Mech.u2;0]);

%% Declare Damping
System.genDamper({'Link Dampers' 2 1}).setProp(P_Sys.b, qqdot);

%% Compile Dynamics
System.Model.Init('Full');
[Default]=System.Model.genNode({'Default'});
System.Model.plotGraph(2);
Default.genFlowEOM({});
System.Model.makeDynamics('num');
