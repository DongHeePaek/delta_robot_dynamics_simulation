% author: Dong-Hee Paek, Yeong-Dae Kim
% e-mail: donghee.paek@kaist.ac.kr, 

clc, clear
close all

%%%%%%%%%%%%%%%%%%% [ Configuration ] %%%%%%%%%%%%%%%%%%%
PRECISE_DT = 0; ANIMATION_FAST_DT = 1; ANIMATION_SLOW_DT = 0;
CONFIG_TF = 20; % [s] Termination time of the simulation
CONFIG_ROT_NUM = 0.2; % [rev/s]
CIRCLE = 0; BUTTERFLY = 1; % Trajectory
PD = 1; PID = 0; % Control Type

%%%%%%%%%%%%%%%%%%% [ Parameter Initialization ] %%%%%%%%%%%%%%%%%%%
% global variables for ode 45 function
global a b L1 L2 m1 m2 mp m_m L1c I1 gm_Ia IA gc fx fy fz;

% Kinematic Parameter of Delta Robot Prototype
a = 0.20/sqrt(3); % [m] Radius of the fixed base
b = 0.06/sqrt(3); % [m] Radius of the moving flatform
L1 = 0.12;        % [m] Link1 length
L2 = 0.30;        % [m] Link2 length

% Mass Properties of Delta Robot Prototype
m1 = 0.200;       % 0.426 [kg] Link1 mass
m2 = 0.050;       % 0.069 [kg] Link2 mass
mp = 0.050;       % 0.096 [kg] Moving platform mass
m_m = mp + 3*m2;  % [kg]
L1c = 0.06;       % 0.06184 [m] Link1 mass center
I1 = 0.00098;     % 0.00398 [kg*m^2] Link1 inertia

gm_Ia = 25^2*(0.24+0.03)*10^(-4); % [kg*m^2] Gear and Motor inertia (Include Reduction Ratio)
IA = gm_Ia + I1 + m2*L1^2;        % [kg*m^2]

% Force Properties of Delta Robot Prototype;
gc = 9.806; % [m/s^2]
% External force to the End Effector: F = [fx, fy, fz]
fx = 0; fy = 0; fz = gc*0.3; % 300g mass (for Chocolate Extruder)

% Time Properties
% use small t for looping
% use ot for solving differential equation
if PRECISE_DT
    dT = 0.001;
elseif ANIMATION_FAST_DT
    dT = 0.01;
elseif ANIMATION_SLOW_DT
    dT = 0.002;
end
Ti = 0; Tf = CONFIG_TF; T = Ti:dT:Tf;
cntT = size(T); N = cntT(1,2);

%%% [ Memory Allocation ] %%%
% i: leg number idx

% Current position, speed, acceleration of End Effector (Cartesian)
% Idx // 1: x,y,z // 2: time
Y = zeros(3,N);
Yd = zeros(3,N);
Ydd = zeros(3,N);

% Dst ~ of End Effector (Cartesian)
% Idx // 1: x,y,z // 2: time
Y_dst = zeros(3,N);
Yd_dst = zeros(3,N);
Ydd_dst = zeros(3,N);

% Base Point (Fixed location)
% Idx // 1: xyz // 2: leg number
Ai = [a*cos(pi/3), -a,  a*cos(pi/3);
      a*sin(pi/3),  0, -a*sin(pi/3);
                0,  0,            0];

% Current M Point
% Idx // 1: xyz // 2: leg number // 3: time
Mi = zeros(3,3,N);

% Current Platform Point
% Idx // 1: xyz // 2: leg number // 3: time
Bi = zeros(3,3,N);

% Current angle, angular velocity, and accleration
% Idx // 1: leg number // 2: time
th1i = zeros(3,N); th2i = zeros(3,N); th3i = zeros(3,N);
th1id = zeros(3,N); th2id = zeros(3,N); th3id = zeros(3,N);
th1idd = zeros(3,N);

% Dst angle ~; derived from Inverse Kinematics
% Idx // 1: leg number // 2: time
th1i_dst = zeros(3,N); th2i_dst = zeros(3,N); th3i_dst = zeros(3,N);
th1id_dst = zeros(3,N); th2id_dst = zeros(3,N); th3id_dst = zeros(3,N);
th1idd_dst = zeros(3,N);

% Current Jacobian Matrix / Jacobian Dot Matrix
% Idx // 1,2: 3x3 matrix // 3: time
Jx = zeros(3,3,N); Jq = zeros(3,3,N); J = zeros(3,3,N);
Jxd = zeros(3,3,N); Jqd = zeros(3,3,N); Jd = zeros(3,3,N);

% Dst Jacobian Matrix / Jacobian Dot Matrix
% Idx // 1,2: 3x3 matrix // 3: time
Jx_dst = zeros(3,3,N); Jq_dst = zeros(3,3,N); J_dst = zeros(3,3,N);
Jxd_dst = zeros(3,3,N); Jqd_dst = zeros(3,3,N); Jd_dst = zeros(3,3,N);

% Current Dynamic Model Matrix
% Idx // 1,2: 3x3 matrix // 3: time
M = zeros(3,3,N); C = zeros(3,3,N);
% Idx // 1,2: 3x1 matrix // 3: time
G = zeros(3,1,N);

% Dst Dynamic Model Matrix
% Idx // 1,2: 3x3 matrix // 3: time
M_dst = zeros(3,3,N); C_dst = zeros(3,3,N);
% Idx // 1,2: 3x1 matrix // 3: time
G_dst = zeros(3,1,N);

% input motor Torque
% Idx // 1: leg number // 2: time
Trqi = zeros(3,N);

%%%%%%%%%%%%%%%%%%% [ Error Dynamics Gain ] %%%%%%%%%%%%%%%%%%%
%kp = 10; kd = 0;
kp = 150; kd = 100;
% 3x3 matrix
I33 = diag([1,1,1]);
Kp = kp*I33; Kd = kd*I33;

%%%%%%%%%%%%%%%%%%% [ Rotation Matrix ] %%%%%%%%%%%%%%%%%%%
% Idx: leg number
phi_i = [pi/3, pi, -pi/3];
% Idx // 1,2: 3x3 matrix // 3: leg number
Rbi = zeros(3,3,3); Rbitr = zeros(3,3,3);
for i = 1:3
    Rbi(:,:,i) = [ cos(phi_i(i)), -sin(phi_i(i)), 0;
                   sin(phi_i(i)),  cos(phi_i(i)), 0;
                                 0,               0, 1 ];
    Rbitr(:,:,i) = [ cos(phi_i(i)), sin(phi_i(i)), 0;
                    -sin(phi_i(i)), cos(phi_i(i)), 0;
                                   0,               0, 1 ];
end

%%%%%%%%%%%%%%%%%%% [ Initial Point & Theta ] %%%%%%%%%%%%%%%%%%%
Y_init = [ 0.01; 0.02; 0.34 ];
% i: leg number
for i = 1:3
    b1i_dst = Y_init(1)*cos(phi_i(i))+Y_init(2)*sin(phi_i(i))+b-a;
    b2i_dst = -Y_init(1)*sin(phi_i(i))+Y_init(2)*cos(phi_i(i));
    b3i_dst = Y_init(3);
    % acos returns 0 ~ pi
    th3i(i,1) = acos(b2i_dst/L2);
    K_dst = (b1i_dst^2+b2i_dst^2+b3i_dst^2-L1^2-L2^2) / (2*L1*L2*sin(th3i(i,1)));
    th2i(i,1) = acos(K_dst);
    g1i_dst = L1 + L2*cos(th2i(i,1))*sin(th3i(i,1));
    g2i_dst = L2*sin(th2i(i,1))*sin(th3i(i,1));
    th1i(i,1) = atan2( -g2i_dst*b1i_dst + g1i_dst*b3i_dst, g1i_dst*b1i_dst + g2i_dst*b3i_dst );
    
    iAiMi = L1*[cos(th1i(i,1)); 0; sin(th1i(i,1))];
    iMiBi = L2*[sin(th3i(i,1))*cos(th1i(i,1)+th2i(i,1));
                cos(th3i(i,1));
                sin(th3i(i,1))*sin(th1i(i,1)+th2i(i,1))];
    iAiBi = iAiMi + iMiBi;
    % set vectors for base coordinate
    bAiMi = Rbi(:,:,i)*iAiMi;
    bAiBi = Rbi(:,:,i)*iAiBi;
    
    % xyz matching, j=1:x, 2:y, 3:z coordinate
    for j = 1:3
        Mi(j,i,1) = Ai(j,i) + bAiMi(j);
        Bi(j,i,1) = Ai(j,i) + bAiBi(j);
    end
end

%%%%%%%%%%%%%%%%%%% [ Dst Caculation = Inverse Dynamics ] %%%%%%%%%%%%%%%%%%%
n = 1;  % time Idx initialization
Rot_num = CONFIG_ROT_NUM;
Omega = 2*pi*Rot_num;
for t=T
    %%%%%%%%%%%%%%%%%%% [ Trajectory ] %%%%%%%%%%%%%%%%%%%
    if CIRCLE
        Y_dst(:,n) = [0.03*sin(Omega*t); -0.03*cos(Omega*t); 0.31];
        Yd_dst(:,n) = [0.03*Omega*cos(Omega*t); 0.03*Omega*sin(Omega*t); 0];
        Ydd_dst(:,n) = [-0.03*Omega*Omega*sin(Omega*t); 0.03*Omega*Omega*cos(Omega*t); 0];
    elseif BUTTERFLY
        Y_dst(:,n) = [0.03*sin(Omega*t) ; -0.03*sin(2*Omega*t); 0.31];
        Yd_dst(:,n) = [0.03*Omega*cos(Omega*t); -0.03*2*Omega*cos(2*Omega*t); - 0.0012];
        Ydd_dst(:,n) = [-0.03*Omega*Omega*sin(Omega*t); 0.03*2*2*Omega*Omega*sin(2*Omega*t); 0];
    end
    
    %%%%%%%%%%%%%%%%%%% [ Inverse Kinematics ] %%%%%%%%%%%%%%%%%%%
    % i: leg number
    for i = 1:3
        b1i_dst = Y_dst(1,n)*cos(phi_i(i))+Y_dst(2,n)*sin(phi_i(i))+b-a;
        b2i_dst = -Y_dst(1,n)*sin(phi_i(i))+Y_dst(2,n)*cos(phi_i(i));
        b3i_dst= Y_dst(3,n);
        % acos returns 0 ~ pi
        th3i_dst(i,n) = acos(b2i_dst/L2);
        K_dst = (b1i_dst^2+b2i_dst^2+b3i_dst^2-L1^2-L2^2) / (2*L1*L2*sin(th3i_dst(i,n)));
        th2i_dst(i,n) = acos(K_dst);
        g1i_dst = L1 + L2*cos(th2i_dst(i,n))*sin(th3i_dst(i,n));
        g2i_dst = L2*sin(th2i_dst(i,n))*sin(th3i_dst(i,n));
        th1i_dst(i,n) = atan2( -g2i_dst*b1i_dst + g1i_dst*b3i_dst, g1i_dst*b1i_dst + g2i_dst*b3i_dst );
        
        %%%%%%%%%%%%%%%%%%% [ Jacobian Components ] %%%%%%%%%%%%%%%%%%%
        % TODO: Inserting current angle or dst angle for Jacobian?
        Jx_dst(i,1,n) = cos(th1i_dst(i,n)+th2i_dst(i,n))*sin(th3i_dst(i,n))*cos(phi_i(i))-cos(th3i_dst(i,n))*sin(phi_i(i));
        Jx_dst(i,2,n) = cos(th1i_dst(i,n)+th2i_dst(i,n))*sin(th3i_dst(i,n))*sin(phi_i(i))+cos(th3i_dst(i,n))*cos(phi_i(i));
        Jx_dst(i,3,n) = sin(th1i_dst(i,n)+th2i_dst(i,n))*sin(th3i_dst(i,n));
        Jq_dst(i,i,n) = L1*sin(th2i_dst(i,n))*sin(th3i_dst(i,n));
        % Inv operation after calculation of all legs
    end
    
    %%%%%%%%%%%%%%%%%%% [ Calculation of Jacobian ] %%%%%%%%%%%%%%%%%%%
    J_dst(:,:,n) = inv(Jq_dst(:,:,n))*Jx_dst(:,:,n); % TODO: try with other operation
    th1id_dst(:,n) = J_dst(:,:,n)*Yd_dst(:,n);
    
    %%%%%%%%%%%%%%%%%%% [ Calculation of Angular Velocities ] %%%%%%%%%%%%%%%%%%%
    % i: leg number
    for i = 1:3
        ipd_dst = Rbitr(:,:,i)*Yd_dst(:,n);
        th3id_dst(i,n) = -ipd_dst(2)/(L2*sin(th3i_dst(i,n)));
        % Temp numerator
        th2id_num_dst = -ipd_dst(1)*cos(th1i_dst(i,n))-ipd_dst(3)*sin(th1i_dst(i,n))+L2*cos(th2i_dst(i,n))*cos(th3i_dst(i,n))*th3id_dst(i,n);
        th2id_den_dst = L2*sin(th2i_dst(i,n))*sin(th3i_dst(i,n));
        th2id_dst(i,n) = th2id_num_dst/th2id_den_dst - th1id_dst(i,n);
        
        %%%%%%%%%%%%%%%%%%% [ Jacobian Dot Components ] %%%%%%%%%%%%%%%%%%%
        th1p2_dst = th1i_dst(i,n)+th2i_dst(i,n);
        thd1p2_dst = th1id_dst(i,n)+th2id_dst(i,n);
        temp_former_dst = sin(th1p2_dst)*thd1p2_dst*sin(th3i_dst(i,n));
        temp_mid_dst = cos(th1p2_dst)*th3id_dst(i,n)*cos(th3i_dst(i,n));
        temp_latter_dst = th3id_dst(i,n)*sin(th3i_dst(i,n));
        Jxd_dst(i,1,n) = -temp_former_dst*cos(phi_i(i)) + temp_mid_dst*cos(phi_i(i)) + temp_latter_dst*sin(phi_i(i));
        Jxd_dst(i,2,n) = -temp_former_dst*sin(phi_i(i)) + temp_mid_dst*sin(phi_i(i)) - temp_latter_dst*cos(phi_i(i));
        Jxd_dst(i,3,n) = cos(th1p2_dst)*thd1p2_dst*sin(th3i_dst(i,n))+sin(th1p2_dst)*th3id_dst(i,n)*cos(th3i_dst(i,n));
        Jqd_dst(i,i,n) = L1*(th2id_dst(i,n)*cos(th2i_dst(i,n))*sin(th3i_dst(i,n))+th3id_dst(i,n)*sin(th2i_dst(i,n))*cos(th3i_dst(i,n)));
    end
    
    %%%%%%%%%%%%%%%%%%% [ Caculation of Angular Acceleration ] %%%%%%%%%%%%%%%%%%%
    invJq_dst = inv(Jq_dst(:,:,n));
    th1idd_dst(:,n) = invJq_dst*(Jxd_dst(:,:,n)*Yd_dst(:,n)+Jx_dst(:,:,n)*Ydd_dst(:,n)-Jqd_dst(:,:,n)*th1id_dst(:,n));
    
    % Temp value for calculating maximum torque
    %%%%%%%%%%%%%%%%%%% [ Caculation of Dynamic Model ] %%%%%%%%%%%%%%%%%%%
    invJtr_dst = inv(J_dst(:,:,n)');
    invJ_dst = inv(J_dst(:,:,n));
    invJx_dst = inv(Jx_dst(:,:,n));
    dt_invJ_dst = -invJx_dst*Jxd_dst(:,:,n)*invJx_dst*Jq_dst(:,:,n)+invJx_dst*Jqd_dst(:,:,n);
    M_dst(:,:,n) = IA*I33+m_m*invJtr_dst*invJ_dst;
    C_dst(:,:,n) = m_m*invJtr_dst*dt_invJ_dst;
    % TODO: f3 + gc dimension?, incorporating gc term or not
    v1_dst = [cos(th1i_dst(1,n)); cos(th1i_dst(2,n)); cos(th1i_dst(3,n))];
    v2_dst = [fx; fy; fz + gc];
    G_dst(:,1,n) = -(m1*L1c+m2*L1)*gc*v1_dst - m_m*invJtr_dst*v2_dst;
    
    n = n + 1; % Iter
end

%%%%%%%%%%%%%%%%%%% [ Foward Kinematics Parameter ] %%%%%%%%%%%%%%%%%%%
% Idx // 1: xyz, 2: leg number
fkBias = [(a-b)*cos(pi/3), -a+b,  (a-b)*cos(pi/3);
          (a-b)*sin(pi/3),    0, -(a-b)*sin(pi/3);
                        0,    0,               0];

r1 = L2; r2 = L2; r3 = L2; thr = 0.00001;

%%%%%%%%%%%%%%%%%%% [ Simulation ] %%%%%%%%%%%%%%%%%%%
n = 1;
for t=T
    %%%%%%%%%%%%%%%%%%% [ Forward Kinematics, th1i ==> P ] %%%%%%%%%%%%%%%%%%%
    % Center of 3 Sphere
    fkC1 = L1*Rbi(:,:,1)*[cos(th1i(1,n)); 0; sin(th1i(1,n))] + fkBias(:,1);
    fkC2 = L1*Rbi(:,:,2)*[cos(th1i(2,n)); 0; sin(th1i(2,n))] + fkBias(:,2);
    fkC3 = L1*Rbi(:,:,3)*[cos(th1i(3,n)); 0; sin(th1i(3,n))] + fkBias(:,3);
    
    % 3 Sphere Algorithm
    x1 = fkC1(1); y1 = fkC1(2); z1 = fkC1(3);
    x2 = fkC2(1); y2 = fkC2(2); z2 = fkC2(3);
    x3 = fkC3(1); y3 = fkC3(2); z3 = fkC3(3);
    a11=2*(x3-x1); a12=2*(z3-z1); a13=2*(y3-y1);
    a21=2*(x3-x2); a22=2*(z3-z2); a23=2*(y3-y2);
    b1=r1^2-r3^2-x1^2-z1^2-y1^2+x3^2+z3^2+y3^2;
    b2=r2^2-r3^2-x2^2-z2^2-y2^2+x3^2+z3^2+y3^2;
    a1=a11/a13-a21/a23; a2=a12/a13-a22/a23; a3=b2/a23-b1/a13;
    a4=-a2/a1; a5=-a3/a1; a6=(-a21*a4-a22)/a23; a7=(b2-a21*a5)/a23;
    
    % Quadratic Coefficients
    qa = a4^2+1+a6^2; qb = 2*a4*(a5-x1)-2*z1+2*a6*(a7-y1);
    qc = a5*(a5-2*x1)+a7*(a7-2*y1)+x1^2+z1^2+y1^2-r1^2;
    quad = qb^2-4*qa*qc; % inside sqrt
    
    % Solution xyz = Complex Number
    sz = (-qb+sqrt(quad))/(2*qa); sx = a4*sz+a5; sy = a6*sz+a7;
    
    % Singularities
    %if(a13 < thr || a23 < thr || a1 < thr || quad < 0)
        %disp("singular point")
    %else
        Y(1,n) = real(sx);
        Y(2,n) = real(sy);
        Y(3,n) = real(sz);
    %end
    
    %%%%%%%%%%%%%%%%%%% [ th1i, P ==> th2i, th3i ] %%%%%%%%%%%%%%%%%%%
    % Idx: leg number
    for i = 1:3
        b1i_fk = Y(1,n)*cos(phi_i(i))+Y(2,n)*sin(phi_i(i))+b-a;
        b2i_fk = -Y(1,n)*sin(phi_i(i))+Y(2,n)*cos(phi_i(i));
        b3i_fk = Y(3,n);
        % acos returns 0 ~ pi
        th3i(i,n) = acos(b2i_fk/L2);
        K_fk = (b1i_fk^2+b2i_fk^2+b3i_fk^2-L1^2-L2^2) / (2*L1*L2*sin(th3i(i,n)));
        th2i(i,n) = acos(K_fk);
        
        %%%%%%%%%%%%%%%%%%% [ Jacobian Components ] %%%%%%%%%%%%%%%%%%%
        Jx(i,1,n) = cos(th1i(i,n)+th2i(i,n))*sin(th3i(i,n))*cos(phi_i(i))-cos(th3i(i,n))*sin(phi_i(i));
        Jx(i,2,n) = cos(th1i(i,n)+th2i(i,n))*sin(th3i(i,n))*sin(phi_i(i))+cos(th3i(i,n))*cos(phi_i(i));
        Jx(i,3,n) = sin(th1i(i,n)+th2i(i,n))*sin(th3i(i,n));
        Jq(i,i,n) = L1*sin(th2i(i,n))*sin(th3i(i,n));
        % Inv operation after calculation of all legs
    end
    
    %%%%%%%%%%%%%%%%%%% [ Calculation of Jacobian ] %%%%%%%%%%%%%%%%%%%
    J(:,:,n) = inv(Jq(:,:,n))*Jx(:,:,n); % TODO: another inv operation
    invJ = inv(J(:,:,n));
    Yd(:,n) = invJ*th1id(:,n);
    
    %%%%%%%%%%%%%%%%%%% [ th1id ==> th2id, th3id ] %%%%%%%%%%%%%%%%%%%
    % i: leg number
    for i = 1:3
        ipd = Rbitr(:,:,i)*Yd(:,n);
        th3id(i,n) = -ipd(2)/(L2*sin(th3i(i,n)));
        % Temp numerator
        th2id_num = -ipd(1)*cos(th1i(i,n))-ipd(3)*sin(th1i(i,n))+L2*cos(th2i(i,n))*cos(th3i(i,n))*th3id(i,n);
        th2id_den = L2*sin(th2i(i,n))*sin(th3i_dst(i,n));
        th2id(i,n) = th2id_num/th2id_den - th1id(i,n);
        
        %%%%%%%%%%%%%%%%%%% [ Jacobian Dot Components ] %%%%%%%%%%%%%%%%%%%
        th1p2 = th1i(i,n)+th2i(i,n);
        thd1p2 = th1id(i,n)+th2id(i,n);
        temp_former = sin(th1p2)*thd1p2*sin(th3i(i,n));
        temp_mid = cos(th1p2)*th3id(i,n)*cos(th3i(i,n));
        temp_latter = th3id(i,n)*sin(th3i(i,n));
        Jxd(i,1,n) = -temp_former*cos(phi_i(i)) + temp_mid*cos(phi_i(i)) + temp_latter*sin(phi_i(i));
        Jxd(i,2,n) = -temp_former*sin(phi_i(i)) + temp_mid*sin(phi_i(i)) - temp_latter*cos(phi_i(i));
        Jxd(i,3,n) = cos(th1p2)*thd1p2*sin(th3i(i,n))+sin(th1p2)*th3id(i,n)*cos(th3i(i,n));
        Jqd(i,i,n) = L1*(th2id(i,n)*cos(th2i(i,n))*sin(th3i(i,n))+th3id(i,n)*sin(th2i(i,n))*cos(th3i(i,n)));
    end
    invJq = inv(Jq(:,:,n));
    invJx = inv(Jx(:,:,n));
    
    %%%%%%%%%%%%%%%%%%% [ Caculation of Dynamic Model ] %%%%%%%%%%%%%%%%%%%
    invJtr = inv(J(:,:,n)');
    invJ = inv(J(:,:,n));
    invJx = inv(Jx(:,:,n));
    dt_invJ = -invJx*Jxd(:,:,n)*invJx*Jq_dst(:,:,n)+invJx*Jqd(:,:,n);
    M(:,:,n) = IA*I33+m_m*invJtr*invJ;
    C(:,:,n) = m_m*invJtr*dt_invJ;
    % f3 + gc dimension
    v1 = [cos(th1i(1,n)); cos(th1i(2,n)); cos(th1i(3,n))];
    v2 = [fx; fy; fz + gc];
    G(:,1,n) = -(m1*L1c+m2*L1)*gc*v1 - m_m*invJtr*v2;
    
    Th = [th1i(1,n); th1i(2,n); th1i(3,n)];
    Thd = [th1id(1,n); th1id(2,n); th1id(3,n)];
    
    E = th1i_dst(:,n) - th1i(:,n);
    Ed = th1id_dst(:,n) - th1id(:,n);
    
    % Control input
    U = th1idd_dst(:,n) + Kd*Ed + Kp*E;
    
    % Input torque
    Trqi(:,n) = M(:,:,n)*U+C(:,:,n)*Thd+G(:,1,n);
    
    th1idd(:,n) = inv(M(:,:,n))*(Trqi(:,n)-C(:,:,n)*Thd-G(:,1,n));
    
    State = [ th1id(1,n); th1idd(1,n); th1id(2,n); th1idd(2,n); th1id(3,n); th1idd(3,n) ];
    
    rk4_1 = @(t, th) State(1,1);
    rk4_2 = @(t, thd) State(2,1);
    rk4_3 = @(t, th) State(3,1);
    rk4_4 = @(t, thd) State(4,1);
    rk4_5 = @(t, th) State(5,1);
    rk4_6 = @(t, thd) State(6,1);
    
    k1 = rk4_1(t, th1i(1,n));
    k2 = rk4_1(t+dT*0.5, th1i(1,n)+dT*k1*0.5);
    k3 = rk4_1(t+dT*0.5, th1i(1,n)+dT*k2*0.5);
    k4 = rk4_1(t+dT, th1i(1,n)+dT*k3);
    th1i(1,n+1) = th1i(1,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    k1 = rk4_3(t, th1i(2,n));
    k2 = rk4_3(t+dT*0.5, th1i(2,n)+dT*k1*0.5);
    k3 = rk4_3(t+dT*0.5, th1i(2,n)+dT*k2*0.5);
    k4 = rk4_3(t+dT, th1i(2,n)+dT*k3);
    th1i(2,n+1) = th1i(2,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    k1 = rk4_5(t, th1i(3,n));
    k2 = rk4_5(t+dT*0.5, th1i(3,n)+dT*k1*0.5);
    k3 = rk4_5(t+dT*0.5, th1i(3,n)+dT*k2*0.5);
    k4 = rk4_5(t+dT, th1i(3,n)+dT*k3);
    th1i(3,n+1) = th1i(3,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    k1 = rk4_2(t, th1id(1,n));
    k2 = rk4_2(t+dT*0.5, th1id(1,n)+dT*k1*0.5);
    k3 = rk4_2(t+dT*0.5, th1id(1,n)+dT*k2*0.5);
    k4 = rk4_2(t+dT, th1id(1,n)+dT*k3);
    th1id(1,n+1) = th1id(1,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    k1 = rk4_4(t, th1id(2,n));
    k2 = rk4_4(t+dT*0.5, th1id(2,n)+dT*k1*0.5);
    k3 = rk4_4(t+dT*0.5, th1id(2,n)+dT*k2*0.5);
    k4 = rk4_4(t+dT, th1id(2,n)+dT*k3);
    th1id(2,n+1) = th1id(2,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    k1 = rk4_6(t, th1id(3,n));
    k2 = rk4_6(t+dT*0.5, th1id(3,n)+dT*k1*0.5);
    k3 = rk4_6(t+dT*0.5, th1id(3,n)+dT*k2*0.5);
    k4 = rk4_6(t+dT, th1id(3,n)+dT*k3);
    th1id(3,n+1) = th1id(3,n)+dT*(k1+2*k2+2*k3+k4)/6;
    
    %%%%%%%%%%%%%%%%%%% [ Caculation of Points ] %%%%%%%%%%%%%%%%%%%
    for i = 1:3
        % Points (dst-> current)
        % set vectors for i coordinate
        iAiMi = L1*[cos(th1i(i,n)); 0; sin(th1i(i,n))];
        iMiBi = L2*[sin(th3i(i,n))*cos(th1i(i,n)+th2i(i,n));
                    cos(th3i(i,n));
                    sin(th3i(i,n))*sin(th1i(i,n)+th2i(i,n))];
        iAiBi = iAiMi + iMiBi;
        % set vectors for base coordinate
        bAiMi = Rbi(:,:,i)*iAiMi;
        bAiBi = Rbi(:,:,i)*iAiBi;
        
        % xyz matching, j=1:x, 2:y, 3:z coordinate
        for j = 1:3
            Mi(j,i,n) = Ai(j,i) + bAiMi(j);
            Bi(j,i,n) = Ai(j,i) + bAiBi(j);
        end
    end
    
    n = n + 1; % Idx
end

%%%%%%%%%%%%%%%%%%% [ Calculation of Power ] %%%%%%%%%%%%%%%%%%%
Power = zeros(3,N);
for n = 1:N
    for i = 1:3
        Power(i,n) = th1id_dst(i,n)*Trqi(i,n);
    end
end

%%%%%%%%%%%%%%%%%%% [ Plot Input Torque, Velocity, Power ] %%%%%%%%%%%%%%%%%%%
P_dyn = figure(1);
set(P_dyn, 'Position', [10,50,390,930])
subplot(313)
plot(T,Power(1,:))
hold on;
grid on;
plot(T,Power(2,:))
plot(T,Power(3,:))
title('Power')
legend('Motor_1','Motor_2','Motor_3')
subplot(312)
plot(T,Trqi(1,:))
hold on;
grid on;
plot(T,Trqi(2,:))
plot(T,Trqi(3,:))
title('Input Torque')
legend('Motor_1','Motor_2','Motor_3')
subplot(311)
Pth1i = zeros(3,N); % one more component for th1id
Pth1id = zeros(3,N);
for pn = 1:N
    for i = 1:3
        Pth1i(i,pn) = th1i(i,pn);
        Pth1id(i,pn) = th1id(i,pn);
    end
end
plot(T,Pth1id(1,:))
hold on;
grid on;
plot(T,Pth1id(2,:))
plot(T,Pth1id(3,:))
title('Motor Velocity')
legend('Motor_1','Motor_2','Motor_3')

%%%%%%%%%%%%%%%%%%% [ Plot Track ] %%%%%%%%%%%%%%%%%%%
P_Trk = figure(2);
set(P_Trk, 'Position', [410,50,600,930])
subplot(3,3,1)
plot(T,th1i_dst(1,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1i(1,:),'-g')
title('Theta_1')
legend('DST', 'Theta_1')
subplot(3,3,2)
plot(T,th1i_dst(2,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1i(2,:),'-g')
title('Theta_2')
legend('DST', 'Theta_2')
subplot(3,3,3)
plot(T,th1i_dst(3,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1i(3,:),'-g')
title('Theta_3')
legend('DST', 'Theta_3')
subplot(3,3,4)
plot(T,th1id_dst(1,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1id(1,:),'-g')
title('ThetaDot_1')
legend('DST', 'ThetaDot_1')
subplot(3,3,5)
plot(T,th1id_dst(2,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1id(2,:),'-g')
title('ThetaDot_2')
legend('DST', 'ThetaDot_2')
subplot(3,3,6)
plot(T,th1id_dst(3,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,Pth1id(3,:),'-g')
title('ThetaDot_3')
legend('DST', 'ThetaDot_3')
subplot(3,3,7)
plot(T,th1idd_dst(1,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,th1idd(1,:),'-g')
title('Theta2Dot_1')
legend('DST', 'Theta2Dot_1')
subplot(3,3,8)
plot(T,th1idd_dst(2,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,th1idd(2,:),'-g')
title('Theta2Dot_2')
legend('DST', 'Theta2Dot_2')
subplot(3,3,9)
plot(T,th1idd_dst(3,:),'-r','LineWidth', 1.2)
hold on; grid on;
plot(T,th1idd(3,:),'-g')
title('Theta2Dot_3')
legend('DST', 'Theta2Dot_3')

%%%%%%%%%%%%%%%%%%% [ Plot 3D Model ] %%%%%%%%%%%%%%%%%%%
P_ani3D = figure(3);
set(P_ani3D, 'Position', [1010,50,890,930])
plot3(Y_dst(1,:),Y_dst(2,:),-Y_dst(3,:), '-b', 'Linewidth', 1)
xlim([-0.25,0.25])
ylim([-0.25,0.25])
zlim([-0.4,0.1])
hold on
grid on
PYdst = plot3(Y_dst(1,1), Y_dst(2,1), -Y_dst(3,1),'ro');
PY = plot3(Y(1,1), Y(2,1), -Y(3,1), 'k*');

% Base (Ai)
plot3([Ai(1,1), Ai(1,2)], [Ai(2,1), Ai(2,2)], [-Ai(3,1), -Ai(3,2)], '-r', 'Linewidth', 3)
plot3([Ai(1,2), Ai(1,3)], [Ai(2,2), Ai(2,3)], [-Ai(3,2), -Ai(3,3)], '-r', 'Linewidth', 3)
plot3([Ai(1,3), Ai(1,1)], [Ai(2,3), Ai(2,1)], [-Ai(3,3), -Ai(3,1)], '-r', 'Linewidth', 3)

% 1st, 2nd Link (Mi, Bi)
PM = zeros(3);
PB = zeros(3);
stT = text(0.1,-0.2,0.1,'');
for i = 1:3
    PM(i) = plot3([Ai(1,i), Mi(1,i,1)], [Ai(2,i), Mi(2,i,1)], [-Ai(3,i), -Mi(3,i,1)], '-k', 'Linewidth', 5);
    PB(i) = plot3([Mi(1,i,1), Bi(1,i,1)], [Mi(2,i,1), Bi(2,i,1)], [-Mi(3,i,1), -Bi(3,i,1)], '-k', 'Linewidth', 2);
end

% Platform; connect 3 points of Bi
P1 = plot3([Bi(1,1,1), Bi(1,2,1)], [Bi(2,1,1), Bi(2,2,1)], [-Bi(3,1,1), -Bi(3,2,1)], '-g', 'Linewidth', 3);
P2 = plot3([Bi(1,2,1), Bi(1,3,1)], [Bi(2,2,1), Bi(2,3,1)], [-Bi(3,2,1), -Bi(3,3,1)], '-g', 'Linewidth', 3);
P3 = plot3([Bi(1,3,1), Bi(1,1,1)], [Bi(2,3,1), Bi(2,1,1)], [-Bi(3,3,1), -Bi(3,1,1)], '-g', 'Linewidth', 3);

xlabel('x');
ylabel('y');
zlabel('z');

%figure(4)
%plot(T,th1i_dst(1,:))
%hold on
%grid on
%plot(T,th1i_dst(2,:))
%plot(T,th1i_dst(3,:))
%legend('Th1', 'Th2', 'Th3')

%%%%%%%%%%%%%%%%%%% [ Animation ] %%%%%%%%%%%%%%%%%%%
while(1)
n = 1;
for t=T
    drawnow;
    for i = 1:3
        set(PM(i), 'xdata', [Ai(1,i), Mi(1,i,n)]...
                  ,'ydata', [Ai(2,i), Mi(2,i,n)]...
                  ,'zdata', [-Ai(3,i), -Mi(3,i,n)]);
        set(PB(i), 'xdata', [Mi(1,i,n), Bi(1,i,n)]...
                  ,'ydata', [Mi(2,i,n), Bi(2,i,n)]...
                  ,'zdata', [-Mi(3,i,n), -Bi(3,i,n)]);
        set(P1, 'xdata', [Bi(1,1,n), Bi(1,2,n)]...
               ,'ydata', [Bi(2,1,n), Bi(2,2,n)]...
               ,'zdata', [-Bi(3,1,n), -Bi(3,2,n)]);
        set(P2, 'xdata', [Bi(1,2,n), Bi(1,3,n)]...
               ,'ydata', [Bi(2,2,n), Bi(2,3,n)]...
               ,'zdata', [-Bi(3,2,n), -Bi(3,3,n)]);
        set(P3, 'xdata', [Bi(1,3,n), Bi(1,1,n)]...
               ,'ydata', [Bi(2,3,n), Bi(2,1,n)]...
               ,'zdata', [-Bi(3,3,n), -Bi(3,1,n)]);
        %temp3(i) = (Bi(i,1,1) + Bi(i,2,1) + Bi(i,3,1))/3;
    end
    set(PYdst, 'xdata', Y_dst(1,n)...
              ,'ydata', Y_dst(2,n)...
              ,'zdata', -Y_dst(3,n));
    set(PY, 'xdata', Y(1,n)...
           ,'ydata', Y(2,n)...
           ,'zdata', -Y(3,n));
    set(stT,'string',sprintf('t = %.2f',t));
    
    pause(0.01)
    
    %xlabel(['x = ', num2str(Y(1,n))]);
    %ylabel(['y = ', num2str(Y(2,n))]);
    %zlabel(['z = ', num2str(Y(3,n))]);
    n = n+1;
end

end