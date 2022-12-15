% This code contains the EKF based Attitude, Heading, and Reference System (AHRS) for a fighter-jet class of aircraft.
% Sensors: GPS, IMU, Magnetometer, Pitot-tube 

clc
clearvars 
close all

%% Load all data

[acc_imu, gps_vel, ang_rates, psi_mag,ang_imu,ang_hyb,gps_pose,imu_vel, dt, t] = flight_data();
%count_dt = 1;
% for ii = 1:length(dt)
%     if dT(ii)~=0
%         dt_gps(count_dt) = dT(ii);
%         gps_vel(1,count_dt) = GPS_vel(1,:);
%         count_dt = count_dt+1;
%     end
% end
        

dt = ones(10000)*0.02;
t = 0:dt(1):1800;

acc_imu_x = acc_imu(1,:);      % Acc Data, z Up to down
acc_imu_y = acc_imu(2,:);
acc_imu_z = acc_imu(3,:);    

p_m = ang_rates(1,:);          % Gyro Data, Z down 
q_m = ang_rates(2,:);
r_m = ang_rates(3,:);

ang_imu_phi = ang_imu(1,:);      % IMU angle Data
ang_imu_theta = ang_imu(2,:);
ang_imu_psi = ang_imu(3,:);

ang_hyb_phi = ang_hyb(1,:);      % hyb angle Data
ang_hyb_theta = ang_hyb(2,:);
ang_hyb_psi = ang_hyb(3,:);


Gps_data_rate = 1e10/4; % every "" seconds
Imu_data_rate = 1/100; % every "" seconds
Mag_data_rate = 1/100; % every "" seconds
Pitot_data_rate = 1e10; % every "" seconds, high value if no data


gps_vel_N(1) = 0;
gps_vel_E(1) = 0;
gps_vel_D(1) = 0;


vel_gps_N = gps_vel(1,:);      % GPS Data, z Up to down
vel_gps_E = gps_vel(2,:); 
vel_gps_D = gps_vel(3,:); 

%
dv_N = diff(vel_gps_N);
dv_E = diff(vel_gps_E);
dv_D = diff(vel_gps_D);

acc_gps_N = dv_N./dt(1);   % Acc derived from GPS
acc_gps_E = dv_E./dt(1);
acc_gps_D = dv_D./dt(1);

g = 9.81;
g_vec = [0;0;g];

count_gps = 1; 
count_mag = 1; 
count_pitot = 1;
count_imu = 1;
count_pred = 1;
t_gps(count_gps) = 0;
t_imu(count_imu) = 0;
t_mag(count_mag) = 0;
t_pitot(count_pitot) = 0;

%% Parameter Initialization

b_p = 0.001;
b_q = 0.001;
b_r = 0.001;

bias = [b_p,b_q,b_r];

[X,P] = EKF_init(bias);           % Filter Initialization

XX(:,1) = X;
PP{1} = P;

phi   = X(1);
theta = X(2);
psi   = X(3);
b_p   = X(4);
b_q   = X(5);
b_r   = X(6); 

for i=1:length(t)
  
%% Prediction
% Continuous Time  x_dot = F_nlin*x + G*w
% Discrete Time  x(k+1) = (1 + F_nlin*dt)*x(k) + G*dt*w(k) = Phi(k)*x(k) + Tau(k)*w(k)

%count_pred = count_pred +1;

F_nlin = [(p_m(i) - b_p) + (q_m(i)-b_q)*sin(phi)*tan(theta) + (r_m(i)-b_r)*cos(phi)*tan(theta);
          (q_m(i)-b_q)*cos(phi) - (r_m(i)-b_r)*sin(phi);
          (q_m(i)-b_q)*sin(phi)/cos(theta) + (r_m(i)-b_r)*cos(phi)/cos(theta);
          0;0;0];
     
F_lin = [((q_m(i)-b_q)*cos(phi)*tan(theta) - (r_m(i)-b_r)*sin(phi)*tan(theta)),  (sec(theta)^2)*((q_m(i)-b_q)*sin(phi) + (r_m(i)-b_r)*cos(phi)),    0, -1,  -sin(phi)*tan(theta),  -cos(phi)*tan(theta);
         (-(q_m(i)-b_q)*sin(phi) - (r_m(i)-b_r)*cos(phi)),  0,        0,        0,     -cos(phi),     sin(phi);
         ((q_m(i)-b_q)*cos(phi)/cos(theta) - (r_m(i)-b_r)*sin(phi)/cos(theta)),   ((q_m(i)-b_q)*sin(phi)*tan(theta)/cos(theta) + (r_m(i)-b_r)*cos(phi)*tan(theta)/cos(theta)),      0,   0,  -sin(phi)/cos(theta),    -cos(phi)/cos(theta);
         0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 0];
     
G = [-1 -sin(phi)*tan(theta)  -cos(phi)*tan(theta) 0 0 0;
      0 -cos(phi)              sin(phi)            0 0 0;
      0 -sin(phi)/cos(theta)  -cos(phi)/cos(theta) 0 0 0;
      0  0                    0                    1 0 0;
      0  0                    0                    0 1 0;
      0  0                    0                    0 0 1 ];

err_bound_phi = deg2rad(1);
err_bound_theta = deg2rad(1);
err_bound_psi = deg2rad(0.5);
eps_phi_rate = err_bound_phi;
eps_theta_rate = err_bound_theta;
eps_psi_rate = err_bound_psi;
eps_b_p = 1e-2;
eps_b_q = 1e-2;
eps_b_r = 1e-2;

Q = diag([eps_phi_rate^2,eps_theta_rate^2,eps_psi_rate^2,eps_b_p^2,eps_b_q^2,eps_b_r^2]);
 
A_stm = eye(6) + F_lin*dt(i);

[X, P] = EKF_prediction(X,P,F_nlin,A_stm,G,Q,dt(i));

phi = X(1); theta = X(2); psi = X(3); b_p = X(4); b_q = X(5); b_r = X(6);

%% Update

del_t_gps = t(i) - t_gps(count_gps);
del_t_pitot = t(i) - t_pitot(count_pitot);
del_t_mag = t(i) - t_mag(count_mag);

if del_t_gps >= Gps_data_rate   % Check if GPS data_received, perform update
    
    count_gps = count_gps + 1;  
    t_gps(count_gps) = t(i);  
    
%     gps_pose_N(count_gps)  = gps_pose(1,i);
%     gps_pose_E(count_gps)  = gps_pose(2,i);
%     gps_pose_D(count_gps)  = gps_pose(3,i);
    
    gps_vel_N(count_gps) = vel_gps_N(i);%(gps_pose_N(count_gps) - gps_pose_N(count_gps-1))/del_t_gps;
    gps_acc_N(count_gps) = acc_gps_N(i);%(gps_vel_N(count_gps) - gps_vel_N(count_gps-1))/del_t_gps;
    gps_vel_E(count_gps) = vel_gps_E(i);%(gps_pose_E(count_gps) - gps_pose_E(count_gps-1))/del_t_gps;
    gps_acc_E(count_gps) = acc_gps_E(i);%(gps_vel_E(count_gps) - gps_vel_E(count_gps-1))/del_t_gps;
    gps_vel_D(count_gps) = vel_gps_D(i);%(gps_pose_D(count_gps) - gps_pose_D(count_gps-1))/del_t_gps;
    gps_acc_D(count_gps) = acc_gps_D(i);%(gps_vel_D(count_gps) - gps_vel_D(count_gps-1))/del_t_gps;

       
    % GPS-IMU Measurement model
  %  psi_m = mod(atan2(vel_gps_y(i),vel_gps_x(i)),2*pi);
     psi_m = mod(atan2(gps_vel_E(count_gps),gps_vel_N(count_gps)),2*pi);
     
%     r_x_m = acc_gps_x(i)*cos(psi_m) + acc_gps_y(i)*sin(psi_m);% gps acc ?
%     r_y_m = -acc_gps_x(i)*sin(psi_m) + acc_gps_y(i)*cos(psi_m);
%     r_z_m = acc_gps_z(i) - g;
      
    r_x_m = gps_acc_N(count_gps)*cos(psi_m) + gps_acc_E(count_gps)*sin(psi_m);% gps acc ?
    r_y_m = -gps_acc_N(count_gps)*sin(psi_m) + gps_acc_E(count_gps)*cos(psi_m);
    r_z_m = gps_acc_D(count_gps) - g;
    exprs_1_m = 1*sqrt(r_x_m^2 + r_z_m^2 - acc_imu_x(i)^2);
      
    theta_m = atan((-r_x_m*r_z_m - acc_imu_x(i)*exprs_1_m)/(acc_imu_x(i)^2 - r_z_m^2));
    
    r_theta_m = r_x_m*sin(theta_m) + r_z_m*cos(theta_m);
    exprs_2_m = 1*sqrt(r_y_m^2 + r_theta_m^2 - acc_imu_y(i)^2);
  
    phi_m = atan((r_theta_m*r_y_m + acc_imu_y(i)*exprs_2_m)/(acc_imu_y(i)^2 - r_theta_m^2));

    Z_gps(:,count_gps) = [phi_m;theta_m;psi_m];
    
    psi_comp = mod(psi,2*pi);
%     r_x = acc_gps_x(i)*cos(psi_comp) + acc_gps_y(i)*sin(psi_comp);
%     r_y = -acc_gps_x(i)*sin(psi_comp) + acc_gps_y(i)*cos(psi_comp);
%     r_z = acc_gps_z(i) - g;
    r_x = gps_acc_N(count_gps)*cos(psi_comp) + gps_acc_E(count_gps)*sin(psi_comp);
    r_y = -gps_acc_N(count_gps)*sin(psi_comp) + gps_acc_E(count_gps)*cos(psi_comp);
    r_z = gps_acc_D(count_gps) - g;    
    exprs_1 = 1*sqrt(r_x^2 + r_z^2 - acc_imu_x(i)^2);
    theta_comp = atan((-r_x*r_z - acc_imu_x(i)*exprs_1)/(acc_imu_x(i)^2 - r_z^2));
    
    r_theta = r_x*sin(theta_comp) + r_z*cos(theta_comp);
    
    exprs_2 = 1*sqrt(r_y^2 + r_theta^2 - acc_imu_y(i)^2);
    phi_comp = atan((r_theta*r_y + acc_imu_y(i)*exprs_2)/(acc_imu_y(i)^2 - r_theta^2));
        
    
    h_gps = [phi_comp;theta_comp;psi_comp];
    H_gps = [eye(3,3) zeros(3,3)];

    sigma_phi_m = deg2rad(2);
    sigma_theta_m = deg2rad(2);
    sigma_psi_m = deg2rad(1);

    R_gps = diag([sigma_phi_m^2,sigma_theta_m^2,sigma_psi_m^2]);
    sensor = 'GPS'; 
    
    [X, P] = EKF_update(X,P,Z_gps(:,count_gps),H_gps,h_gps,R_gps,sensor);
    
    phi = X(1); theta = X(2); psi = X(3); b_p = X(4); b_q = X(5); b_r = X(6); 
    
end

if del_t_mag >= Mag_data_rate   % Check if MAG data_received, perform update
    
    count_mag = count_mag + 1; 
    t_mag(count_mag) = t(i);  
    Z_mag(count_mag) = psi_mag(i);
    h_mag = psi;
    H_mag = [0 0 1 0 0 0];
    
    sigma_psi_mag = deg2rad(1.5);
    R_mag = sigma_psi_mag^2;
    
    sensor = 'MAG';
    [X, P] = EKF_update(X,P,Z_mag(count_mag),H_mag,h_mag,R_mag,sensor);
    
    phi = X(1); theta = X(2); psi = X(3); b_p = X(4); b_q = X(5); b_r = X(6); 

end

if del_t_pitot >= Pitot_data_rate   % Check if Pitot data_received, perform update
    
    disp('inside_pitot')
    pause(2)
    count_pitot = count_pitot + 1; 
    t_pitot(count_pitot) = t(i);   
    h_pitor = psi_mag(i);
    H_pitot = [0 0 1 0 0 0];
    
    sigma_velx_pitot = rand;
    sigma_vely_pitot = rand;
    R_pitot = diag([sigma_velx_pitot, sigma_vely_pitot]);
    Z_pitot(count_pitot) = [0,0];
    sensor = 'Pitot';
    [X, P] = EKF_update(X,P,Z_pitot(count_pitot),H_pitot,h_pitot,R_pitot,sensor);
    
    phi = X(1); theta = X(2); psi = X(3); b_p = X(4); b_q = X(5); b_r = X(6); 
end
 

XX(:,i+1) = X;
PP{i+1} = P;
end
 %gps_acc_x_f = medfilt1(gps_acc_E,3);
% acc_gps_x_f  = medfilt1(acc_gps_x,2);

figure(200)

grid on
subplot(3,1,1)
plot(t,rad2deg(ang_imu_phi(1:length(t))),'k',t,rad2deg(ang_hyb_phi(1:length(t))),'c.',t,rad2deg(XX(1,1:length(t))),'b');
xlim([0 1600]);
legend('IMU', 'HYB','COMPUTED')
%title('Velocity East')
xlabel('Time (s)')
ylabel('Roll: \phi (\circ)')

subplot(3,1,2)
plot(t,rad2deg(ang_imu_theta(1:length(t))),'k',t,rad2deg(ang_hyb_theta(1:length(t))),'c.',t,rad2deg(XX(2,1:length(t))),'b');
xlim([0 1600]);
legend('IMU', 'HYB','COMPUTED')
xlabel('Time (s)')
ylabel('Pitch: \theta (\circ)')

subplot(3,1,3)
plot(t,rad2deg(ang_imu_psi(1:length(t))),'k',t_mag,rad2deg(Z_mag),'c.',t,rad2deg(XX(3,1:length(t))),'b');
legend('IMU', 'Sensor','COMPUTED')
xlim([0 1600]);
xlabel('Time (s)')
ylabel('Yaw: \psi (\circ)')

print(figure(200),'Angles','-depsc');

figure(111)
subplot(3,1,1)
plot(t,rad2deg(XX(4,1:length(t))),'k')
xlim([0 1600]);
xlabel('Time (s)')
ylabel('Roll-rate-bias:b_p(\circ/s)')

subplot(3,1,2)
plot(t,rad2deg(XX(5,1:length(t))),'k')
xlim([0 1600]);
xlabel('Time (s)')
ylabel('Pitch-rate-bias:b_q(\circ/s)')

subplot(3,1,3)
plot(t,rad2deg(XX(6,1:length(t))),'k')
xlim([0 1600]);
xlabel('Time (s)')
ylabel('Yaw-rate-bias:b_r(\circ/s)')

print(figure(111),'bias','-depsc')
count = 0;
for kk=1:13:length(t)
    count = count+1;
    
    R_N_B = [cos(ang_imu_psi(kk))*cos(ang_imu_theta(kk))  sin(ang_imu_psi(kk))*cos(ang_imu_theta(kk))  -sin(ang_imu_theta(kk));
        cos(ang_imu_psi(kk))*sin(ang_imu_theta(kk))*sin(ang_imu_phi(kk))-sin(ang_imu_psi(kk))*cos(ang_imu_phi(kk))  cos(ang_imu_psi(kk))*cos(ang_imu_phi(kk))+sin(ang_imu_psi(kk))*sin(ang_imu_theta(kk))*sin(ang_imu_phi(kk)) cos(ang_imu_theta(kk))*sin(ang_imu_phi(kk))
        sin(ang_imu_psi(kk))*sin(ang_imu_phi(kk))+cos(ang_imu_psi(kk))*sin(ang_imu_theta(kk))*cos(ang_imu_phi(kk)) sin(ang_imu_psi(kk))*sin(ang_imu_theta(kk))*cos(ang_imu_phi(kk))-cos(ang_imu_psi(kk))*sin(ang_imu_phi(kk)) cos(ang_imu_theta(kk))*cos(ang_imu_phi(kk))];

    a_g(:,count) = R_N_B*[gps_acc_N(count);gps_acc_E(count);gps_acc_D(count)];
end

figure(100)
subplot(3,1,1)
plot(t,acc_gps_N(1:length(t)),'c.',t,acc_imu_x(1:length(t)),'k')%,t_gps,gps_acc_N,'r')%, t,acc_gps_x(1:length(t)),'m')
subplot(3,1,2)
plot(t,acc_imu_y(1:length(t)),'k',t_gps(1:length(a_g)),a_g(2,:),'r')%,t,acc_gps_E(1:length(t)),'m')
subplot(3,1,3)
plot(t,acc_imu_z(1:length(t)),'k',t_gps(1:length(a_g)),a_g(3,:),'r')%,t,acc_gps_D(1:length(t)),'m')

figure(101) % 3d trajectory
grid on
plot(t, imu_vel(1,1:length(t)),'k',t_gps,gps_vel_N,'b',t,vel_gps_N(1,1:length(t)),'c.');%  imu_vel(2,1:13:end),imu_vel(3,1:13:end),'b');                 
% hold on
% plot3(gps_vel_N,gps_vel_E,gps_vel_D,'k');