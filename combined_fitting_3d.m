% Take Pitch -> x, Roll -> y, Yaw -> z
% Data: 1-3 columns are voltage, 4-6 are Magnetic Field

d_p = xlsread('qual_model.xlsx','Pitch_Raw');
d_p = sortrows(d_p, 4); % Sort using 4th column in ascending
dim_p = size(d_p);

d_r = xlsread('qual_model.xlsx','Roll_Raw');
d_r = sortrows(d_r, 4); % Sort using 4th column in ascending
dim_r = size(d_r);

d_y = xlsread('qual_model.xlsx','Yaw_Raw');
d_y = sortrows(d_y, 4); % Sort using 4th column in ascending
dim_y = size(d_y);

% Make Linear Regression Individually for each axis to get intial guess
% Make intitial guess. We guess that the cross axis terms are 0 and use
% linear fitting for individual axis to get sensitivity (slope) and offset
v0 = zeros(1, 12);
[r_p, m_p, c_p] = regression(d_p(:,1), d_p(:,4), 'one');
[r_r, m_r, c_r] = regression(d_r(:,2), d_r(:,4), 'one');
[r_y, m_y, c_y] = regression(d_y(:,3), d_y(:,4), 'one');
v0(1) = m_p;
v0(5) = m_r;
v0(9) = m_y;
v0(10) = -c_p/m_p;  % Offsets
v0(11) = -c_r/m_r;
v0(12) = -c_y/m_y;

% Add more data with Fields having values in all axes (optional).
d_com = xlsread('qual_model.xlsx','Combined_Raw');
dim_com = size(d_com);

sz = dim_r(1) + dim_p(1) + dim_y(1) + dim_com(1);
d_f = zeros(sz, 6); %  Full data of V_Measured, B_Cal

for k = 1:dim_p(1)
    d_f(k, 1:3) = d_p(k, 1:3);
    d_f(k, 4) = d_p(k, 4);
end
for k = dim_p(1)+1 : dim_r(1)+dim_p(1)
    d_f(k, 1:3) = d_r(k-dim_p(1), 1:3);
    d_f(k, 5) = d_r(k-dim_p(1), 4);
end
for k = dim_r(1)+dim_p(1)+1 : dim_r(1)+dim_p(1)+dim_y(1)
    d_f(k, 1:3) = d_y(k-(dim_r(1)+dim_p(1)), 1:3);
    d_f(k, 6) = d_y(k-(dim_r(1)+dim_p(1)), 4);
end
for k = dim_r(1)+dim_p(1)+dim_y(1)+1 : sz
    d_f(k, 1:6) = d_com(k-(dim_r(1)+dim_p(1)+dim_y(1)), 1:6);
end

% Define the function with 12 parameters
fun = @(v,vdata)[v(1)*(vdata(:,1)-v(10)) + v(2)*(vdata(:,2)-v(11)) + v(3)*(vdata(:,3)-v(12)), ...
                 v(4)*(vdata(:,1)-v(10)) + v(5)*(vdata(:,2)-v(11)) + v(6)*(vdata(:,3)-v(12)), ...
                 v(7)*(vdata(:,1)-v(10)) + v(8)*(vdata(:,2)-v(11)) + v(9)*(vdata(:,3)-v(12))];
% Intital guess is based on linear single axis fitting as defined above
% optimset('lsqcurvefit')
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', ...
    'TolFun',1e-12,'MaxIter',1000);
lb = [];
ub = [];
v = lsqcurvefit(fun,v0,d_f(:,[1 2 3]),d_f(:,[4 5 6]),lb,ub,options);
% Inputs to lsqcurvefit: function_handle, intial_guess,
% independent_variable (3D), dependent_variable (3D), parameter_upper_bound
% , parameter_lower_bound, special_options

% Substitute obtained parameters into voltage dataset and find residuals
d_fit = [d_f, zeros(sz,1)];
for k = 1:sz
    d_fit(k,4) = v(1)*(d_f(k,1)-v(10))+v(2)*(d_f(k,2)-v(11))+v(3)*(d_f(k,3)-v(12));
    d_fit(k,5) = v(4)*(d_f(k,1)-v(10))+v(5)*(d_f(k,2)-v(11))+v(6)*(d_f(k,3)-v(12));
    d_fit(k,6) = v(7)*(d_f(k,1)-v(10))+v(8)*(d_f(k,2)-v(11))+v(9)*(d_f(k,3)-v(12));
    d_fit(k,7) = sqrt(sum((d_f(k,[4 5 6])-d_fit(k,[4 5 6])).^2));
end

figure, subplot(3, 1, 1),
stem(d_fit(:,7));
xlabel('Bit Counts'); ylabel('Residuals (nT)'); 
grid on;
title('Interpolated Data Difference For 3D Linear Fitting');

% Substitute original parameters into voltage dataset and find residuals
d_fit2 = [d_f, zeros(sz,1)];
for k = 1:sz    
    d_fit2(k,4) = v0(1)*(d_f(k,1)-v0(10))+v0(2)*(d_f(k,2)-v0(11))+v0(3)*(d_f(k,3)-v0(12));
    d_fit2(k,5) = v0(4)*(d_f(k,1)-v0(10))+v0(5)*(d_f(k,2)-v0(11))+v0(6)*(d_f(k,3)-v0(12));
    d_fit2(k,6) = v0(7)*(d_f(k,1)-v0(10))+v0(8)*(d_f(k,2)-v0(11))+v0(9)*(d_f(k,3)-v0(12));
    d_fit2(k,7) = sqrt(sum((d_f(k,[4 5 6])-d_fit2(k,[4 5 6])).^2));
end

subplot(3, 1, 2),
stem(d_fit2(:,7));
xlabel('Bit Counts'); ylabel('Residuals Single Axis Fit (nT)'); 
title('Interpolated Data Difference For Single Axis Linear Fitting');
grid on;

subplot(3, 1, 3),
plot(d_fit2(:,7)-d_fit(:,7));
xlabel('Bit Counts'); ylabel('Residuals Difference (nT)'); 
title('Difference in Residuals for Both Fits');
grid on;

%%
% Refer Section 9.5 of Gunter Musmann for the breakdown of the formulas
E = zeros(3,3); % Sensitivity Matrix
O = E;          % Othogonalization matrix
M_2 = [v(1:3); v(4:6); v(7:9)]; % Calibration matrix

% Sensitivity Matrix (Purely Diagonal)
E(1,1) = sqrt(v(1)^2 + v(4)^2 + v(7)^2);
E(2,2) = sqrt(v(2)^2 + v(5)^2 + v(8)^2);
E(3,3) = sqrt(v(3)^2 + v(6)^2 + v(9)^2);

% Othogonalization matrix
% We assume that sensor x and facility X is aligned, sensor y is in the
% same plane as facility X-Y plane.
cos_xy = (v(1)*v(2) + v(4)*v(5) + v(7)*v(8))/(E(1,1)*E(2,2));
cos_yz = (v(2)*v(3) + v(5)*v(6) + v(8)*v(9))/(E(3,3)*E(2,2));
cos_xz = (v(1)*v(3) + v(4)*v(6) + v(7)*v(9))/(E(1,1)*E(3,3));
sin_xy = sqrt(1-cos_xy^2);

O(1,1) = 1;
O(1,2) = cos_xy;
O(1,3) = cos_xz;
O(2,2) = sin_xy;
O(2,3) = (cos_yz - cos_xy*cos_xz)/sin_xy;
O(3,3) = sqrt(1-cos_xz^2 - ((cos_yz - cos_xy*cos_xz)/sin_xy)^2);
% Other elements are 0

% Transformation Matrix (Euler angles)
% D_Euler = M_2*inv(E)*inv(O);
D_Euler = M_2/(O*E);                % Faster than commented statement
a_y = asin(D_Euler(3,1));           % y Euler Angle (b/w coil and sensor)
a_x = asin(D_Euler(3,2)/cos(a_y));
a_z = asin(-D_Euler(2,1)/cos(a_y));

% _d subscript denotes angle in degrees
a_y_d = a_y*180/pi;         % y Euler Angle (b/w coil and sensor) in degree
a_x_d = a_x*180/pi;
a_z_d = a_z*180/pi;

xy_d = acos(cos_xy)*180/pi; % xy orthogonal angle in degrees
xz_d = acos(cos_xz)*180/pi;
yz_d = acos(cos_yz)*180/pi;

gamma = asin(cos_xz);
eta_d = asin(O(2,3)/cos(gamma))*180/pi; 
beta_d = asin(cos_xy)*180/pi; 
gamma_d = gamma*180/pi;

offsets = v(10:12); % Offsets
final_cal = O*E;    % Final Calibration matrix