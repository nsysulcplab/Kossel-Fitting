function X = Kossels(x)
% 
% Lattice parameters: a, b, c, strain angle, and strain factors
% a = x(1); % first variable
% b = x(2); % second variable
% c=1;
% 
% % Define triclinic cell parameters
% alpha = x(3); 
% beta = x(4);
% gamma = x(5) ;

a = 0.6936					;		
b =0.9473   ;
c = 1;

% Define triclinic cell parameters
alpha = 90.0274; 
beta = 89.9161; 
gamma = 90.3722	 ;

% Define triclinic cell parameters
% alpha = 90; beta = 90; gamma = 90 ;
% Convert angles to radians
alpha_rad = deg2rad(alpha);
beta_rad = deg2rad(beta);
gamma_rad = deg2rad(gamma);
% Calculate transformation matrix
cos_alpha = cos(alpha_rad);
cos_beta = cos(beta_rad);
cos_gamma = cos(gamma_rad);
sin_gamma = sin(gamma_rad);
strain_ang_d = 90-gamma ;

% a=0.582411830099909;
% b=0.742887305361581	;
% strain_ang_d = 2.3457;
strain_ang = strain_ang_d*(pi/180); % unit in radiance; shear strain angle
% 
% a = 1/sqrt(2);
% b = 1;
factor_l = sqrt((a-0)^2 + (a*tan(strain_ang)+b/2)^2 + (c-c/2)^2)/sqrt((a-0)^2 + (b/2-0)^2 + (c-c/2)^2);
factor_s = sqrt((a-0)^2 + (a*tan(strain_ang)-b/2)^2 + (c/2-c)^2)/sqrt((a-0)^2 + (0-b/2)^2 + (c/2-c)^2);


%%
% Angle-dependent reflection spectrum of BP
%% Construct simulattion Kossel
% Angle-dependent reflection spectrum of BP
i=sqrt(-1);

lc= 525.50 ; % center wavelength of the normal reflectance  ***************************************
p = lc/(1.5179+10740/lc^2); % 1.5179+10740/𝜆^2 ,for 𝜆=lc 

n = 0.02; % difference of refractive index 
n_avg = 1.583; % 1.5179+10740/𝜆^2 ,for 𝜆=405nm
n_o = n_avg - n/2 ; % ordinary refractive index
n_e = n_avg + n/2 ; % extraordinary refractive index

alpha = 0.5*(n_e^2 - n_o^2); % high-frequency dielectric anisotropy
Beta = 0.5*(n_e^2 + n_o^2); % high-frequency dielectric anisotropy
sig = -1; % chirality/handedness
p_l = p*factor_l; % strained pitch 
p_s = p*factor_s; % shrunk pitch
N = 40; % number of layers
%P = 0.0011; %theoretical focal length

w = [7]; % bandwidth of the light source 
xc = 405; % center wavelength of the light source 

z = 765; % CCD pixel #

%% CCD pixel vs. angle 
pixel = 1:z ;
% P = -2e-9*pixel.^2 + 2e-6*pixel + 0.0007; % original one from Po-Chang
P = -5.04e-9*pixel.^2 + 5e-6*pixel + 0.0004; % modified one
% Real_theta = atan(((2*pixel)*1.2e-6)./P); % Pixels vs. angle in real space
Real_theta = atan(((2*pixel)*1.7e-6)./P);
p_p = (p*sqrt((a*p-0)^2 + (a*p*tan(strain_ang)+b*p/2)^2 + ...
    (c*p-c*p/2)^2)/sqrt((a*p-0)^2 + (b*p/2-0)^2 + (c*p-c*p/2)^2)* ...
    sqrt((a*p)^2 + (b*p)^2 + (c*p)^2)/sqrt(2.5*(c*p)^2))*cos(Real_theta); % angle-dependent pitch 
p_p_l = (p_l*sqrt((a*p_l-0)^2 + (a*p_l*tan(strain_ang)+b*p_l/2)^2 + ...
    (c*p_l-c*p_l/2)^2)/sqrt((a*p_l-0)^2 + (b*p_l/2-0)^2 + (c*p_l-c*p_l/2)^2)* ...
    sqrt((a*p_l)^2 + (b*p_l)^2 + (c*p_l)^2)/sqrt(2.5*(c*p_l)^2))*cos(Real_theta); % angle-dependent pitch 
p_p_s = (p_s*sqrt((a*p_s-0)^2 + (a*p_s*tan(strain_ang)+b*p_s/2)^2 + ...
    (c*p_s-c*p_s/2)^2)/sqrt((a*p_s-0)^2 + (b*p_s/2-0)^2 + (c*p_s-c*p_s/2)^2)* ...
    sqrt((a*p_s)^2 + (b*p_s)^2 + (c*p_s)^2)/sqrt(2.5*(c*p_s)^2))*cos(Real_theta); % angle-dependent pitch 
% p_p_l = p_l*cos(Real_theta);
% p_p_s = p_s*cos(Real_theta);
q = (2*pi)./p_p; % chirality
q_l = (2*pi)./p_p_l;
q_s = (2*pi)./p_p_s;
L = N*p_p; % angle-dependent thickness
L_l = N*p_p_l;
L_s = N*p_p_s;

%% Reflectance
x = 400:0.01:700; 
z1 = length(x);
k_0 = 2*pi./x; % vacuum propagation constant
k = sqrt(Beta*(k_0.^2)); % propagation constant
kapa = (k_0.^2).*alpha./(2.*k); % coupling constant
delta_k = 2*(repmat(k,z,1)-repmat(q.',1,z1));
delta_k_l = 2*(repmat(k,z,1)-repmat(q_l.',1,z1));
delta_k_s = 2*(repmat(k,z,1)-repmat(q_s.',1,z1));
s = sqrt(repmat(kapa.^2,z,1)-((delta_k).^2./4));
s_l = sqrt(repmat(kapa.^2,z,1)-((delta_k_l).^2./4));
s_s = sqrt(repmat(kapa.^2,z,1)-((delta_k_s).^2./4));
sL = s.*repmat(L.',1,z1);
sL_l = s_l.*repmat(L_l.',1,z1);
sL_s = s_s.*repmat(L_s.',1,z1);
r = 1i*repmat(kapa,z,1).*sinh(sL/2)./(s.*cosh(sL/2)-1i.*(delta_k/2).*sinh(sL/2));
r_l = 1i*repmat(kapa,z,1).*sinh(sL_l/2)./(s_l.*cosh(sL_l/2)-1i.*(delta_k_l/2).*sinh(sL_l/2));
r_s = 1i*repmat(kapa,z,1).*sinh(sL_s/2)./(s_s.*cosh(sL_s/2)-1i.*(delta_k_s/2).*sinh(sL_s/2));
y = abs(r).^2;
y_l = abs(r_l).^2;
y_s = abs(r_s).^2;

% figure;
% plot(x,y(1,:))

%%
%Monoutility light source
y2 = 1./(1+0.67801.*(2.*(x-xc)./w).^2+-1.71844.*(2.*(x-xc)./w).^4+3.49548.*(2.*(x-xc)./w).^6);
kossel = y2*(y.');
kossel_l = y2*(y_l.');
kossel_s = y2*(y_s.');
n_kossel = kossel/max(kossel);
n_kossel_l = kossel_l/max(kossel_l);
n_kossel_s = kossel_s/max(kossel_s);
mz = 1:z;

%figure(2);
%plot(mz,n_kossel_1)

%% (1 1 0) Lattice plane 

% Lattice constant 
% a=1;
% b=1.1;


r = z;% number of circles
t = 360;% around
R=linspace(0,1.79,r); % radius
theta = (0:t)/t*2*pi; %theta 0 to 2*pi

%%
% coordinate  
% coordinate 
A=[a,b/2,0;...
    a,0,c/2;...
    a,0,c/2;...
    a,b/2,c]; 
B=[0,b/2,c;...
    0,b,c/2;...
    0,b,c/2;...
    0,b/2,0];
C=[0,0,c/2;...
    0,b/2,0;...
    0,b/2,c;...
    0,0,c/2];

T = [1, 1*cos_gamma, 1*cos_beta;
     0, 1*sin_gamma, 1*(cos_alpha - cos_beta*cos_gamma)/sin_gamma;
     0, 0, 1*sqrt(1 - cos_alpha^2 - cos_beta^2 - cos_gamma^2 + 2*cos_alpha*cos_beta*cos_gamma)/sin_gamma];

% A=[a,(b/2 + a*tan(strain_ang)),0;...
%     a,(0 + a*tan(strain_ang)),c/2;...
%     a,(0 + a*tan(strain_ang)),c/2;...
%     a,(b/2 + a*tan(strain_ang)),c]; 
% B=[0,b/2,c;...
%     0,b,c/2;...
%     0,b,c/2;...
%     0,b/2,0];
% C=[0,0,c/2;...
%     0,b/2,0;...
%     0,b/2,c;...
%     0,0,c/2];

%% Coordinate vector 
% % Convert each matrix separately
A_tri =  (T*A')' ;
B_tri =  (T*B')' ;
C_tri =  (T*C')' ;


A = A_tri;
B = B_tri;
C = C_tri;
% Coordinate vector 
C_A=A-C; 
C_B=B-C;

N = cross(C_A,C_B,2); % Normal vector
N = N./sqrt(sum(N.^2,2));
ConeCent = N ; % Coordinates of the center of the circle 

VecA = cross(N,repmat([1 0 0],4,1),2); % vector B = n cross i 
VecA = VecA./sqrt(sum(VecA.^2,2)); %normalized A
VecB = cross(N,VecA,2); % vector B = n cross A 
VecB = VecB./sqrt(sum(VecB.^2,2)); % normalized B

%XX 361 by 765 by 4
%ConeCent 4 by 3 
%VecA  4 by 3
%R 1 by 765
%theta
CCr = reshape(ConeCent.',1,3,4); 
VAr = reshape(VecA.',1,3,4);
VBr = reshape(VecB.',1,3,4);
Rr = repmat(R,t+1,1,4); 
Thetar = repmat(theta.',1,r,4);

XX = CCr(1,1,:)+Rr.*cos(Thetar).*VAr(1,1,:)+Rr.*sin(Thetar).*VBr(1,1,:);
YY = CCr(1,2,:)+Rr.*cos(Thetar).*VAr(1,2,:)+Rr.*sin(Thetar).*VBr(1,2,:);
ZZ = CCr(1,3,:)+Rr.*cos(Thetar).*VAr(1,3,:)+Rr.*sin(Thetar).*VBr(1,3,:);
Nkossel = n_kossel .* ones(t+1,1);

Nkossel_l = n_kossel_l .* ones(t+1,1);
Nkossel_s = n_kossel_s .* ones(t+1,1);

%% Plot
fig = figure(3)

mesh(XX(:,:,1),YY(:,:,1),Nkossel_s )
hold on
mesh(XX(:,:,2),YY(:,:,2),Nkossel_l )
hold on
mesh(XX(:,:,3),YY(:,:,3),Nkossel_l )
hold on
mesh(XX(:,:,4),YY(:,:,4),Nkossel_s )
colormap gray
axis equal
grid off
axis off

set(gca,'xlim',[-0.6,0.6],'xtick', [-0.6:0.6:0.6]);
set(gca,'ylim',[-0.6,0.6],'ytick', [-0.6:0.6:0.6]);

view(0,90)

set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1167/100 875/100]);
print(fig, 'origin', '-dpng', '-r100');
% print(fig, 'origin', '-dpng')
close (3)
K_I = imread('origin.png');
K_ideal = double(rgb2gray(K_I));
width = 712;
K_ideal_c = imcrop(K_ideal, [249 67 width-1 width-1]);

% ang_I = 45 - rad2deg(atan(a/b)); ang_II = 45 + rad2deg(atan(a/b));
% ang_III = 315 - rad2deg(atan(a/b)); ang_IV = 315 + rad2deg(atan(a/b));  
% TwinK_I = imrotate(K_ideal_c, ang_I, 'crop');
% TwinK_II = imrotate(K_ideal_c, ang_II, 'crop');
% TwinK_III = imrotate(K_ideal_c, ang_III, 'crop');
% TwinK_IV = imrotate(K_ideal_c, ang_IV, 'crop');
% Sim_K = (TwinK_I/max(max(TwinK_I)) + TwinK_II/max(max(TwinK_II)) + TwinK_III/max(max(TwinK_III)) + TwinK_IV/max(max(TwinK_IV)))/4;
Sim_K = K_ideal_c/max(max(K_ideal_c));

%ini_ang = 45;
%ang_I = ini_ang - rad2deg(atan(a/b)); ang_II = ini_ang + rad2deg(atan(a/b));
%TwinK_I = imrotate(Sim_K, ang_I, 'crop');
%TwinK_II = imrotate(fliplr(Sim_K), ang_II, 'crop');
%Tsim_K = (TwinK_I/max(max(TwinK_I)) + TwinK_II/max(max(TwinK_II)))/2;

% rot_deg = linspace(0, 90, 360);
% rot_fit = zeros(1, 360);
% for r = 1:360;
%     roK = imrotate(Sim_K, rot_deg(r), 'crop');
%     roK_fl = fliplr(roK);
%     rot_fit(r) = corr2(roK, roK_fl);
% end
% rot_ind = find(rot_fit == max(rot_fit));
% rot = rot_deg(rot_ind);
Sim_K = imrotate(Sim_K, strain_ang_d, 'crop');


mask = ones(1001, 1001);
width = 711;
g_n = 459; g_v = linspace(501-(g_n-1)/2, 501+(g_n-1)/2, g_n); % the length of rectangular mask; g_n must be an odd value
h_n = 431; h_v = linspace(501-(h_n-1)/2, 501+(h_n-1)/2, h_n); % the length of rectangular mask; h_n must be an odd value
mask_x = 1:712; mask_x0 = 356;
mask_y = 1:712; mask_y0 = 356;
R0 = 363;  %r0 = 262;
[mask_X, mask_Y] = meshgrid(mask_x, mask_y);
cir = sqrt((mask_X - mask_x0).^2 + (mask_Y - mask_y0).^2);
for g = 1:g_n;
    for h = 1:h_n;
        mask(g_v(g), h_v(h)) = 0;
    end
end
mask = imrotate(mask, 45);
width_II = size(mask); width_II = width_II(1);
mask = imcrop(mask, [(width_II/2+1)-(width+1)/2, (width_II/2+1)-(width+1)/2, width, width]);
for e = mask_x;
    for f = mask_y;
        if cir(e, f) >= R0
            mask(e,f) = 0;
%         elseif cir(e,f) <= r0
%             mask(e,f) = 0;
%         else
%             mask(e,f) = 1;
        end    
    end
end

Sim_K = Sim_K.*mask;

Exp_K_o = rgb2gray(imread('Exp-Kossel.tif'));
Exp_K_c = imcrop(Exp_K_o, [114 38 739 739]);
Exp_Kr = imresize(Exp_K_c, 712/740);
Exp_K = double(Exp_Kr).*mask;
Exp_K = Exp_K/(max(max(Exp_K)));

figure(1)    
imshow(Sim_K)
figure(2)  
imshow(Exp_K)

X = sum(sum(abs(Sim_K - Exp_K))); % objectivefunction




