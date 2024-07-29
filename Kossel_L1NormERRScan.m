clear
close all

%% Read the experimental figure
Exp_K = im2double(imread('Exp-Kossel.tiff'));

%% Design two graphical masks to filter out noise
% When creating the mask, confirm that the parameters only target the noise and do not disturb the main Kossel lines.
[height, width] = size(Exp_K);

% Rectangular mask
nExp_K = zeros(height, height);
mask_h = height/2;
mask_w = height/2;

mask = ones(height, height);
g_n = mask_h+1; 
g_v = linspace(mask_h-(g_n-1)/2, mask_h+(g_n-1)/2, g_n); % the length of rectangular mask; g_n must be an odd value
h_n = mask_w+1; 
h_v = linspace(mask_w-(h_n-1)/2, mask_w+(h_n-1)/2, h_n); % the length of rectangular mask; h_n must be an odd value
mask_x = 1:height; mask_x0 = height/2;
mask_y = 1:height; mask_y0 = height/2;

[mask_X, mask_Y] = meshgrid(mask_x, mask_y);
cir = sqrt((mask_X - mask_x0).^2 + (mask_Y - mask_y0).^2);
for g = 1:g_n
    for h = 1:h_n
        mask(g_v(g), h_v(h)) = 0;
    end
end

mask = imrotate(mask, 45.5);
width_II = size(mask); width_II = width_II(1);
mask = imcrop(mask, [(width_II/2+1)-(height)/2, (width_II/2+1)-(height)/2, height-1, height-1]);

% Circular mask
R0 = (height/2)*1.02; 

for e = mask_x
    for f = mask_y
        if cir(e, f) >= R0
            mask(e,f) = 0;
        end    
    end
end

nExp_K = Exp_K.*mask;

% figure
% imshow(nExp_K);

%% Construct simulated Kossel diagram 
% Angle-dependent reflection spectrum of BP
i=sqrt(-1);

lc= 525.5 ; % center wavelength of the normal reflectance [Modify according to experimental data.]
p = lc/(1.5179+10740/lc^2); % liquid crystal dispersion [Modify according to experimental data.]

n = 0.02; % difference of refractive index  [Modify according to experimental data.]
n_avg = 1.583; % [Modify according to experimental data.]
n_o = n_avg - n/2 ; % ordinary refractive index 
n_e = n_avg + n/2 ; % extraordinary refractive index

alpha = 0.5*(n_e^2 - n_o^2); % high-frequency dielectric anisotropy
Beta = 0.5*(n_e^2 + n_o^2); % high-frequency dielectric anisotropy
sig = -1; % chirality/handedness

N = 40; % number of layers [Modify according to experimental data.]

w = [7]; % bandwidth of the light source [Modify according to experimental data.]
xc = 405; % center wavelength of the light source [Modify according to experimental data.]

z = 765; % CCD pixel [Modify according to experimental data.]

%% Calibrate imaging system using a reflective amplitude grating
pixel = 1:z ;
P = -5.04e-9*pixel.^2 + 5e-6*pixel + 0.0004; % [Modify according to experimental data.]
Real_theta = atan(((2*pixel)*1.7e-6)./P); % [Modify according to pixel size.]

%% Set the scanning range for a/c and b/c, and divide the range into m parts.
num_scan = 20; % the m numberb [Modify according to experimental data.]
a_space = linspace(0.67, 0.71, num_scan); % the bounds of a/c [Modify according to experimental data.]
b_space = linspace(0.93, 0.97, num_scan); % the bounds of b/c [Modify according to experimental data.]

[A_dev, B_dev] = meshgrid(a_space, b_space);
X_dev = zeros(num_scan, num_scan);

x = 400:0.01:700; % spectrum range and resolution
c = 1; % c is in the BCC symmetry and normalized to 1
r1 = z;% number of circles
t = 360;% around
R=linspace(0,1.79,r1); % radius
theta = (0:t)/t*2*pi; %theta 0 to 2*pi

for aa = 1:num_scan
    for bb = 1:num_scan

%% (Q1, Q2, and Q3) on the lattice plane used to calculate the plane normal.        
Q1=[a_space(aa),b_space(bb)/2,0;...
    a_space(aa),0,c/2;...
    a_space(aa),0,c/2;...
    a_space(aa),b_space(bb)/2,c]; 
Q2=[0,b_space(bb)/2,c;...
    0,b_space(bb),c/2;...
    0,b_space(bb),c/2;...
    0,b_space(bb)/2,0];
Q3=[0,0,c/2;...
    0,b_space(bb)/2,0;...
    0,b_space(bb)/2,c;...
    0,0,c/2];

% Coordinate vector 
C_A=Q1-Q3; 
C_B=Q2-Q3;

N1 = cross(C_A,C_B,2); % Normal vector
N1 = N1./sqrt(sum(N1.^2,2));
ConeCent = N1 ; % Coordinates of the center of the circle 

%% Reflectance
a = a_space(aa);
b = b_space(bb);
p_p = (p*sqrt((a*p-0)^2 + (b*p/2)^2 + ...
    (c*p-c*p/2)^2)/sqrt((a*p-0)^2 + (b*p/2-0)^2 + (c*p-c*p/2)^2)* ...
    sqrt((a*p)^2 + (b*p)^2 + (c*p)^2)/sqrt(2.5*(c*p)^2))*cos(Real_theta); % angle-dependent pitch 
q = (2*pi)./p_p; % chirality
L = N*p_p; % angle-dependent thickness

z1 = length(x);
k_0 = 2*pi./x; % vacuum propagation constant
k = sqrt(Beta*(k_0.^2)); % propagation constant
kapa = (k_0.^2).*alpha./(2.*k); % coupling constant
delta_k = 2*(repmat(k,z,1)-repmat(q.',1,z1));
s = sqrt(repmat(kapa.^2,z,1)-((delta_k).^2./4));
sL = s.*repmat(L.',1,z1);
r = 1i*repmat(kapa,z,1).*sinh(sL/2)./(s.*cosh(sL/2)-1i.*(delta_k/2).*sinh(sL/2));
y = abs(r).^2;

%% The narrowband probe light source for observing the Kossel diagram
% [Modify according to experimental data.]
% Make sure it corresponds to the previously defined spectrum range and resolution
y2 = 1./(1+0.67801.*(2.*(x-xc)./w).^2+-1.71844.*(2.*(x-xc)./w).^4+3.49548.*(2.*(x-xc)./w).^6);

kossel_1 = y2*(y.');
n_kossel = kossel_1/max(kossel_1);
mz = 1:z;

VecA = cross(N1,repmat([1 0 0],4,1),2); 
VecA = VecA./sqrt(sum(VecA.^2,2)); 
VecB = cross(N1,VecA,2); 
VecB = VecB./sqrt(sum(VecB.^2,2));

CCr = reshape(ConeCent.',1,3,4); 
VAr = reshape(VecA.',1,3,4);
VBr = reshape(VecB.',1,3,4);
Rr = repmat(R,t+1,1,4); 
Thetar = repmat(theta.',1,r1,4);

XX = CCr(1,1,:)+Rr.*cos(Thetar).*VAr(1,1,:)+Rr.*sin(Thetar).*VBr(1,1,:);
YY = CCr(1,2,:)+Rr.*cos(Thetar).*VAr(1,2,:)+Rr.*sin(Thetar).*VBr(1,2,:);
ZZ = CCr(1,3,:)+Rr.*cos(Thetar).*VAr(1,3,:)+Rr.*sin(Thetar).*VBr(1,3,:);
Nkossel = n_kossel .* ones(t+1,1);

%% Plot
fig = figure(1);

mesh(XX(:,:,1),YY(:,:,1),Nkossel )
hold on
mesh(XX(:,:,2),YY(:,:,2),Nkossel )
hold on
mesh(XX(:,:,3),YY(:,:,3),Nkossel )
hold on
mesh(XX(:,:,4),YY(:,:,4),Nkossel )
colormap gray
axis equal
grid off
axis off

set(gca,'xlim',[-0.6,0.6],'xtick', [-0.6:0.6:0.6]);
set(gca,'ylim',[-0.6,0.6],'ytick', [-0.6:0.6:0.6]);

view(0,90)
imshow(Exp_K, 'Border', 'tight'); 
set(gca, 'LooseInset', get(gca, 'TightInset')); 
print(fig, 'origin', '-dpng')
sK_I = imread('origin.png');
sK_ideal = double(rgb2gray(sK_I));
[sK_height, sK_width] = size(sK_ideal);
crop_x = round((sK_height-height)/2);
crop_y = round((sK_width-width)/2);
sK_ideal_c = imcrop(sK_ideal, [crop_x crop_y height-1 height-1]);   
Sim_K = sK_ideal_c/max(max(sK_ideal_c)); 

Sim_K = Sim_K.*mask;

% Cost function is the L1 norm of the difference between simulated and measured Kossel diagrams.
X_dev(bb,aa) =  sum(sum(abs(Sim_K - nExp_K))); 

    end
end

%% Kossel-diagram fitting for extracting lattice constants a and b
figure
surf(A_dev, B_dev, X_dev)
shading interp

tar_ind = find(X_dev == min(min(X_dev)));

% The optimum solution for a and b (relative to c)
a2c = A_dev(tar_ind);
b2c = B_dev(tar_ind);

fprintf('a/c=%f\n', a2c);
fprintf('b/c=%f\n', b2c);