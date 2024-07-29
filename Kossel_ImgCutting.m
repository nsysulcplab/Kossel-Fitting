clear
close all

%% Read in the Kossel image captured by CMOS
K_I = imread('Org-Kossel.tif');
K_ideal = rgb2gray(K_I);

%% Find the centroid of Kossel diagram using weighted average
kossel_width = 400; % Determine the pixel dimensions of the Kossel field within the image.
sim_size = 800; % Determine the pixel size of the simulated Kossel diagram

[height, width] = size(K_ideal);
width = min(height, width) - 1;
K_ideal_c = imcrop(K_ideal, [1 1 width width]);
K_ideal_f = imgaussfilt(K_ideal_c, 3);
K_ideal_d = double(K_ideal_f); 
IwKossel_x = zeros(width+1, width+1);
IwKossel_y = zeros(width+1, width+1);
Ix = linspace(1, width+1, width+1);
Iy = linspace(1, width+1, width+1);
for u = 1:width+1
    IwKossel_x(:,u) = K_ideal_d(:,u).* Ix.';
    IwKossel_y(u,:) = K_ideal_d(u,:).* Iy;
end
I_ss = zeros(1,2);
I_ss(1) = round(sum(sum(IwKossel_x))/sum(sum(K_ideal_d))); % Calculate the x-axis of the centroid
I_ss(2) = round(sum(sum(IwKossel_y))/sum(sum(K_ideal_d))); % Calculate the y-axis of the centroid

%%  Confirm the cropped region is the Kossel diagram field of view
% Because the brightness of the Kossel diagram is uneven, this part relies on the previously found centroid position and uses the correlation coefficient to obtain a new centroid position. For most crystal structures (e.g., cubic, orthorhombic, tetragonal), the Kossel diagram displays line symmetry, with the correlation coefficient calculated using the flipped image.
K_ideal_new = zeros(kossel_width, kossel_width);
Ihalf_width_II = floor(kossel_width/2);
Itarget_w = 150 + 1;
Itarget_x = I_ss(2) - floor(Itarget_w/2): I_ss(2) + floor(Itarget_w/2); 
Itarget_y = I_ss(1) - floor(Itarget_w/2): I_ss(1) + floor(Itarget_w/2); 
IFit = zeros(Itarget_w, Itarget_w);
for m = 1: Itarget_w
    for n = 1:Itarget_w
     K_ideal_new = imcrop(K_ideal_d, [Itarget_x(m)-Ihalf_width_II, Itarget_y(n)-Ihalf_width_II, kossel_width, kossel_width]);
     K_ideal_flip = fliplr(flipud(K_ideal_new));
     IFit(n, m) = corr2(K_ideal_new, K_ideal_flip); 
    end
end
ICx_ind = find(IFit == max(max(IFit)));
ICx_new = Itarget_x(ceil(ICx_ind/Itarget_w));
ICy_ind = find(IFit == max(max(IFit)));
ICy_new = Itarget_y(rem(ICx_ind,Itarget_w));

% figure(1)
% imagesc(IFit); title('IFit')

%% Recrop the image based on the calculated centroid position
K_ideal_new = imcrop(K_ideal_d, [ICx_new-Ihalf_width_II, ICy_new-Ihalf_width_II, kossel_width, kossel_width]); 

% figure(2)
% imshow(K_ideal_new, [0 max(max(K_ideal_new))]); hold on; plot(Ihalf_width_II + 1, Ihalf_width_II + 1, 'r*'); title('Centered')

%% Rotate the image based on calculating the correlation coefficient.
rot_deg = linspace(0, 90, 360);
rot_fit = zeros(1, 360);
for r = 1:360;
    roK = imrotate(K_ideal_new, rot_deg(r), 'crop');
    roK_fl = fliplr(roK);
    rot_fit(r) = corr2(roK, roK_fl);
end
rot_ind = find(rot_fit == max(rot_fit));
rot = rot_deg(rot_ind); 
% rot=rot-45; % Adjusting the Kossel orientation: Rotate the lattice plane [xxx] to the vertical direction.
roK = imrotate(K_ideal_new, rot , 'crop');

% figure(3)
% imshow(roK, [8 25]); title('Rotated') % Set the grayscale range of the image

mag = sim_size/(kossel_width + 1); 
re_K = imresize(roK, mag);
Exp_K = re_K/max(max(re_K));

% figure(4)
% imshow(Exp_K, [0 1]); title('Normalized')

%% Save the cropped and rotated image
figure();
imshow(Exp_K);
imwrite(Exp_K, 'Exp-Kossel.tiff' );