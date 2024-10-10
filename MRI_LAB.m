%% P2 
load MRIdata.mat

X_dim = nifti_info.ImageSize(1);
Y_dim = nifti_info.ImageSize(2);
Z_dim = nifti_info.ImageSize(3);

% Field of view in the different dimensions 
FOVx = X_dim * nifti_info.PixelDimensions(1);
FOVy = Y_dim * nifti_info.PixelDimensions(2);
FOVz = Z_dim * nifti_info.PixelDimensions(3);


disp("FOVx equals " + FOVx + " mm ")
disp("FOVy equals " + floor(FOVy) + " mm ")
disp("FOVz equals " + floor(FOVz) + " mm ")

%% P3 

% Look at the function created: save_nifti_image.m

%% A1, A2 and A3 Basic MR reconstruction

load MRIdata.mat

% MR image reconstruction (see function for detail) 
reconstructed_image = k2im(kspace_full);

% saves it as a magnitude image using the Nifti file format
save_nifti_image(reconstructed_image,'A3',nifti_info)

% Plots one slice 
imshow(abs(reconstructed_image(:,:,100)),[])
title('Basic MR reconstruction')

%% B1 

load MRIdata.mat

% K space data from first coil 
k_coil1 = kspace_coils(:,:,:,1); 

% MR image reconstruction (see function for detail) 
reconstructed_image_1 = k2im(k_coil1);

% saves it as a magnitude image using the Nifti file format
save_nifti_image(reconstructed_image_1.*255,'B1_coil1',nifti_info)

% K space data from second coil 
k_coil2 = kspace_coils(:,:,:,2);

% MR image reconstruction (see function for detail) 
reconstructed_image_2 = k2im(k_coil2);

% saves it as a magnitude image using the Nifti file format
% Multiplies with 255 to be able to view in 3D-slicer 
save_nifti_image(reconstructed_image_2.*255,'B1_coil2',nifti_info)

%  compare one slice for each coil-image pair 

% Plot the pair for coil1
subplot(2,2,1)
imshow(real(abs(reconstructed_image_1(:,:,150))),[])
title('Recontructed Image from coil 1')

subplot(2,2,2)
imshow(coil_sensitivities(:,:,150,1),[])
title('Image from coil sensitivity 1')

% And for pair coil2 
subplot(2,2,3)
imshow(real(abs(reconstructed_image_2(:,:,150))),[])
title('Recontructed Image from coil 2')

subplot(2,2,4)
imshow(coil_sensitivities(:,:,150,2),[])
title('Image from coil sensitivity 2')

%% B2 Combine the different coil images into a single image using a  'sum-of-squares reconstruction':

% Step one: Square the pixels in the images
[kx,ky,slice] = size(reconstructed_image_1);

% Step two: Sum the squared images 
sum_squared_images_recon = zeros(kx,ky,slice);

% Step two: Sum the squared images
% Step three: Square root 
for s = 1:slice

    % Square of each pixel for each slice 
    squared_slice1 = reconstructed_image_1(:,:,s).^2;
    squared_slice2 = reconstructed_image_2(:,:,s).^2;

    % Sum the squared images 
    sum_slices = (squared_slice1 + squared_slice2);

    % Square root of the squared images 
    sum_squared_images_recon(:,:,s) = sqrt(sum_slices);

end

% Multiplies with 255 to be able to view in 3D-slicer 
save_nifti_image(sum_squared_images_recon.*255,'B2 Combine Coils',nifti_info)

imshow(abs(sum_squared_images_recon(:,:,100)),[])
title('Sum of squares')


%% B3 Combine the different coil images into a single image using a modified version 

% Step one: Square the pixels in the images
[kx,ky,slice] = size(reconstructed_image_1);

% Step two: Sum the squared images 
sum_images = zeros(kx,ky,slice);

for s = 1:slice
    coil_sense1_in = 1./coil_sensitivities(:,:,s,1);
    coil_sense2_in = 1./coil_sensitivities(:,:,s,2);

    % Helps to remove the variations 
    slice1 = (reconstructed_image_1(:,:,s).^2).*coil_sense1_in;
    slice2 = (reconstructed_image_2(:,:,s).^2).*coil_sense2_in;
    sum = sqrt(slice1 + slice2);
    sum_images(:,:,s) = sum;
end

save_nifti_image(sum_images,'B3',nifti_info)

subplot(1,2,1)
imshow(real(abs(sum_images(:,:,150))),[])
title('Modified sum of squares')

subplot(1,2,2)
imshow(abs(sum_squared_images_recon(:,:,150)),[])
title('Sum of squares')

%% C1 Reconstruct images from data from each coil separately and save them as Nifti files. 
% Compare them with the images obtained from part B

load MRIdata.mat

% Lets edit the meta data in niftiinfo 
new_nifti_info = nifti_info; % Copy of nifti_info
new_nifti_info.ImageSize = [96, 129, 192]; % size for folded 

% K space data from first coil 
k_coil1 = kspace_subsampled(:,:,:,1); 

% MR image reconstruction (see function for detail) 
folded1 = k2im(k_coil1);

% saves it as a magnitude image using the Nifti file format
save_nifti_image(folded1.*255,'C1_folded1',new_nifti_info)

% K space data from second coil 
k_coil2 = kspace_subsampled(:,:,:,2);

% MR image reconstruction (see function for detail) 
folded2 = k2im(k_coil2);

% saves it as a magnitude image using the Nifti file format
save_nifti_image(folded2.*255,'C1_folded2',new_nifti_info)

% Coil 1 
subplot(1,2,1)
imshow(real(abs(folded1(:,:,150))),[])
title('Folded Image from coil 1')

% Coil 2
subplot(1,2,2)
imshow(real(abs(folded2(:,:,150))),[])
title('Folded Image from coil 2')


%% C2 reconstruction using sense 

load MRIdata.mat

% Reconstruction using SENSE 
reconstructed_sense = SENSE(kspace_subsampled,coil_sensitivities);

subplot(1,2,1)
imshow(abs(reconstructed_sense(:,:,150)),[])
title('Reconstruction using sense')


% To plot this you need to run section A1-3 first 
subplot(1,2,2)
imshow(abs(reconstructed_image(:,:,150)),[])
title('Reconstruction of kspacefull')

%save_nifti_image(reconstructed_sense.*255,'C2 Sense',nifti_info)










