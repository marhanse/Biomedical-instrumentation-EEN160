function  im = k2im(k_space)

    % Dimensions of the k_space 
    [kx,ky,slice] = size(k_space);
    
    % Creating an 3D array with zeros 
    im = zeros(kx,ky,slice);

    % Calculating each image for each slice and adding to the zerro
    % 3D-array
    for i = 1:1:slice
        one_slice_image = ifftshift(ifft2(fftshift(k_space(:,:,i))));
        im(:,:,i) = one_slice_image;
    end

end 

