function im = SENSE(kspace_subsampled, coil_sensitivities)

    load MRIdata.mat 

    % Folded for coil 1
    k_coil1 = kspace_subsampled(:,:,:,1); 
    folded_im1 = k2im(k_coil1);

    % Folded for coil 2 
    k_coil2 = kspace_subsampled(:,:,:,2); 
    folded_im2 = k2im(k_coil2);
    
    % Size of the full coil_sensitivities 
    [FOVx,FOVy,slice,coil] = size(coil_sensitivities);
    
    % The zero matrix were we will build the reconstructed image 
    im = zeros(FOVx,FOVy,slice); 
    
    % One loop for each slices 
    for s = 1:slice
       
        % One loop for each row in the folded image 
        for row = 1:size(folded_im1,1)
            
            % One more loop for each column in the row 
            for col = 1:FOVy
                
                % Getting the position for the corresponding point in 
                % Coil_sensitivity
                posA = (FOVx/2)-(FOVx / 4) + row;
    
                % Getting the other position for the corresponding point in 
                % Coil_sensitivity
    
                % Mod is to wraparound so we dont go out of range
                if mod((FOVx / 2) + (FOVx / 4) + row,FOVx) == 0
                    posB = 1;
                else
                    posB = mod((FOVx / 2) + (FOVx / 4) + row,FOVx);
                end
                
                % The folded pixel values in coil1 and coil2 
                P1 = folded_im1(row,col,s);
                P2 = folded_im2(row,col,s);
    
                % Picking out the pixel in the coil sens. for posA and posB
                % Coil1
                S1A = coil_sensitivities(posA,col,s,1);
                S2A = coil_sensitivities(posB,col,s,1);
                
                % Picking out the pixel in the coil sens. for posA and posB
                % Coil2
                S1B = coil_sensitivities(posA,col,s,2);
                S2B = coil_sensitivities(posB,col,s,2);
    
                % Creating the matrix with sensitivitys 
                sensitivity_matrix = [S1A S2A; S1B S2B];
    
                % Creating one with pixelvalues from folded 
                pixel_matrix = [P1; P2];
    
                % Reconstructed points is given by solving the matrix equation
                reconstructed_points = sensitivity_matrix\pixel_matrix;
                
                % Building the reconstructed image pixel by pixel 
                im(posA, col,s) = reconstructed_points(1);
                im(posB, col,s) = reconstructed_points(2);
    
            end
        end
    end


end

