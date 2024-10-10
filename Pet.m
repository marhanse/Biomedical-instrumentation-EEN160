%% A1 Forward matrix 

load PETdata.mat 

% Calculated forward projected matrix 
% For the script check the function compute_forward_matrix
A = compute_forward_matrix(1:180,87,87,87);


%% A2 NAIVE RECONSTRUCTION 

load PETdata.mat

% Function to calculate the reconstructed image using NAIVE 
F_NAIVE = NAIVE(A,G,[87,87,10]);

%% A2 Plot 

% Comutes an inverted grayscale image of the reconstructed 
naive_reconstructed = imcomplement(F(:,:,10));

imshow(naive_reconstructed,[])
title('A2 PET Reconstruction using NAIVE reconstruction')

%% B2 ITERATIVE reconstruction for different itterations 

load PETdata.mat
load A.mat

p = 0;
for iter = [5,10,20,50, 70, 100]
    p = p + 1;
    [F_iter, error] = MLEM(A,G,[87,87,10],iter); 
    subplot(2,3,p)
    imshow(imcomplement(F_iter(:,:,10)),[])
    title("MLEM reconstruction " + iter + " Iterations")
end

%% B2 plot the residuels vs iterations 

% Plotting only one plot for residuals vs iterations and marking a line on the
% different number of iterations to see how much it converged. 

plot(0:69,error(2:71))
xline([5 10 20 50],'-',{'5 iterations','10 iterations','20 iterations', '50 iterations'})
title('Sum-of-squared residuals vs iteration number')
ylabel('Sum-of-squared Residuals')
xlabel('Iteration number')




%% C1 Using NAIVE to reconstrcut G-noisy 

load PETdata.mat
load F.mat

% Function to calculate the reconstructed image using NAIVE 
% Applied on sinogram G_noisy 

subplot(1,2,1)
F_noisy = NAIVE(A,G_noisy,[87,87,10]);
imshow(imcomplement(F_noisy(:,:,10)),[])
title('Reconstruction of sinogram G_noisy, slice 10, image')

% Function to calculate the reconstructed image using NAIVE 
% Applied sinogram G

subplot(1,2,2)
imshow(imcomplement(F_NAIVE(:,:,10)),[])  
title('Reconstruction of sinogram G ,slice 10, image')

%% C2 Using FBP to reconstruct G-noisy 

load PETdata.mat
load A.mat


subplot(1,2,1)
F_FBP1 = FBP(G_noisy,0:179,[87,87,26]);
imshow(imcomplement(F_FBP1(:,:,10)),[]);
title('FBP of sinogram G noicy ,slice 10, image')

subplot(1,2,2)
F_FBP = FBP(G,0:179,[87,87,26]);
imshow(imcomplement(F_FBP(:,:,10)),[]);
title('FBP of sinogram G ,slice 10, image')

%% C3 Subplotting reconstruction of G-noisy with different numbers of iterations 


subplot(2,3,1)
[F, ~] = MLEM(A,G_noisy,[87,87,10],10); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 10 Iterations')

subplot(2,3,2)
[F, ~] = MLEM(A,G_noisy,[87,87,10],30); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 30 Iterations')

subplot(2,3,3)
[F, ~] = MLEM(A,G_noisy,[87,87,10],50); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 50 Iterations')

subplot(2,3,4)
[F, ~] = MLEM(A,G_noisy,[87,87,10],70); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 70 Iterations')

subplot(2,3,5)
[F, ~] = MLEM(A,G_noisy,[87,87,10],90); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 90 Iterations')

subplot(2,3,6)
[F, ~] = MLEM(A,G_noisy,[87,87,10],100); 
imshow(imcomplement(F(:,:,10)),[])
title('MLEM reconstruction of G noisy 100 Iterations')

%% C4 Find optimal number of iterations 

load PETdata.mat

% This was used to find the optimal number of iterations 
% See the report for a more detailed explanation 

% F = ones(87,87,10);    
% r_error = zeros(200+1, 1);
% for slice = 10
%     g = G_noisy(:,:,slice);
%     g = g(:);
%                                                                        
%     for k = 1:200
%         f_guess = F(:,:,slice);
%         f_guess = f_guess(:);
%         g_bar = A*f_guess;
%         F(:,:,slice) = reshape((f_guess./sum(A)').* (A' * (g./g_bar)),[87,87]);
%         r_error(k) = sum((F_iter(:)-F(:)).^2);
%       
%     end
%     r_error(200+1) = sum((F_iter(:)-F(:)).^2);
% end


plot(0:199,r_error(2:201))
xline(19,'-',{'x = 19: Minimum SSD'})
title('SSD between good images and reconstructed noisy')
xlabel('Iterations')
ylabel('SSD')

%% C4 Using the optimal number of iterations to reconstruct the image and also comparing with FBP 

subplot(1,2,1)
[F_opt, ~] = MLEM(A,G_noisy,[87,87,10],19); 
imshow(imcomplement(F_opt(:,:,10)),[])
title('Iterative reconstruction of G noisy using optimal number of iterations')

subplot(1,2,2)
F_FBP= FBP(G_noisy,0:179,[87,87,26]);
imshow(imcomplement(F_FBP(:,:,10)),[]);
title('FBP of sinogram G noicy ,slice 10, image')

%%
function A = compute_forward_matrix(thetas, L, R, C)

    A = zeros(L*length(thetas), R*C);
    for j = 1:R*C
        % Matix with zeros to calculate probabiliti of pixel on a line 
        start_matrix = zeros(R,C); 
        % Inserting a one at each pixel at a time 
        % To see how that pixel is affected 
        start_matrix(j) = 1;
        % Radontransform of the matrix with one 1
        radtransform = radon(start_matrix);
        cut = radtransform(20:106,:);
        A(:, j) = cut(:);
    end


end

function F = NAIVE(A, G, output_size)
    
    % Naive reconstruction 
    F = zeros(output_size);
    for slice = output_size(3)
        one_slice = G(:,:,slice);
        g = one_slice(:);
        f  = A\g; % Solving the equation to get our flattend image 
        F(:,:,slice) = reshape(f,[87,87]);
    end
end

function [F,r_error] = MLEM(A,G,output_size,nbr_iter)

    % Creating an initial guess matrix with ones 
    F = ones(output_size); 
    

    % A zero vector to store SSD
    r_error = zeros(nbr_iter+1, 1);

    % For loop to reconstrcut each slice 
    for slice = 1:size(F,3)
        
        % Number of iterations 
        for k = 1:nbr_iter
            f_guess = F(:,:,slice);
            f_guess = f_guess(:);
            g = G(:,:,slice);
            g = g(:);
            % Calculating gbar 
            g_bar = A*f_guess;
            % Calculating SSD
            SSD = sum((g(:)-g_bar(:)).^2);
            r_error(k) = SSD;
            % Updating our guess 
            F(:,:,slice) = reshape((f_guess./sum(A)').* (A' * (g./g_bar)),[87,87]);
            
        end
        r_error(nbr_iter+1) = sum((g(:)-g_bar(:)).^2);
        
    end
end

function F = FBP(G,theta,outputsize)

% Function to use iradon with a Hamming filter 
F = zeros(outputsize);

for slice = 1:outputsize(3)
    
    F(:,:,slice) = iradon(G(:,:,slice),theta,'linear',"Hamming",outputsize(1));

end
end














