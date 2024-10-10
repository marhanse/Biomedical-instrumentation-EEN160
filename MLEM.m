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