function [F,r_error] = MLEM1(A,G1,F_good,output_size,nbr_iter)

    % ONLY A COPY OF MLEM% 
    % ONLY DIFFERENCE I THE WAY OF CUMPUTING SSD%
    % SSD BETWEEN GOOD IMAGE AND RECONSTRUCTION OF NOISY%
    % Creating an initial guess matrix with ones 
    F = ones(output_size);    

    r_error = zeros(nbr_iter+1, 1);
    for slice = 1:size(F,3)
        g1 = G1(:,:,slice);
        g1 = g1(:);
                                                                           
        for k = 1:nbr_iter
            f_guess = F(:,:,slice);
            f_guess = f_guess(:);
            g_bar = A*f_guess;
            F(:,:,slice) = reshape((f_guess./sum(A)').* (A' * (g1./g_bar)),[87,87]);
            r_error(k) = sum((F_good(:)-F(:)).^2);
          
        end
        r_error(nbr_iter+1) = sum((F_good(:)-F(:)).^2);
    end
end

