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
