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