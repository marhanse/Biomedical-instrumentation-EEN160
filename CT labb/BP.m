function BP(sinogram, nproj,step)
    
    % Row and Column of the sinogram 
    [l, theta] = size(sinogram); 
    
    % Size of the reconstructed image 
    image_size = 2*floor(size(sinogram,1)/(2*sqrt(2))); 
    
    % Sum_matrix used to create the image 
    sum_matrix = zeros(image_size);
    
    % Creating a matrix with the desiered x and y coordinates
    [X,Y] = meshgrid((-(image_size)/2 + 1):(image_size)/2);
    
        for a = 1:floor(theta/nproj):theta
    
           % Calculating the lines for each value of thetas depending on
           % step size! 

           if step == 0.5
                line = X*cosd(a/2) + Y*sind(a/2) + floor(l/2);
           else
               line = X*cosd(a) + Y*sind(a) + floor(l/2);
           end

           % Interpolating value of the line to the image.
           x = 1:l; % sample points
           v = sinogram(:,a); % corresponding values, v(x)
           % Adding up for each theta to get the BP image 
           sum_matrix = sum_matrix + interp1(x, v, line);
    
        end
    
        % Multiplying factor 
        f = pi/(2*nproj);
        sum_matrix = sum_matrix * f;
    
        imshow(sum_matrix,[])
end

