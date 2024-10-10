%% A2: one backprojection at theta = 30 
load CTdata.mat

sinogram = g;

% Row and Column of the sinogram 
[l, theta] = size(sinogram); 

% Size of the reconstructed image 
image_size = 2*floor(size(sinogram,1)/(2*sqrt(2))); 

% Sum_matrix used to create the image 
sum_matrix = zeros(image_size);

% Creating a matrix with the desiered x and y coordinates
[X,Y] = meshgrid((-(image_size)/2 + 1):(image_size)/2);

% For theta equals 30 degrees 
a = 30;

% Calculating the lines for each value of thetas 
line = X*cosd(a) + Y*sind(a) + floor(l/2);

% Interpolating value of the line to the image.
x = 1:l; % sample points
v = sinogram(:,a);  % corresponding values, v(x)
sum_matrix = sum_matrix + interp1(x, v, line);

% Showing the image
imshow(sum_matrix,[])
title('One backprojection at theta equals 30')


%% A3 Sum of backprojections: 6,12,30,60,180
load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Number of projections
N = [6,12,30,60,180];

% Plotting the result 
p=0;
for nproj = N
    p = p + 1;
    subplot(2,3,p)
    % Implemented backprojection function (step size 1) 
    BP(sinogram,nproj,step)
    title("N = " + nproj)
end 
%% B2 Filtered backprojection FBP
load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Row and Column of the sinogram 
[l, theta] = size(sinogram); 

% Create the normalized ramlak filter
% One for even l sizes and one for odd 
if mod(l,2) == 0
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:1];
else
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:0];
end

% Multiplying the ramlak with fft(g) and taking the real part of the ifft
filtered_sinogram = fft(sinogram).* ramlak';
g_filtered = real(ifft(filtered_sinogram));

% Number of projections 6, 12, 30, 60, 180 
N = [6,12,30,60,180];

% Plotting the result 
p = 0; 
for nproj = N
    p = p + 1;
    subplot(2,3,p)
    BP(g_filtered,nproj,step)
    title("FBP N = " + nproj)
end 

%% B3 FBP with hamming-window 
load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Row and Column of the sinogram 
[l, theta] = size(sinogram); 

% Create the normalized ramlak filter
% One for even l sizes and one for odd 
if mod(l,2) == 0
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:1];
else
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:0];
end


% Creating the hamming window 
omega = 0:1:l-1; c= 0.54;
hamming_window = c + (c-1) .* cos(2.*pi*omega./l);

% Creating the filter 
filter = fftshift(hamming_window) .* ramlak;

% Filtering the sinogram in the fourier domain 
filtered_sinogram = fft(sinogram).*filter'; 

% The filtered sinogram
g_filtered = real(ifft(filtered_sinogram));

% Number of projections
N = [6,12,30,60,180];

% Plotting the result 
p = 0; 
for nproj = N
    p = p + 1;
    subplot(2,3,p)
    BP(g_filtered,nproj,step)
    title("CBP (Hamming) N = " + nproj)
end 

%% C3 Convolution backprojection

load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2; step= 1;
else
    sinogram = g3; step = 0.5; 
end

% Row and Colum n of the sinogram 
[l, theta] = size(sinogram); 


% Create the normalized ramlak filter
% One for even l sizes and one for odd 
if mod(l,2) == 0
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:1];
else
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:0];
end

% Creating the hamming window 
omega = 0:1:l-1; c= 0.54;
hamming_window = c + (c-1) .* cos((2*pi*omega)/l);

% Filter in fourier space 
filter1 = fftshift(hamming_window) .* ramlak;

% inverse FT on the filter 
filter = real(ifftshift((ifft(filter1))));

% Using convolution to filter the sinogram 
filtered_sinogram = conv2(sinogram,filter','same'); 

% Number of projections
N = [6,12,30,60,180];

% Plotting the result 
p = 0; 
for nproj = N
    p = p + 1;
    subplot(2,3,p)
    BP(filtered_sinogram,nproj,step)
    title("CBP N = " + nproj)
end 

%% E1 g2 

%% Normal backprojection with g2,g3 run A3 and with input 2,3 

load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Number of projections
% For g2 we used N = 180
% For g3 we used N = 360
disp('Enter number of projections wanted to reconstruct');
N=input('Sinogram will be reconstructed with N projections');

BP(sinogram,N,step)
title("Backprojection with N = " + N)

%% FBP (hamming window) with g2,g3, run B3 and with input 2,3

load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Row and Column of the sinogram 
[l, theta] = size(sinogram); 

% Create the normalized ramlak filter
if mod(l,2) == 0
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:1];
else
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:0];
end

% Creating the hamming window 
omega = 0:1:l-1; c= 0.54;
hamming_window = c + (c-1) .* cos(2.*pi*omega./l);

% Creating the filter 
filter = fftshift(hamming_window) .* ramlak;

% Filtering the sinogram in the fourier domain 
filtered_sinogram = fft(sinogram).*filter'; 

% The filtered sinogram
g_filtered = real(ifft(filtered_sinogram));

% Number of projections (180 for g2, 360 for g3)
disp('Enter number of projections wanted to reconstruct');
N=input('Sinogram will be reconstructed with N projections');

% Plotting the result 
BP(g_filtered,N,step)
title("FBP (Hamming) N = " + N)

%% Convolution BP with g2, run C3 and with input 2
load CTdata.mat

disp('Enter (1) for g (2) for g2 (3) for g3');
s=input('This sinogram will be reconstructed ');

if s == 1
    sinogram = g; step = 1;
elseif s == 2 
    sinogram = g2;step=1;
else
    sinogram = g3;step = 0.5; 
end

% Row and Colum n of the sinogram 
[l, theta] = size(sinogram); 


% Create the normalized ramlak filter
if mod(l,2) == 0
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:1];
else
    ramlak = 1/(l/2)*[0:1:floor(l/2)-1,floor(l/2):-1:0];
end

% Creating the hamming window 
omega = 0:1:l-1; c= 0.54;
hamming_window = c + (c-1) .* cos((2*pi*omega)/l);

% Filter in fourier space 
filter1 = fftshift(hamming_window) .* ramlak;

% inverse FT on the filter 
filter = real(ifftshift((ifft(filter1))));

% Using convolution to filter the sinogram 
filtered_sinogram = conv2(sinogram,filter','same'); 

% Number of projections (180 for g2, 360 for g3)
disp('Enter number of projections wanted to reconstruct');
N=input('Sinogram will be reconstructed with N projections');

% Plotting the result 
BP(filtered_sinogram,N,step)
title("Convolution BP N = " + N)









