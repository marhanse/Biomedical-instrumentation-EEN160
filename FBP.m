function F = FBP(G,theta,outputsize)

% Function to use iradon with a Hamming filter 
F = zeros(outputsize);

for slice = 1:outputsize(3)
    
    F(:,:,slice) = iradon(G(:,:,slice),theta,'linear',"Hamming",outputsize(1));

end

