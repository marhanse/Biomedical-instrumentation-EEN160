function save_nifti_image(complex_im, filename, meta_data)
    
    % Compute the magnitude image
    magnitude_image = abs(complex_im);
    
    % saves it as a image using the Nifti file format
    niftiwrite(int16(magnitude_image),filename,meta_data)

end

