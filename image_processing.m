% This code is contributed by Miss Chatgpt and Mohitha......

% The spatial frequency directly relates with the brightness of the image. The magnitude of the sinusoid directly relates with the contrast. 
% f(x,y) denotes the image , and F(u,v) denotes the discrete Fourier transform.
% Read input image
inputImage = imread('sp.jpg');
inputImage = rgb2gray(inputImage); % Convert to grayscale if necessary

% Perform DFT on the image
[M, N] = size(inputImage);
dftImage = zeros(M, N); % empty matrix same as size of image

for u = 0:M-1
    for v = 0:N-1
        for x = 0:M-1
            for y = 0:N-1
                dftImage(u+1, v+1) = dftImage(u+1, v+1) + double(inputImage(x+1, y+1)) * exp(-1i*2*pi*(u*x/M + v*y/N)); 
                %This represents an element of the DFT image matrix at position (u+1, v+1)
            end
        end
    end
end

% Shift the zero-frequency component to the center of the spectrum
% This shift allows us to observe the low-frequency components in the center of the spectrum, while the high-frequency components are spread towards the edges.
shiftedDftImage = zeros(M, N); %This line initializes a new matrix shiftedDftImage with the same size as the DFT image

for u = 0:M-1
    for v = 0:N-1
        shiftedDftImage(u+1, v+1) = dftImage(u+1, v+1) * (-1)^(u+v);
        %The value of the element in the original DFT image is multiplied by -1 raised to the power of u+v. 
        % This operation applies a phase shift to the corresponding frequency component.
        %When u+v is even, (-1)^(u+v) evaluates to 1, and the value is unchanged.
      %When u+v is odd, (-1)^(u+v) evaluates to -1, effectively applying a phase shift of 180 degrees (or pi radians).
    end
end

% Apply some image processing operation in the frequency domain
% For example, let's enhance high-frequency components by multiplying by a high-pass filter
cutoffFrequency = 0.1 * min(M, N); % Adjust the cutoff frequency as desired
% The cutoff frequency is set as 10% of the minimum dimension of the shifted DFT image
highPassFilter = ones(M, N);
highPassFilter(M/2-cutoffFrequency:M/2+cutoffFrequency, N/2-cutoffFrequency:N/2+cutoffFrequency) = 0;
processedDftImage = shiftedDftImage .* highPassFilter;

% Shift the spectrum back to the original position
shiftedBackDftImage = zeros(M, N);

for u = 0:M-1
    for v = 0:N-1
        shiftedBackDftImage(u+1, v+1) = processedDftImage(u+1, v+1) * (-1)^(u+v);
    end
end

% Perform inverse DFT to obtain the processed image
processedImage = zeros(M, N);

for x = 0:M-1
    for y = 0:N-1
        for u = 0:M-1
            for v = 0:N-1
                processedImage(x+1, y+1) = processedImage(x+1, y+1) + shiftedBackDftImage(u+1, v+1) * exp(1i*2*pi*(u*x/M + v*y/N));
            end
        end
        processedImage(x+1, y+1) = processedImage(x+1, y+1) / (M*N);
    end
end

% Display the original and processed images
subplot(1, 2, 1);
imshow(inputImage);
title('Original Image');
subplot(1, 2, 2);
imshow(processedImage, []);
title('Processed Image');
