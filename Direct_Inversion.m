clc
clear all
close all

load('modifiedshep.mat');
FOV=256;
ph=phantom('modified shepp-logan',FOV);
Nc = 8;
Nx =  FOV;
Ny =  FOV;
rate = 8;
figure(1) ;
imshow(ph,[])

%% Images of each coils

for n=1:Nc
    c_img1(:,:,n) = ph.*c_sens(:,:,n); 
end

c_raw=fftshift(fft2(fftshift(c_img1)));

mask = sort([Ny/2+rate:rate:Ny, Ny/2:-rate:1 ]);
for n=1:Nc
k_space_undersampling(:,:,n) = c_raw(mask,:,n); 
end
figure(5);
for n=1:Nc
    subplot(2,ceil(Nc/2),n);
    imshow(abs(k_space_undersampling(:,:,n)),[]);
    title('Under-Sampled kspace');
end


%% 2D inverse Discrete Fourier Transform in Seperated directions
k_space = fftshift(k_space_undersampling);
inner = zeros(Ny/rate,Nx,Nc);
inner_sum = 0;
% Inverse Discrete Fourier Transform along the x-direction
for n=1:Nc
    for ky=0:Ny/rate-1
        for x=0:Nx-1 % index in spatial domain
            for kx=-Nx/2:Nx/2-1 % index in K space
                 inner_sum = inner_sum + k_space(ky+1,kx+Nx/2+1,n)*exp(1i*2*pi*(kx/Nx)*x); 
            end 
            inner(ky+1,x+1,n) = inner_sum;
            inner_sum = 0;
        end  
    end
end
S_iDFT_x_direction = ifftshift(inner/Nx);

permuted_S = permute(S_iDFT_x_direction, [1,3,2]);
S_iDFT_x_direction_group = reshape(permuted_S, [], Nx);

%%  Reconstructing image
BB_hat = zeros(Ny,Ny,Nc);
ReImage = zeros(Nx,Ny); 
for x=1:Nx
    B_hat_undersampling = [];
    for n=1:Nc
        for  ky=-Ny/2:Ny/2-1
            for  y=-Ny/2:Ny/2-1
                BB_hat(ky+1+Ny/2,y+1+Ny/2,n) = c_sens(y+1+Ny/2,x,n).*exp(-1i*2*pi*ky*y/Ny);
                
            end  
        end
    end
    for n=1:Nc
    temp = BB_hat(1+(rate-1):rate:end,:,n);
    B_hat_undersampling = [B_hat_undersampling;temp];
    end 
     total_B(:,:,x)=B_hat_undersampling;
end
for x=1:256
ReImage(:,x) = pinv(total_B(:,:,x))*S_iDFT_x_direction_group(:,x); 
end
figure,
imshow(abs(ReImage),[])


%% SOS of images
squared_img = power(abs(c_img1), 2);
sum_img = sum(squared_img, 3);
rsos = sqrt(sum_img);


%% Error 

error = (abs(ph)-abs(sqrt(ReImage)).^2);
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)



