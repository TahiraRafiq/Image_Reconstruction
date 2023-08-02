clc
clear all
close all

load('modifiedshep.mat');
FOV=256;
ph=phantom('modified shepp-logan',FOV);
Nc = 8;
Nx =  FOV;
Ny =  FOV;
rate = 2;
figure(1) ;
imshow(ph,[])

%% Images of each coils

for n=1:Nc
    c_img1(:,:,n) = ph.*c_sens(:,:,n); 
end

c_raw=fftshift(fft2(fftshift(c_img1)));

%% undersampling k

M=rate;
mask=zeros(Nx,Ny);
mask(1:M:end,:)=1;
k_space_undersampling=zeros(Nx,Ny,Nc);

for n=1:Nc
k_space_undersampling = mask.*c_raw; 
end
figure(5);
for n=1:Nc
    subplot(2,ceil(Nc/2),n);
    imshow(abs(k_space_undersampling(:,:,n)),[]);
%     title('Under Sampled k space');
end




C=zeros(Ny,Nc)
for n=1:Nc
    C(:,n)=c_sens(:,1,n);
end

y=(1:256)';
delta_ky=(2*pi)/256;
Weight=zeros(Nc,M)
for n=1:M
    Weight(:,n)=pinv(C)*exp(1i*-(n-1)*delta_ky*y);
    
end


harmonics=zeros(Ny/M,Nx,M);
for n=1:M
    for c=1:Nc
        harmonics(:,:,n)=harmonics(:,:,n)+(Weight(c,n)*k_space_undersampling(1:M:end,:,c));
    end
end


%% Final Image

full_kspace=zeros(Nx,Ny);
for i=1:M
full_kspace(i:M:end,:)=harmonics(:,:,i);
end
image=ifftshift(ifft2(ifftshift(full_kspace)),2);

figure,
imshow(abs(image),[])
 
%% Difference Image
figure,
imshow(abs(abs(ph)-abs(image)),[])


error = (abs(ph)-abs(image)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)