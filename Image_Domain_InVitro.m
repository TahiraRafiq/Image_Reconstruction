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
% 
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
%% 2D IDFT

S_along_x=zeros(Nx,Ny,Nc);
for n=1:Nc
S_along_x(:,:,n) = (ifftshift(ifft(ifftshift(k_space_undersampling(:,:,n),2),[],2),2)); 
end

S_along_y=zeros(Nx,Ny,Nc);
for n=1:Nc
S_along_y(:,:,n) = (ifftshift(ifft(ifftshift(S_along_x(:,:,n),1),[],1),1)); 
end

figure
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(S_along_y(:,:,n)),[])
end


delta=Ny/M;
recon_img=zeros(Nx,Ny);

for x=1:Nx
     for y=1:delta
          for L=1:Nc
              B(L,1:M)=c_sens(y:delta:end,x,L);
              pixel_vector(L,1)=S_along_y(y,x,L);
          end
          invB=pinv(B);
          recon_img(y:delta:end,x)=invB*pixel_vector;
     end
end

figure,
imshow(abs(recon_img),[])


figure,
imshow(abs(abs(ph)-abs(recon_img)),[])

%% SOS of images
squared_img = power(abs(c_img), 2);
sum_img = sum(squared_img, 3);
rsos = sqrt(sum_img);

%% Error 

error = (abs(ph)-abs(recon_img)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)
