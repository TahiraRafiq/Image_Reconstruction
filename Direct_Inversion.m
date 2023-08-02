clc
clear all
close all

%% phantom image
FOV=256;
ph=phantom('modified shepp-logan',FOV);
 
Ny=FOV;
Nx=FOV;
Nc=8;
M=1; 
mask=zeros(Ny,Nx);
mask(1:M:Ny, :)= 1;

%% coil's senstivity (y-direction)
c_sens=zeros(Ny,Nx,Nc);
for n=1:Nc
    for i=1:Ny
        for j=1:Nx
 
            mean=(FOV/(Nc+1)).*n;
            var=30;
 
            a = 1/(var*(2*pi)^(0.5));
            b = (i-mean)^2;
            d = 2*((var)^2);
            k = (-1)*(b/d);
 
            c_sens(i,j,n)=a*exp(k);
 
        end
    end
end

% figure1 - c_sens
figure(1),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(c_sens(:,:,n),[])
end


%% coil's image
c_img=zeros(Ny,Nx,Nc);
for n=1:Nc
    c_img(:,:,n)=c_sens(:,:,n) .* ph;
end
 
% figure2 - c_img
 
figure(2),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(c_img(:,:,n),[])
end


%% k-space data
for n=1:Nc
  raw(:,:,n)=fftshift(fft2(fftshift(c_img(:,:,n))));
end 
% figrue3 - raw(k-space)
figure(4),
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(raw(:,:,n)),[])
end


%% Zero filling in k space

for n=1:Nc
    reduced_k(:,:,n) = mask.*raw(:,:,n);
end

figure,
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(reduced_k(:,:,n)),[])
end


%% Undersample coil images

for n=1:Nc
  c_undersample(:,:,n)=ifftshift(ifft2(ifftshift(reduced_k(:,:,n))));
end

figure,
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(c_undersample(:,:,n)),[])
end

%% IDFT along X

s_hat = ifftshift(ifft(ifftshift(reduced_k(:,:,:)), [], 2));
permuted_S = permute(s_hat, [1,3,2]);
reshaped_S = reshape(permuted_S, [], Nx);

figure,
for n=1:Nc
    subplot(2,ceil(Nc/2),n)
    imshow(abs(s_hat(:,:,n)),[])
end

figure,
imshow(abs(reshaped_S),[])


%% Final image
   
final_img=zeros(Nx,Ny);

y= -Ny/2:(Ny/2-1);
B=zeros(Ny,Nx,Nc);

for x=1:Nx
     for n = 1:Nc
        for ky=0:Ny-1      
            B(ky+1,:,n)= c_sens(:,x,n)'.* exp(2*pi*-1i*ky*y/Ny);
        end 
        B(:,:,n) = mask.*B(:,:,n);
     end
     
     permuted_B = permute(B, [1,3,2]);
     reshaped_B = reshape(permuted_B, [], Nx);
     
     final_img(:,x) = pinv(reshaped_B) * reshaped_S(:,x);
end

figure, 
imshow(abs(final_img),[])


figure,
imshow(abs(reshaped_B),[])

%% SOS of images
squared_img = power(abs(c_img), 2);
sum_img = sum(squared_img, 3);
rsos = sqrt(sum_img);


%% Error 

error = (abs(rsos)-abs(final_img)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)
