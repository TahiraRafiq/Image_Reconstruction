clc
clear all
close all

%% Load file
load('Brain2D');

%% Parameters
FOV=256;
Nc = 12;
Nx =  FOV;
Ny =  FOV;
rate = 8;
M=rate;
mask=zeros(Nx,Ny);
mask(1:M:end,:)=1;
block_size=256;
Total_lines=(block_size/2)*Nc;
delta_ky=(2*pi)/256;
num_harmonics=block_size;
E=block_size/M;
EM=E*M;

%% Coil images
coil_img=ifftshift(ifft2(ifftshift(DATA)));


%% Combined SOS image
 
squared_img = power(abs(coil_img), 2);
sum_img = sum(squared_img, 3);
rsos = sqrt(sum_img);
 figure,
 imshow(abs(rsos),[])

%% Coil sensitivity
[Ny,Nx,Nc]= size(DATA);
rawDataKs = DATA; 
%rawDataKs = normalizeData(DATA);
imshowCoilData(1,rawDataKs,'kx','ky','rawDataKs:');
% sum of quare image to calculate optimizing image
rawDataIm=zeros(Ny,Nx,Nc);
for n = 1:Nc
rawDataIm(:,:,n) = fftshift(ifft2(ifftshift(rawDataKs(:,:,n))));
end
% imshowCoilData(2,rawDataIm,'x','y','rawDataIm:');
% sum of square of coil Imaging 
rmsOptImg  = sqrt(sum(abs(rawDataIm).^2,3)); 
figure(3),
imshow(abs(rmsOptImg),[]),title(['Image(Reference)'],'FontSize',12);
% get low spatial resolution for  estimation sens 
nK =20;
litLineKy = (129-nK):1:(129+nK);
litLineKx = (129-nK):1:(129+nK);
rawDataKs_LowRes=zeros(Ny,Nx,Nc);
for n = 1:Nc
rawDataKs_LowRes(litLineKy,litLineKx,n) = rawDataKs(litLineKy,litLineKx,n);
end
imshowCoilData(4,rawDataKs_LowRes,'kx','ky','rawDataKs-LitK:');
rawDataEstIm=zeros(Ny,Nx,Nc);
for n = 1:Nc
rawDataEstIm(:,:,n) = fftshift(ifft2(ifftshift(rawDataKs_LowRes(:,:,n))));
end
estCoilSenSmth = estimateCoilSens_FilterKs(rawDataKs_LowRes,nK);
imshowCoilData(5,estCoilSenSmth,'x','y','estCoilSenSmth:');
coilSens = estCoilSenSmth;

%% undersampling k-space

k_space_undersampling=zeros(Nx,Ny,Nc);

for n=1:Nc
   k_space_undersampling(:,:,n) = mask.*DATA(:,:,n); 
end

S_along_x=zeros(Ny,Nx,Nc);
for n=1:Nc
S_along_x(:,:,n) = (ifftshift(ifft(ifftshift(k_space_undersampling(:,:,n),2),[],2),2)); 
end
% permute individual k-space 
% permuted_raw=permute(S_along_x,[3,1,2]);
% reshaped_raw=reshape(permuted_raw,[],Ny);

%% Target harmonics
y=[-Ny/2:Ny/2-1];
ideal_harmonics=zeros(Ny,block_size);
for n=0:block_size-1
ideal_harmonics(:,n+1)=exp(1i*n*delta_ky*y);
end


%% 
% calculate Weight
nWeightFactor =  zeros(EM,E*Nc,Nx); % m,E*n,x  
permuted_coilSens = permute(coilSens(:,:,:),[3 1 2]);% y,x,n -> n,y,x

%calculate BsubMat matrix co
B = zeros(E*Nc,Ny,Nx); %(E*n,y,x) 

for x = 1:Nx
    for e=0:(E-1)
        for y =-(Ny/2):1:(Ny/2-1)
            B(1+(e*Nc):Nc+(e*Nc),y+(Ny/2)+1,x) = permuted_coilSens(:,y+(Ny/2)+1,x)*exp(-2*pi*i*(e*M)*y/Ny);
        end
    end
end

for x = 1:Nx
    smaller_B = B(:,:,x); % E*n,y,x 
    nWeightFactor(:,:,x) = ideal_harmonics(:,:)'*pinv(smaller_B); % (m,E*n) = (m,y)*(y,E*n) 
end

% compute the missing point ky 
permuted_k_space = permute(S_along_x(:,:,:),[3 1 2]); % (ky,x,n)->(n,ky,x) 

%Calculate matrix S_combine
Composite_k_space = zeros(E*Nc,Ny,x); 
ky=1:EM:Ny;
for e = 0:(E-1)
    Composite_k_space(1+(e*Nc):Nc+(e*Nc),ky,:) = permuted_k_space(:,ky+e*M,:);
end

S_Comp = zeros(Ny,Nx);
threshold=0.05;
kySampLoc = 1:EM:Ny; 

for x= 1:Nx
    % calulcate nweightFactor
    C1 = B(:,:,x); % E*n,y,x 
    
    %numerical coditioning 
    [Cl_U,Cl_S,Cl_V] = svd(C1);
    Cl_S_MaxValue = max(Cl_S,[],[1 2]);  
    Cl_S(Cl_S<0.05*Cl_S_MaxValue)=0;
    Cl_S(Cl_S<threshold*Cl_S_MaxValue)=0;
    C2 = Cl_U*Cl_S*ctranspose(Cl_V);
    %
    nWeightFactor(:,:,x) = ideal_harmonics(:,:)'*pinv(C2); % (m,E*n) = (m,y)*(y,E*n) 
    
 
    for mShift=0:(EM-1)
        synthdata=squeeze(nWeightFactor(mShift+1,:,x))*(Composite_k_space(:,kySampLoc,x)); % (m,ky_:)= (m,E*n_:)*(E*n_:,ky_:)
        S_Comp(kySampLoc+mShift,x) = synthdata;
    end
end



recon_image=ifftshift(ifft(ifftshift(S_Comp,1),[],1),1); 
figure,
imshow(abs(recon_image),[])

error = (abs(rsos)-abs(recon_image)).^2;
RMSE = sqrt(sum(error(:))/(Nx * Ny));
NRMSE = RMSE/(Nx*Ny)

diff = (abs(rsos)-abs(recon_image));
figure,
imshow(abs(diff),[])