clear all
close all

%% image de départ
N=24;
%N=40;
Im=zeros(N,N);
Im(8:16,8:16)=1;
%Im(10:30,10:30)=1;
%Im=phantom('Modified Shepp-Logan',N);
F=Im(:);

%% calcul de la TR

%radon kernel
R=ker_radon(N);
g_ker=R*F;
g_ker_reshap=reshape(g_ker,ceil(N*sqrt(2))+3,180);

%radon matlab function
theta=(0:179);
[g,xp]=radon(Im,theta);
g_mat=g(:);

erreur_gkerVSgmat=immse(g,g_ker_reshap); %error between matlab and 

%%%%% Reconstruction %%%%%

%iradon() matlab function
x_rad=iradon(g,theta); %Reconstructed with g-rad iradon() matlab
x_ker=iradon(g_ker_reshap,theta); %Reconstructed with g-ker iradon() matlab

erreur_irad_mat=immse(Im,x_rad(2:end-1,2:end-1));%with skip 1 pixel on each dimension
erreur_irad_ker=immse(Im,x_ker(2:end-1,2:end-1));%iradon matlab reconstruction with ker RT for construction

%%Pour la suite on prend g_ker_reshap/g_ker pour faire la reconstruction 

%% Probleme en R_carre solution des MC/LS
g_inv=(R')*g_mat;
R_carre=(R')*R;

%% Inv généralisée
%x_gene0=inv(R_carre)*g_inv; inv matlab compliqué warning conditionement 
%x_rgene0=reshape(x_gene0,N,N);
%erreur_invgenezero=immse(Im,x_rgene0);

x_gene=pinv(R)*g_mat;
x_rgene=reshape(x_gene,N,N);
erreur_invgene=immse(Im,x_rgene);

%SVD 
[U,S,V]=svd(R);
sigma=svd(R);
Spinv=pinv(S);

%SVD not troncated equivalent to Generalized Inverse
F_svd=V*Spinv*(U')*g_mat;
x_svd=reshape(F_svd,N,N);
erreur_svd=immse(Im,x_svd);

%% Troncated SVD <=> ACP ?!
K=330; %K=575;  
F_re3=V(:,1:K)*Spinv(1:K,:)*(U')*g_mat;
x_tsvd=reshape(F_re3,N,N);
erreur_tsvd=immse(Im,x_tsvd);

%% ART : long 2-3 mins
%%f_art=ART(R,g_mat,N);
%%x_art=reshape(f_art(:,end),N,N);
%%erreur_art=immse(Im,x_art);

%% descente de grad

x_grad=grad(R_carre,g_inv,N);
x_rgrad=reshape(x_grad(:,end),N,N);

erreur_grad=immse(Im,x_rgrad);

%% grad conjugué
cond_Rcarre=cond(R_carre);

x_gradC=gradC(R_carre,g_inv,N);
x_rgradC=reshape(x_gradC(:,end),N,N);
erreur_gradC=immse(Im,x_rgradC);

%% regularisation de tikhonov
alpha=10;
R_tik=R_carre+alpha*eye(N*N);

cond_Rtik=cond(R_tik);
x_tik=grad(R_tik,g_inv,N);
x_rtik=reshape(x_tik(:,end),N,N);

erreur_tik=immse(Im,x_rtik);

%% Plot

figure(1)
subplot(211)
imshow(g,[])
colormap(hot)
colorbar
title("Radon Transform from matlab")
subplot(212)
imshow(g_ker_reshap,[])
colormap(hot)
colorbar()
title("Radon Transform from discreet construction with matrix R")

figure(2)
subplot(221)
imshow(x_rgene,[])
colorbar
title("Generalized inverse")
subplot(222)
imshow(x_svd,[])
colorbar
title("SVD")
subplot(223)
imshow(x_tsvd,[])
colorbar
title("TSVD")
subplot(224)
%imshow(x_art,[])
%colorbar
title("ART reconstruction")

figure(3)
subplot(221)
imshow(Im,[])
colorbar
title("Image ref")
subplot(222)
imshow(x_rgrad,[])
colorbar
title("Gradient Descente")
subplot(223)
imshow(x_rgradC,[])
colorbar
title("Conjugated Gradient")
subplot(224)
imshow(x_rtik,[])
colorbar
title("Conjugated Gradient + Tik")

% figure(4)
% imshow(R,[])% imshow(Im,[])
% title("Matrix R")

% figure(4)
% subplot(131)
% imshow(Im,[])
% title("Image Ref")
% subplot(132)
% imshow(x_rad,[])
% title("g-matlab mes")
% subplot(133)
% imshow(x_ker,[])
% title("g-ker mes")

% figure(4)
% plot(sigma(1:K))
% title("Singular values of R")

% figure(5)
% subplot(121)
% imshow(x_svd,[])
% title("SVD")
% subplot(122)
% imshow(x_tsvd,[])
% title("TSVD")

% figure(5)
% subplot(121)
% imshow(x_svd,[])
% title("Gradient Descent")
% subplot(122)
% imshow(x_tsvd,[])
% title("TSVD")









