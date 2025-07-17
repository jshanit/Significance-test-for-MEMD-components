
% % Simulation to test performance of the method in the presence of various
%%% noise

tic
clear

t=0.001:0.001:1;

a1=sin(2*pi*12*t);
a2=sin(2*pi*26*t);
a3=sin(2*pi*50*t);

c=[0.2 0.5 0.8 1 1.2 1.5 1.8 2];%the standard deviation of the input data

for p=1:length(c)
    En_w_all=[];
    En_all=[];
    imf_all=[];
for pp=1:10000
     pp
     imfx=[];
%%%%%%%%%%%%%%%generate data according to H1%%%%%%%%%%
% %     x11=a2+a3+c(p)*randn(size(a1))+a1;
% %     x22=a3+c(p)*randn(size(a1))+a1;
% %     x33=a2+c(p)*randn(size(a1))+a3;

% %%%%%%%%%%%%generate data according to H0%%%%%%%%%%
         x11=c(p)*randn(size(a1));
         x22=c(p)*randn(size(a1));
         x33=c(p)*randn(size(a1));
         x1=[x11;x22;x33];
         x1=x1';
        [N,M]=size(x1);
    
    cn=0.000025; % amount of noise added,
    %% Add noise channels
    x=[x1 cn*randn(N,1)];
 
    %% MEMD in work
    imfx=memd_fast(x,150);
    IMF_all=imfx(1:M,:,:);
    mm=size(IMF_all,1);
    En_IMF=[];En_w=[];w=[];
    for i=1:mm
        En_IMF(i,:)=mean(IMF_all(i,:,:).^2,3);
    end
    IMF_noise_all=imfx(end,:,:);IMF_noise_all=permute(IMF_noise_all,[2 3 1]);
    m1=size(IMF_noise_all,1);
    for j=1:m1
        w(j,:)=mean(IMF_noise_all(j,:).^2,2);
    end
    for ii=1:mm
        En_w(ii,:)=(w/(cn.^2))*(std(x1(:,ii)).^2);%%%%%adjust the noise energy%%%%%%%
    end
    En_w_all(pp).En_w=En_w;
    En_all(pp).En_IMF=En_IMF;
    imf_all(pp).imfx=imfx;
end
 % %  save(['test_power_c',num2str(p)]);%%%%%%%%calculate power%%%%
    save(['test_typeI_error_c',num2str(p)]);%%%%%%%%calculate type I error%%%%
end
toc
