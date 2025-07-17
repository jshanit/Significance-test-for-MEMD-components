
% % Simulations of signal identification via Noise-aided MEMD. 
% %     Identify the informatic components in the synthetic data [x1,x2,x3] and identify the rossler series [x, y]
% % Steps: cc
% %     1.a noise is added as the additonal channel to synthetic data
% %     2.Noise-aided MEMD is applied to decomposed the composite data
% %     3.Significance test based on Bonferroni correction is applied to identify the informatic IMFs

clear
clc;
% %% Import Data
t=0.001:0.001:1;
a1=sin(2*pi*12*t);
a2=sin(2*pi*26*t);
a3=sin(2*pi*50*t);
c=0.5;%the standard deviation of the input data 
x11=a2+a3+c*randn(size(a1))+a1;
x22=a3+c*randn(size(a1))+a1;
x33=a2+c*randn(size(a1))+a3;
x1=[x11;x22;x33];
x1=x1';
%%%%%%%%%%%%%%%%%%generate rossler series%%%%%%%%
% a=0.2;b=0.2;c=3.5;		
% h=0.075;N=15499;
% [data1,data2,data3] = Rossler(a,b,c,h,N);
% Xr1 = normrnd(0,0.5,[N+1,1]);
% Xr2 = normrnd(0,0.5,[N+1,1]);
% XX2 = data1+Xr1;YY2 = data2+Xr2;
% X=XX2(14501:15500);Y=YY2(14501:15500);
% x1 = [X Y];
%%%%%%%%%%%%import real data for case study%%%%%%%%%%%
% % load SN_SF_CR_1958_2009_month_data

%%%%%%%%%%%%%%%calculate test statistics%%%%%%%%%% 

[N,M]=size(x1);
imf_all=[];
cn=0.000025; % amount of noise added
for pp=1:5
%% Add noise channels
pp
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
        En_w(ii,:)=(w/(cn.^2))*(std(x1(:,ii)).^2);%%%%%w/白噪声的方差再*SN或SF或CR的方差做调整
    end
    En_w_all(pp).En_w=En_w;
    En_all(pp).En_IMF=En_IMF;
    imf_all(pp).imfx=imfx;
end

%%%%%%%construct null distribution%%%%%%%%%%%%%
for aa=1:length(En_all)
    bb1=size(En_all(aa).En_IMF,2);
    bb2=size(En_w_all(aa).En_w,2);
    cc1(aa)=bb1;cc2(aa)=bb2;
end
  L_use1=find(cc1==10);L_use2=find(cc2==10);%calculate power
%    L_use1=find(cc1==11);L_use2=find(cc2==11);%calculate type1
for q=1:length(L_use1)
    E_IMF=En_all(L_use1(q)).En_IMF;
    E_IMF_all(q).E_IMF=E_IMF;
end

for q=1:length(L_use2)
    Ew_IMF=En_w_all(L_use2(q)).En_w;
    Ew_IMF_all(q).Ew_IMF=Ew_IMF;
end

for cf=1:length(Ew_IMF_all)
    p1=Ew_IMF_all(cf).Ew_IMF';
%     p1=p1(1:10,:);
    p2=E_IMF_all(cf).E_IMF';
%     p2=p2(1:10,:);
    pn(:,cf)=p1(:);%%%%噪声的能量
    px(:,cf)=p2(:);%%%三个变量的能量
end
% % %%% Identification by weighted multiple one-tailed test%%%%%%%%%%%%%%%%%

% mu=mean(pn,2);sigma=std(pn',1);
% for jj=1:length(pn(1,:))
%     pn_norm(:,jj)=(pn(:,jj)-mu)./sigma';  
% end
% for ii=1:length(px(1,:))
% px_norm(:,ii)=(px(:,ii)-mu)./sigma';
% end
% neta=max(pn_norm,[],1);
% n_x=max(px_norm,[],1);
% alp=0.05; %% significance level, with alpha = 0.01,0.05,0.1,0.15,0.2,0.25
% En=sort(neta);
% q_cv=En(fix(length(En)*(1-alp))); %% threshold for test
% th_upper=mu+q_cv*sigma';
% px1=mean(px,2);
% r=th_upper<px1; % IMFs are identified in 'r' variable (1-significant; 0 - otherwise)

%%% Identification by Bonferroni correction%%%%%%%%%%%%%%%%%
[M1,N1]=size(px);
alp=0.05;
En=sort(pn,2);
q_cv=En(:,fix(size(En,2)*(1-alp/M1)));%% threshold for test
px1=mean(px,2);
r=px1>q_cv;%%%IMFs are identified in 'r' variable (1-significant; 0 - otherwise)

