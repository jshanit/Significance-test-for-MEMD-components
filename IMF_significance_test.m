
% % Simulations of signal identification via Noise-aided MEMD. 
% %     Identify the informatic components in the synthetic data [x1,x2,x3]
% % Steps: cc
% %     1.Several noises are added as the additonal channels to synthetic data
% %     2.Noise-aided MEMD is applied to decomposed the composite data
% %     3.Significance test of EEMD is applied to identify the informatic IMFs

clear
clc;
% % load SN_SF_CR_1958_2009_month_data
% %load rossler_1000
% %% Import Data
% % t=0.001:0.001:1;
% % a1=sin(2*pi*12*t);
% % a2=sin(2*pi*26*t);
% % a3=sin(2*pi*50*t);
% % c=0.5;%the standard deviation of the input data 
% % x11=a2+a3+c*randn(size(a1))+a1;
% % x22=a3+c*randn(size(a1))+a1;
% % x33=a2+c*randn(size(a1))+a3;
% % x1=[x11;x22;x33];
% % % x1=[x11/std(x11);x22/std(x22);x33/std(x33)];
% % x1=x1';
%%%%%%%%%%%%%%%%%%generate rossler series
% a=0.2;b=0.2;c=3.5;		
% h=0.075;N=15499;
% [data1,data2,data3] = Rossler(a,b,c,h,N);
% Xr1 = normrnd(0,0.5,[N+1,1]);
% Xr2 = normrnd(0,0.5,[N+1,1]);
% % XX1 = data1+Xr;YY1 = data2+Xr;
% % X = XX1(15001:15500);Y = YY1(15001:15500);
% XX2 = data1+Xr1;YY2 = data2+Xr2;
% X=XX2(14001:15500);Y=YY2(14001:15500);
% x1 = [X Y];
% % load SN_SF_CR_1958_2009_month_data
% [N,M]=size(x1);
% imf_all=[];
% % for i=1:M
% % xx1(:,i)=(x1(:,i)-mean(x1(:,i)))./std(x1(:,i)); 
% % end

% cn=0.000025; % amount of noise added,范围在2%-10%之间
% for pp=1:5000
% %% Add noise channels
% pp
% x=[x1 cn*randn(N,1)];
% 
% %% MEMD in work
% imfx=memd_fast(x,150);
% % imfx=memd(x,64,'stop',[0.05 0.5 0.05]);
% % imfx=memd(x,64,'fix_h',10);
%     IMF_all=imfx(1:M,:,:);
%     mm=size(IMF_all,1);
%     En_IMF=[];En_w=[];w=[];
%     for i=1:mm
%         En_IMF(i,:)=mean(IMF_all(i,:,:).^2,3);
%     end
%     IMF_noise_all=imfx(end,:,:);IMF_noise_all=permute(IMF_noise_all,[2 3 1]);
%     m1=size(IMF_noise_all,1);
%     for j=1:m1
%         w(j,:)=mean(IMF_noise_all(j,:).^2,2);
%     end
%     for ii=1:mm
%         En_w(ii,:)=(w/(cn.^2))*(std(x1(:,ii)).^2);%%%%%w/白噪声的方差再*SN或SF或CR的方差做调整
%     end
%     En_w_all(pp).En_w=En_w;
%     En_all(pp).En_IMF=En_IMF;
%     imf_all(pp).imfx=imfx;
% end

% % % load SN_SF_CR_1958_2009_month_data
% load simu_series_imf_cf100
load rossler_test_dec_result_M_5000

for aa=1:length(En_all)
    bb1=size(En_all(aa).En_IMF,2);
    bb2=size(En_w_all(aa).En_w,2);
    cc1(aa)=bb1;cc2(aa)=bb2;
end
  L_use1=find(cc1==9);L_use2=find(cc2==9);%calculate power
%    L_use1=find(cc1==11);L_use2=find(cc2==11);%calculate type1
for q=1:length(L_use1)
    E_IMF=En_all(L_use1(q)).En_IMF;
    E_IMF_all(q).E_IMF=E_IMF;
end

for q=1:length(L_use2)
    Ew_IMF=En_w_all(L_use2(q)).En_w;
    Ew_IMF_all(q).Ew_IMF=Ew_IMF;
end

%%%%%%%%construct null distribution%%%%%%%%%%%%%
for cf=1:length(Ew_IMF_all)
    p1=Ew_IMF_all(cf).Ew_IMF';
%     p1=p1(1:10,:);
    p2=E_IMF_all(cf).E_IMF';
%     p2=p2(1:10,:);
    pn(:,cf)=p1(:);%%%%噪声的能量
    px(:,cf)=p2(:);%%%三个变量的能量
end
% % %%% Identification by multiple testing%%%%%%%%%%%%%%%%%

% mu=mean(pn,2);sigma=std(pn',1);
% for jj=1:length(pn(1,:))
%     pn_norm(:,jj)=(pn(:,jj)-mu)./sigma';  
% end
% for ii=1:length(px(1,:))
% px_norm(:,ii)=(px(:,ii)-mu)./sigma';
% end
% neta=max(pn_norm,[],1);
% n_x=max(px_norm,[],1);
% alp=0.05; %% significance level, with alpha = 0.001,0.005,0.01,0.05,0.1
% En=sort(neta);
% q_cv=En(fix(length(En)*(1-alp))); %% threshold for test
% th_upper=mu+q_cv*sigma';
% % px1=mean(px,3);px1=px1';
% % px1=px1(:);
% % r=th_upper<px1; % IMFs are identified in 'r' variable (1-significant; 0 - otherwise)
%%% Identification by single testing with Bonferroni correction%%%%%%%%%%%%%%%%%
[M1,N1]=size(px);
alp=0.05;
En=sort(pn,2);
% th=wn(:,fix(size(wn,2)*(1-alp/2))); %% threshold for test
q_cv=En(:,fix(size(En,2)*(1-alp/M1)));
% q_cv=En(:,fix(size(En,2)*(1-alp/30)));%%calculate power
px1=mean(px,2);
r=px1>q_cv;

% %%%%%%%%%%%%%%%%%%%%%decomposition figure%%%%%%%%%%%%%%%%
% % % 
% allmode=permute(IMF_all,[3 2 1]);%allmode(数据长度：IMF的个数：变量个数)
% [m1,n,nn]=size(allmode);
% for i=1:nn
% figure(i); 
% allmode1 = allmode(:,:,i);nimf = size(allmode1,2); 
% for m = 1:nimf
%   subplot(nimf,1,m); plot(1:m1,allmode1(:,m),'k','LineWIdth',1.5);
% end
% subplot(nimf,1,1); title('IMFs');
% end
%%%%%%%%%%%%%%%%%%%case study figure2%%%%%%%%%%%%%%%%%%
% de=IMF_all;
% de=imfx;
% figure
% subplot(9,4,1);plot(squeeze(de(1,1,:)));ylabel('C^(1)');set(gca,'XTick',[])
% subplot(9,4,2);plot(squeeze(de(2,1,:)));set(gca,'XTick',[]);
% subplot(9,4,3);plot(squeeze(de(3,1,:)));set(gca,'XTick',[]);
% subplot(9,4,4);plot(squeeze(de(4,1,:)));set(gca,'XTick',[]);
% subplot(9,4,5);plot(squeeze(de(1,2,:)));ylabel('c^2');set(gca,'XTick',[])
% subplot(9,4,6);plot(squeeze(de(2,2,:)));set(gca,'XTick',[]);
% subplot(9,4,7);plot(squeeze(de(3,2,:)));set(gca,'XTick',[]);
% subplot(9,4,8);plot(squeeze(de(4,2,:)));set(gca,'XTick',[]);
% subplot(9,4,9);plot(squeeze(de(1,3,:)));ylabel('c^3');set(gca,'XTick',[])
% subplot(9,4,10);plot(squeeze(de(2,3,:)));set(gca,'XTick',[]);
% subplot(9,4,11);plot(squeeze(de(3,3,:)));set(gca,'XTick',[]);
% subplot(9,4,12);plot(squeeze(de(4,3,:)));set(gca,'XTick',[]);
% subplot(9,4,13);plot(squeeze(de(1,4,:)));ylabel('c^4');set(gca,'XTick',[])
% subplot(9,4,14);plot(squeeze(de(2,4,:)));set(gca,'XTick',[]);
% subplot(9,4,15);plot(squeeze(de(3,4,:)));set(gca,'XTick',[]);
% subplot(9,4,16);plot(squeeze(de(4,4,:)));set(gca,'XTick',[]);
% subplot(9,4,17);plot(squeeze(de(1,5,:)));ylabel('c^5');set(gca,'XTick',[])
% subplot(9,4,18);plot(squeeze(de(2,5,:)));set(gca,'XTick',[]);
% subplot(9,4,19);plot(squeeze(de(3,5,:)));set(gca,'XTick',[]);
% subplot(9,4,20);plot(squeeze(de(4,5,:)));set(gca,'XTick',[]);
% subplot(9,4,21);plot(squeeze(de(1,6,:)),'r');ylabel('c^6');set(gca,'XTick',[])
% subplot(9,4,22);plot(squeeze(de(2,6,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,23);plot(squeeze(de(3,6,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,24);plot(squeeze(de(4,6,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,25);plot(squeeze(de(1,7,:)),'r');ylabel('c^7');set(gca,'XTick',[])
% subplot(9,4,26);plot(squeeze(de(2,7,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,27);plot(squeeze(de(3,7,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,28);plot(squeeze(de(4,7,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,29);plot(squeeze(de(1,8,:)),'r');ylabel('c^8');set(gca,'XTick',[])
% subplot(9,4,30);plot(squeeze(de(2,8,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,31);plot(squeeze(de(3,8,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,32);plot(squeeze(de(4,8,:)),'r');set(gca,'XTick',[]);
% subplot(9,4,33);plot(squeeze(de(1,9,:)));ylabel('c^9');
% subplot(9,4,34);plot(squeeze(de(2,9,:)));
% subplot(9,4,35);plot(squeeze(de(3,9,:)));
% subplot(9,4,36);plot(squeeze(de(4,9,:)));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%simulation_test_figure1%%%%%%%%%%%%
% de=[IMF_all;imfx(4,:,:)];
% subplot(10,4,1);plot(squeeze(de(1,1,:)),'k');ylabel('c^1');set(gca,'XTick',[])
% subplot(10,4,2);plot(squeeze(de(2,1,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,3);plot(squeeze(de(3,1,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,4);plot(squeeze(de(4,1,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,5);plot(squeeze(de(1,2,:)),'k');ylabel('c^2');set(gca,'XTick',[])
% subplot(10,4,6);plot(squeeze(de(2,2,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,7);plot(squeeze(de(3,2,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,8);plot(squeeze(de(4,2,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,9);plot(squeeze(de(1,3,:)),'k');ylabel('c^3');set(gca,'XTick',[])
% subplot(10,4,10);plot(squeeze(de(2,3,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,11);plot(squeeze(de(3,3,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,12);plot(squeeze(de(4,3,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,13);plot(squeeze(de(1,4,:)),'r');ylabel('c^4');set(gca,'XTick',[])
% subplot(10,4,14);plot(squeeze(de(2,4,:)),'r');set(gca,'XTick',[]);
% subplot(10,4,15);plot(squeeze(de(3,4,:)),'r');set(gca,'XTick',[]);
% subplot(10,4,16);plot(squeeze(de(4,4,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,17);plot(squeeze(de(1,5,:)),'r');ylabel('c^5');set(gca,'XTick',[])
% subplot(10,4,18);plot(squeeze(de(2,5,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,19);plot(squeeze(de(3,5,:)),'r');set(gca,'XTick',[]);
% subplot(10,4,20);plot(squeeze(de(4,5,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,21);plot(squeeze(de(1,6,:)),'r');ylabel('c^6');set(gca,'XTick',[])
% subplot(10,4,22);plot(squeeze(de(2,6,:)),'r');set(gca,'XTick',[]);
% subplot(10,4,23);plot(squeeze(de(3,6,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,24);plot(squeeze(de(4,6,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,25);plot(squeeze(de(1,7,:)),'k');ylabel('c^7');set(gca,'XTick',[])
% subplot(10,4,26);plot(squeeze(de(2,7,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,27);plot(squeeze(de(3,7,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,28);plot(squeeze(de(4,7,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,29);plot(squeeze(de(1,8,:)),'k');ylabel('c^8');set(gca,'XTick',[])
% subplot(10,4,30);plot(squeeze(de(2,8,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,31);plot(squeeze(de(3,8,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,32);plot(squeeze(de(4,8,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,33);plot(squeeze(de(1,9,:)),'k');ylabel('c^9');set(gca,'XTick',[])
% subplot(10,4,34);plot(squeeze(de(2,9,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,35);plot(squeeze(de(3,9,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,36);plot(squeeze(de(4,9,:)),'k');set(gca,'XTick',[]);
% subplot(10,4,37);plot(squeeze(de(1,10,:)),'k');ylabel('c^{10}');
% subplot(10,4,38);plot(squeeze(de(2,10,:)),'k');
% subplot(10,4,39);plot(squeeze(de(3,10,:)),'k');
% subplot(10,4,40);plot(squeeze(de(4,10,:)),'k');


  


