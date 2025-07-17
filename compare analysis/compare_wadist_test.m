clear
clc
 load rossler_dec_M_5000
% % % %% Import Data
% % % % t=0.001:0.001:1;
% % % % a1=sin(2*pi*12*t);
% % % % a2=sin(2*pi*26*t);
% % % % a3=sin(2*pi*50*t);
% % % % c=0.1;%the standard deviation of the input data 
% % % % x11=a2+a3+c*randn(size(a1))+a1;
% % % % x22=a3+c*randn(size(a1))+a1;
% % % % x33=a2+c*randn(size(a1))+a3;
% % % % x1=[x11;x22;x33];
% % %
% % 
% a=0.2;b=0.2;c=3.5;		%改变C可以生成混沌序列
% h=0.075;N=15499;
% [data1,data2,data3] = Rossler(a,b,c,h,N);
% Xr1 = normrnd(0,0.5,[N+1,1]);
% Xr2 = normrnd(0,0.5,[N+1,1]);
% XX2 = data1+Xr1;YY2 = data2+Xr2;
% X=XX2(14501:15500);Y=YY2(14501:15500);
% x1 = [X Y];

% % % x1=x1';
% [N,M]=size(x1);
% imf_all=[];
% for i=1:M
% xx1(:,i)=(x1(:,i)-mean(x1(:,i)))./std(x1(:,i)); 
% end
% cn=0.15; % amount of noise added,范围在2%-10%之间
% for pp=1:30
% %% Add noise channels
% pp
% x=[xx1 cn*randn(N,1)];
% %% MEMD in work
% imfx=memd_fast(x,150);
% % imfx=memd(x,64,'stop',[0.05 0.5 0.05]);
% % imfx=memd(x,64,'fix_h',10);
% imf_all(pp).imfx=imfx;
% end
cn=0.000025;
for cf=1:100
    IMF(:,:,:,cf)=imf_all(cf).imfx;
end
[ch_n, nimf_n, tim_n, cf_n]=size(IMF);
IMF_ori=mean(IMF(1:ch_n-1,:,:,:),4);
IMF_noise=IMF(end,:,:,:);IMF_noise_all=permute(IMF_noise,[4 2 3 1]);
[ch1, nimf1, tim1]=size(IMF_noise_all);
for i=1:ch1
    for i1=1:nimf1
    IMF_noise_norm(i,i1,:)=IMF_noise_all(i,i1,:)/cn;
    end
end
IMF_all=[IMF_ori;IMF_noise_norm];
[ch, nimf, tim]=size(IMF_all);%ch表示变量数，前三个表示模拟序列X,Y,Z,第4:18个表示高斯白噪声
w=[]; % Wasserstein distance between signal and noise
 intv=101:900; % remove end effect of MEMD,about10%
% intv=51:450;
for j=1:nimf
    for i=1:ch_n-1
        for ii=ch_n:ch 
            y1=squeeze(IMF_all(i,j,intv));%删除数组中维度为1的维度
            y2=squeeze(IMF_all(ii,j,intv));
            w(j,i,ii)=wadist(y1,y2);        
        end
    end
end
w=w(:,:,ch_n:end);

wn=[];% Wasserstein distance between noises
for j=1:nimf
    tmp=[];
    for i=ch_n:ch-1
        for ii=i+1:ch
            y1=squeeze(IMF_all(i,j,intv));
            y2=squeeze(IMF_all(ii,j,intv));
            tmp=[tmp wadist(y1,y2)];        
        end
    end
    wn(j,:)=tmp;  
end

%% Identification%%%%%%%%%%%5
alp=0.05; %% significance level, with alpha = 0.05
    
ww=squeeze(mean(w,3)); %% distance between signal and noise
wn=sort(wn,2);%% distance between noise and noise
th=wn(:,fix(size(wn,2)*(1-alp/2))); %% CI=[wn(:,fix(size(wn,2)*(alp/2))) wn(:,fix(size(wn,2)*(1-alp/2)))]
th=repmat(th,1,ch_n-1);
r=th<ww; % IMFs are identified in 'r' variable (1-significant; 0 - otherwise) 