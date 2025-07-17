clear;
clc;

load test_typeI_error_c2
%%load test_power_c2
%%%%%%%%%%determine null distribution%%%%%%%%%
for aa=1:length(En_all)
    bb1=size(En_all(aa).En_IMF,2);
    bb2=size(En_w_all(aa).En_w,2);
    cc1(aa)=bb1;cc2(aa)=bb2;
end
   L_use1=find(cc1==10);L_use2=find(cc2==10);
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
    p1=p1(1:9,:);
    p2=E_IMF_all(cf).E_IMF';
    p2=p2(1:9,:);
    pn(:,cf)=p1(:);%%%%噪声的能量
    px(:,cf)=p2(:);%%%三个变量的能量
end

%%% choose the sufficient size M%%%%%%%%%%%%%%%%%
[M1,N1]=size(px);M=1000;%%%M=[1000 3000 5000 7000 9000 10000];
alp=0.05;%%%%%%%%%%alp=[0.01 0.05 0.1 0.15 0.2 0.25];
pn_1=pn(:,1:M);
En=sort(pn_1,2);
q_cv=En(:,fix(size(En,2)*(1-alp/M1)));%% threshold for test
th=repmat(q_cv,1,N1);
r=px>th;
% % %%%%%%%%%%%%%%%%%%%calculate type 1 error%%%%

n11_all=0;
n12_all=0;
for jj=1:N1
     px11=r(:,jj);
     n1=sum(px11);
     if n1>0;
        n11_all=n11_all+1;
     else
     end
end
err1=n11_all/N1
 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate power%%%%%%%%

% n11_all=0;
% n12_all=0;n22_all=0;
% for jj=1:N1
%     px11=r(:,jj);
%     nimf=9;
%     indxsig=[4 5 6 4+nimf 6+nimf 4+2*nimf 5+2*nimf];
%     n1=sum(px11(indxsig(1:3))); n2=sum(px11(indxsig(4:5)));
%     n3=sum(px11(indxsig(6:7)));
%     if n1>0&&n2>0&&n3>0;
%         n11_all=n11_all+1;
%     else
%     end
% end
% power2=n11_all/N1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


