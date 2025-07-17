clear;
clc;

load rossler_dec_M_5000
% % load rossler_1000
% % tran_num=1000;Y = x1(1:tran_num,3);
% % % Y=y2;t = length(Y);
% % % Y=Y/std(Y);
% % t = length(Y);
% % allmode = eemd(Y,0.3,200);%EEMD decompose
% % % % allmode = eemd(Y,0,1);
% % allmode_no_X = allmode(:,2:end);
% % allmode_no_X_T = allmode(:,2:end-1);
% % % % allmode_trend = allmode(:,end);
% % % % save EEMD_1095_1 allmode_trend
% % % % allmode = emd(Y);allmode=allmode';%EMD decompose
% % figure(1);  
% % nimf = size(allmode_no_X,2); 
% % for m = 1:nimf
% %   subplot(nimf,1,m); plot(1:t,allmode_no_X(:,m),'k','LineWIdth',1.5);
% % end
% % subplot(nimf,1,1); title('IMFs');
for cf=1:length(imf_all)
    IMF(:,:,:,cf)=imf_all(cf).imfx;
end
IMF_all=mean(IMF(1:2,:,:,:),4);
IMF_all1=permute(IMF_all,[3 2 1]);
allmode_no_X=IMF_all1(:,:,2);
[Npt,Nimf]=size(allmode_no_X);
logep = significanceIMF(allmode_no_X(:,1:Nimf-1));
sigline95_up = confidenceLine(0.05,Npt);
sigline95_low = confidenceLine(0.95,Npt);
[m1,n1] = size(logep);[m2,n2] = size(sigline95_up);
 figure
 b_text = {'c1','c2','c3','c4','c5','c6','c7','c8','c9'};
 plot(sigline95_low(:,1),sigline95_low(:,2),'r',sigline95_up(:,1),sigline95_up(:,2),'b');
% % legend('90% significance','95% significance','Location', 'SouthWest');
hold on
plot(logep(:,1),logep(:,2),'g*');
text(logep(:,1)+0.04,logep(:,2)+0.04,b_text)
xlabel('Log2(Mean period)','FontSize',8,'VerticalAlignment','middle');ylabel('log2(Mean Normalized Energy)','FontSize',8);
