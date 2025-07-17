
clear;clc;

%%% time=0.001:0.001:1;xlim = [0.001,1];dt=0.001;
% time=0.001:0.001:1;xlim = [0.001,1];dt=0.001;
% a1=sin(2*pi*12*time);
% a2=sin(2*pi*26*time);
% a3=sin(2*pi*50*time);
% c=0.5;%the standard deviation of the input data 
% x11=a2+a3+c*randn(size(a1))+a1;
% x22=a3+c*randn(size(a1))+a1;
% x33=a2+c*randn(size(a1))+a3;
% sst=x33;
% n = length(sst);

% load rossler_1000
% time=1:1000;xlim = [1,1000];dt=1.5;%%%%%%dt=1.5 for red noise test
% sst=x1(:,2);
% n = length(sst);

% load rossler_1000
% time=1:1000;xlim = [1,1000];dt=1.5;%%%%%%dt=1 for white noise test
% sst=x1(:,2);
% n = length(sst);

% load SN_SF_CR_1958_2009_month_data
% time=1:620;xlim = [1,620];
% dt=1.8;%%%%%%%%%%dt=1.8 for red noise test
% sst=x1(:,3);
% n = length(sst);

load SN_SF_CR_1958_2009_month_data
time=1:620;xlim = [1,620];
dt=1;%%%%%%%%%%for white noise test
sst=x1(:,3);
n = length(sst);
%------------------------------------------------------ Computation

variance = std(sst)^2;
sst = (sst - mean(sst))/sqrt(variance) ;

% n = length(sst);
% dt = 0.25 ;
% time = [0:length(sst)-1]*dt + 1871.0 ;  % construct time array
% xlim = [1870,2000];  % plotting range
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 2*dt;    % this says start at a scale of 6 months
j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag1=0;
% lag1 = 0.72;  % lag-1 in range of [0,1], autocorrelation for red noise background,if lag1=0 for white noise
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% % % Global wavelet spectrum & significance levels:
% global_ws = variance*(sum(power')/n);   % time-average over all times
% dof = n - scale;  % the -scale corrects for padding at edges
% global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);


% %%%%--- Contour plot wavelet power spectrum
subplot('position',[0.1 0.37 0.65 0.28])
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
contour(time,log2(period),log2(power),log2(levels));  %*** or use 'contourfill'
%imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
xlabel('t')
ylabel('Period')
title('Wavelet Power Spectrum')
set(gca,'XLim',xlim(:))
set(gca,'YLim',log2([min(period),max(period)]), ...
    'YDir','reverse', ...
	'YTick',log2(Yticks(:)), ...
	'YTickLabel',Yticks)
colorbar
% 95% significance contour, levels at -99 (fake) and 1 (95% signif)
hold on
contour(time,log2(period),sig95,[-99,1],'k');
hold on
% cone-of-influence, anything "below" is dubious
plot(time,log2(coi),'k')
hold off

% subplot('position',[0.1 0.37 0.65 0.28])
% levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16] ;
% Yticks = 2.^(fix(log2(min(1./period))):fix(log2(max(1./period))));
% contour(time,log2(1./period),log2(power),log2(levels));  %*** or use 'contourfill'
% %imagesc(time,log2(period),log2(power));  %*** uncomment for 'image' plot
% xlabel('t')
% ylabel('Frequency')
% 
% set(gca,'XLim',xlim(:))
% set(gca,'YLim',log2([min(1./period),max(1./period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks)
% colorbar
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% contour(time,log2(1./period),sig95,[-99,1],'k');
% hold on
% % cone-of-influence, anything "below" is dubious
% plot(time,log2(1./coi),'k')
% hold off


% end of code

