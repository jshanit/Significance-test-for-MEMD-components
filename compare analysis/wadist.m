function w1=wadist(x1,x2)

%%% Implementation of Wasserstein distance between x1 and x2

%% input: x1, x2  - two vectors
%% output: w1 : Wasserstein distance

%%% Meng Hu @ Drexel University, 2013

x1=zscore(x1);
x2=zscore(x2);

% Wasserstein distance
[f1,ix1]=ecdf(x1); 
[f2,ix2]=ecdf(x2); 

if length(f1)==length(f2)
    w1=sum(abs(ix1-ix2))*(1/length(f1));
elseif length(f1)>length(f2)
    ix2p=interp1(f2,ix2,f1); 
    w1=sum(abs(ix1-ix2p))*(1/length(f1));
else
    ix1p=interp1(f1,ix1,f2); 
    w1=sum(abs(ix1p-ix2))*(1/length(f2));
end


end