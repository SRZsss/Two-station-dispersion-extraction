function [dp,Wx1]=GpstTir1(data,dt,f,v,r,m)
[Wx,fs,dWx,as] = GPST_fwm(data,dt,m,f);
% [Wx,fs,dst] =S_transfom(data,dt);
T = dt:dt:dt*length(data);
Ti = r./v;
Wx1 = imag(Wx);
dp1 = zeros(length(v),length(fs));
dp = zeros(length(v),length(f));
for i=1:length(fs)
    Ti = r./v+1/4/fs(i);
    dp1(:,i) = interp1(T,Wx1(i,:),Ti);
end
for i = 1:length(Ti)
    dp(i,:) = interp1(fs,dp1(i,:),f);
end
for i=1:length(f)
    dp(:,i) = (dp(:,i)+abs(min(dp(:,i))))/max(dp(:,i)+abs(min(dp(:,i))));
end
end