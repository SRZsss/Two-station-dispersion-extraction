function [E, freq, v] = NBFMeasureSWPV(CCF,r, dt,FreqRange,VelocityRange, h,N)
%   Summary of this function goes here.
%   [E, freq, v] = NBFMeasureSWPV(CCF,r, dt,FreqRange,VelocityRange, h,N)
%   Detailed explanation goes here.
%   The function is for measuring the phase velocity of dispersive waves
%   from interstation cross-correlation function (CCF) based on the narrow
%   band-pass filter.
%
%   IN      
%           CCF: the cross-correlation or Green's function of two traces. 
%             r: the propagation distance (m) of the incoming waves for
%                the two-station.
%            dt: the sampling interval in time domain (s).
%     FreqRange: the range of frequency (Hz) in the dispersive energy, such
%                as [1 10 0.25].
% VelocityRange: the parameter of the phase velocity (m/s), first is the 
%                minimal, second is the maximal, third is the interval, such 
%                as [500 2500 5].
%             h: the half-width (Hz) of the narrow band.
%             N: the order of the filter, Preferably greater than 100.
%
%  OUT   
%             E: the matrix of the normalized dispersive energy.
%          freq: the frequency (Hz) vecotor of the normalized dispersive energy.
%             v: the phase velocity (m/s) vector of the normalized dispersive energy.
%
%  References:
%  Yao, H., van Der Hilst, R. D., & De Hoop, M. V. (2006). Surface-wave array 
%  tomography in SE Tibet from ambient seismic noise and two-station analysis¡ªI. 
%  Phase velocity maps. Geophysical Journal International, 166(2), 732¨C744.
%  
%  Author(s): Yan Yingwei
%  Copyright: 2020-2025 
%  Revision: 1.0  
%  Date: 10/5/2020
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology.

if nargin==6
    N = 1000;
elseif nargin==5
    N = 1000;
    h = FreqRange(3)/25;
end

vmin = VelocityRange(1);
vmax = VelocityRange(2);
dv = VelocityRange(3);
v = vmax:-dv:vmin;
lv = length(v);

fmin = FreqRange(1);
fmax = FreqRange(2);
df = FreqRange(3);
M = length(CCF);

freq = fmin:df:fmax;
lf = length(freq);

E = zeros(lv,lf);

t = dt:dt:M*dt;
for j=1:lf
    f0 = freq(j);
    [s,~] = bpf(CCF,dt,N,[f0-h f0+h]);
    ti = r./v;
    E(:,j) = interp1(t,s,ti,'PCHIP');
    E(:,j) = E(:,j).^(2).* E(:,j);
    E(:,j) = E(:,j)./max(abs(E(:,j)));
end
for i=1:length(freq)
    E(:,i) = (E(:,i)+abs(min(E(:,i))))/max(E(:,i)+abs(min(E(:,i))));
end
% E = flip(E,1);
end

