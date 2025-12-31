function [E,freq,v,wtG] =  MeasurePVByTwoPTCWT_G(CCF, r,dt,FreqRange, VelocityRange, wavename)
%   Summary of this function goes here.
%   [E,freq,v,wtG] =  MeasurePVByTwoPTCWT_G(CCF, r,dt,FreqRange, VelocityRange, wavename)
%   Detailed explanation goes here.
%   The function is for measuring the phase velocity of dispersive waves from two traces.
%
%   IN      
%           CCF: the cross-correlation or Green's function of two traces. 
%             r: the distance (m) of A and B.
%            dt: the sampling interval in time domain (s).
%     FreqRange: the range of frequency (Hz) in the dispersive energy, such
%                as [5 100 1].
% VelocityRange: the parameter of the phase velocity (m/s), first is the 
%                minimal, second is the maximal, third is the interval, such 
%                as [50 800 1].
%      wavename: the type of wave, 'cmor1-1', 'db1','harr', and so on.
%
%  OUT   
%             E: the matrix of the normalized dispersive energy.
%          freq: the frequency (Hz) vecotor of the normalized dispersive energy.
%             v: the phase velocity (m/s) vector of the normalized dispersive energy.
%           wtA: the cwt result of Green's function.
%  
%  References: 
%    Kijanka P , Ambrozinski L , Urban M W . Two Point Method For Robust Shear Wave 
%    Phase Velocity Dispersion Estimation of Viscoelastic Materials[J]. Ultrasound 
%    in Medicine & Biology, 2019, 45(9):2540-2553.
%
%  Author(s): Yan Yingwei
%  Copyright: 2020-2025 
%  Revision: 1.0  
%  Date: 10/5/2020
%
%  Department of Earth and Space Sciences, Southern University of Science 
%  and Technology.
%'haar': Haar 小波
% 'dbN': Daubechies 小波，其中 N 是小波的长度。例如，'db2'、'db4' 等。
% 'symN': Symlets 小波，其中 N 是小波的长度。例如，'sym2'、'sym4' 等。
% 'coifN': Coiflets 小波，其中 N 是小波的长度。例如，'coif1'、'coif2' 等。
% 'bior': Biorthogonal 小波，其中可能包括不同的子波对。例如，'bior1.3'、'bior2.2' 等。
% 'dmey': Demeyer 小波（也称为 Meyer 小波）
% 'morl': Morlet 小波
% 'cgau': Cauchy 小波


% 输入参数判断
if nargin==5
    wavename = 'cmor1-1';
elseif nargin==4
    wavename = 'cmor1-1';
    VelocityRange = [50 2000 5];
elseif nargin==3
    wavename = 'cmor1-1';
    VelocityRange = [50 2000 5];
    FreqRange = [5 100 1];
end

% 确定基本参数
vmin = VelocityRange(1);
vmax = VelocityRange(2);
dv = VelocityRange(3);
v = vmax:-dv:vmin;
lv = length(v);

fmin = FreqRange(1);
fmax = FreqRange(2);
df = FreqRange(3);
freq = fmin:df:fmax;
lf = length(freq);

E = zeros(lv, lf);
M  = length(CCF);

Fs = 1/dt;      % 信号的采样率

t = dt:dt:M*dt;  % 信号的时间轴向量
ti = r./v;       % 由速度范围确定的互相关函数的时间范围

% 对信号作小波变换
wave_centfreq =centfrq(wavename);
scales = Fs*wave_centfreq./freq;
wtG = cwt(CCF,scales,wavename);

for i=1:lf
    E(:,i) = interp1(t, real(wtG(i,:)),ti,'PCHIP');
    E(:,i) = E(:,i).^(2).* E(:,i);
    E(:,i) = E(:,i)./max(abs(E(:,i)));
%     E(:,i) = (E(:,i)-min(E(:,i)))./max((E(:,i)-min(E(:,i))));
end
end

