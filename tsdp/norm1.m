function data = norm1(data)
%% Normalized function
% Calculate the result of vector normalization for each column
% 
% Input :
% data : Input vector
% 
% Output : 
% data : Vector after normalization
% 
% Author(s) : CSIM
% Copyright : Open source
% Revison : 1.0  Date : 16/11/2022
% 
% CSIM,Key Laboratory of Applied Geophysics, Ministry of Natural Resources
% College of geoexploration science and technology,Jilin university of China
% 
for i = 1:length(data(1,:))
    if max(data(:,i))==0
        data(:,i) = data(:,i)/1;
    else
        data(:,i) = data(:,i)/(max(abs(data(:,i)))+1e-9);
    end
end
end