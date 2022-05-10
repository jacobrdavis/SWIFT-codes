function [f,Y,Ma,Ph] = fast_fourier_transform(y,fs)

Y = fft(y);
L = length(y);

% magnitude
Ma2 = abs(Y/L);
Ma = Ma2(1:round(L/2+1));
Ma(2:end-1) = 2*Ma(2:end-1);

% phase
Ph2 = angle(Y);
Ph = Ph2(1:round(L/2+1));

% frequency
f = fs*(0:round(L/2))/L;
end