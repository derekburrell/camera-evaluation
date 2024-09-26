function[u0]=propF_TF(u1,L,wvl,Z)
% propagation: transfer-function approach
% Fresnel transfer function from Voelz - Eq. (5.2)
% assumes same x and y side lengths and uniform sampling
% u1: source-plane field
% L: source- and observation-plane side length [m]
% wvl: optical wavelength [m]
% Z: propagation distance [m]
% u0: observation-plane field

[M,N] = size(u1);                       % input-field array size [px]
k = 2*pi/wvl;                           % angular wavenumber [rad/m]

fx = (-M/2:M/2-1)*1/L;                  % x-frequency coordinates
fy = (-N/2:N/2-1)*1/L;                  % y-frequency coordinates
[FX,FY] = meshgrid(fx,fy);              % frequency coordinates

H = exp(1j*k*Z) ...                     % define transfer function
    *exp(-1j*pi*wvl*Z*(FX.^2+FY.^2));
H = fftshift(H);                        % shift transfer function
U1 = fft2(fftshift(u1));                % shift FFT source field
U0 = H.*U1;                             % multiply
u0 = ifftshift(ifft2(U0));              % inv. FFT, center obs. field

end