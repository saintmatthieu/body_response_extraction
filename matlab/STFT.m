% Short-Time Fourier Transform (STFT)
%     [X_,m_] = STFT(x_,w_,R,M) produces the STFT X_ of vector x_, where w_
%     is the analysis window, R the hop size, and M, the FFT length. m_ is
%     the vector of time indices corresponding to the first sample of each
%     short-time spectrum.

%     The function was written to support multi-channel input, the number
%     of columns in x_ being interpreted as the number of channels.
%     The size of X_ will therefore be (M,chans,U), where chans is the
%     number of channels, and U, the minimal number of windows for the
%     transparent reconstruction of x_ upon inverse STFT - provided that
%     w_ and R were chosen in a way to make transparent reconstruction
%     possible.

function [X_,m_] = STFT(x_,w_,R,M)

[in_N,chans] = size(x_);
ww_ = repmat(w_,1,chans);
N = length(w_);
if nargin<4, M = 2^ceil(log2(N)); end

O = N/R;
if mod(O,1), error('Window length must be multiple of hop size'); end
U = ceil(in_N/R)+O-1;
u_ = (1:U)';
m_ = (u_-O)*R;
l_Z = (O-1)*R;
r_Z = (O+floor(in_N/R))*R-in_N;

zxz_ = [zeros(l_Z,chans);x_;zeros(r_Z,chans)];
ZNZ = in_N+l_Z+r_Z;

M2pi = M/2/pi;
X_ = zeros(M,chans,U);
for u=u_'
    n_ = m_(u)+(0:N-1)';
    N_x_ = zxz_(n_+l_Z+1,:);
    X_(:,:,u) = fft(N_x_.*ww_,M);
end
