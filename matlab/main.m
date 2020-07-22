close all, clear all

% Input waveform
[x_,SR]=wavread('sfname.wav');

chans = size(x_,2);

% Onset sample index, determined manually
% (3.2.2)
nu = 0;
x_ = x_(nu+1:end,:);

% Autocorrelation-based preliminary pitch estimate
% (Section 3.1.1 of my thesis)
in_N = 2^floor(log2(size(x_,1)));
xx_ = real(ifft(abs(fft(x_(:,1),in_N)).^2));
xx_(~([0;xx_(1:end-2)<xx_(2:end-1)&xx_(2:end-1)>xx_(3:end);0])) = 0;
[dump,T] = max(xx_(1:in_N/2));
T = T-1;

% 2-nd order continuous window
% (3.1.3)
P = 3;
i = 1;
O = P+i;
B = 2*(P+1);
j = 2;
R = round(j*B*T/O);
N = O*R;

% Pth-oder continuous window and its coefficients
% (2.2.3)
[w_,a_] = cosineWindow(N,P,1);

% Short-Time Fourier Transforms
% (2.4.3)
M = 2^ceil(log2(N));
[X_,m_] = STFT(x_(1:end-1,:),w_,R,M);
Y_ = STFT(x_(2:end,:),w_,R,M);

% String Extraction
u0 = O-1;
u_ = u0:size(X_,3)-O;
E_ = FDSE(X_(:,:,u_),Y_(:,:,u_),m_(u_),a_,N,M/T);
EE_ = FDSEE(X_(:,:,u_),Y_(:,:,u_),m_(u_),a_,N,M/T);

e_ = ISTFT(cat(3,X_(:,:,1:u0-1),E_),m_(1:end-O),N);
ee_ = ISTFT(cat(3,X_(:,:,1:u0-1),EE_),m_(1:end-O),N);