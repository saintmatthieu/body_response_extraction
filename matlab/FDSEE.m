% Frequency-Domain String Extraction (FDSE)
%     E_ = FDSEE(X_,Y_,m_,a_,N,AC_FF) produces the Short-Time Fourier
%     Transform (STFT) E_ of the residual after FDSE of the STFT X_.
%
%     Y_ is the STFT of x_ forward one sample, i.e. X_=STFT(x_(1:end-1)),
%     while Y_=STFT(x_(2:end)). m_ is the vector of time indices where the
%     short-time spectra of X_ were evaluated. a_ and N are the coefficient
%     vector and the length of the cosine window used in the STFT. AC_FF is
%     a preliminary estimate of the pitch of the initial time-domain
%     signal.
% 
%     The function FDSEE is a simplified version of FDSE. FDSEE is
%     computationally cheaper and is more concise, but phantom partials, if
%     any, are ignored. Low-pitch tones should generally use FDSE, but
%     high-pitched tones, where phantom partials are negligible, should be
%     treated with FDSEE.

function X_ = FDSEE(X_,Y_,m_,a_,N,AC_FF)

% Peak detection threshold
thre_dB = -96;

[M,chans,U] = size(X_);
P = length(a_)-1;

% Construction of the matrices for computational-efficient spectral synthesis
% (4.3)
B = 2*(P+1)*M/N;
B2 = floor(B/2);
b_ = (-B2:B2)';
p_ = (-P:P);
OMG = -1i*2*pi*(1/M*b_*ones(1,length(p_))+1/N*ones(length(b_),1)*p_);
aa_ = .5*(-1).^p_.*a_(abs(p_)+1)';
aa_(P+1) = a_(1);
A = ones(length(b_),1)*aa_;

M2pi = 2*pi/M;

for u = 1:U
    % To monitor progression
    disp([num2str(u) '/' num2str(U)])
    
    % Median-Adjustive Trajectories (MAT) for the Fundamental Frequency and Inharmonicity Coefficient estimations
    % (3.3.1)
    [FF,IC] = MAT(conj(X_(:,1,u)).*Y_(:,1,u),AC_FF);
    
    for c = 1:chans
        
        % Magnitude and curvature spectra for peak and bulge detections
        % (3.3.3)
        M_ = abs(X_(:,c,u));
        C_ = [NaN;diff(M_,2);NaN];
        Mp_ = find(M_(1:end-2)<M_(2:end-1)&M_(2:end-1)>M_(3:end))+1;
        MP = length(Mp_);
        Mpi = 1;
        thre = max(find(M_(1:M/2)>N*10^(.05*thre_dB)))-1;
        
        for k=1:200
            
            % Estimated frequencies of next transverse partial
            % (1.2.7)
            fk = k*FF*sqrt(1+IC*k^2);
            
            if fk>thre, break; end
            
            % Find closest peak
            [p,Mpi]=closest(Mp_,fk+1,Mpi);            
            
            if M_(p)>N*10^(.05*thre_dB)
                
                % CSPME measurement of 1st-order terms
                % (2.4.3)
                C = log(Y_(p,c,u)/X_(p,c,u));
                omg = imag(C);
                gam = real(C);
                
                % Cancelation of partial
                % (2.4.3), (3.2.2), (4.3)
                bb_ = p+b_-1;
                nu = max(m_(u),0)-m_(u);
                if nu, gam=0; end
                zet = gam+1i*(omg-(p-1)*M2pi);
                GHV_ = sum(A.*(exp(N*(OMG+zet))-exp(nu*(OMG+zet)))./(OMG+zet),2);
                X_(bb_+1,c,u) = X_(bb_+1,c,u)-GHV_*X_(p,c,u)/GHV_(B2+1);
            end
        end
    end
end

X_ = cat(1,X_(1:M/2+1,:,:),conj(X_(flipud(2:M/2),:,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Xi]=closest(X_,x,Xi)
% Returns value X of sorted vector X_ closest to x,
% starting search in X_ from index Xi and going up.

while Xi<length(X_) & X_(Xi)<x, Xi = Xi+1; end
if Xi>1, Xi = Xi - (X_(Xi)+X_(Xi-1)>2*x); end
X = X_(Xi);
