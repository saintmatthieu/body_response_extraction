% Frequency-Domain String Extraction (FDSE)
%     E_ = FDSE(X_,Y_,m_,a_,N,AC_FF) produces the Short-Time Fourier
%     Transform (STFT) E_ of the residual after FDSE of the STFT X_.

%     Y_ is the STFT of x_ forward one sample, i.e. X_=STFT(x_(1:end-1)),
%     while Y_=STFT(x_(2:end)). m_ is the vector of time indices where the
%     short-time spectra of X_ were evaluated. a_ and N are the coefficient
%     vector and the length of the cosine window used in the STFT. AC_FF is
%     a preliminary estimate of the pitch of the initial time-domain
%     signal.

function E_ = FDSE(X_,Y_,m_,a_,N,AC_FF)

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

E_ = X_;
M2pi = 2*pi/M;

for u = 1:U
    % To monitor progression
    disp([num2str(u) '/' num2str(U)])
    
    % Median-Adjustive Trajectories (MAT) for the Fundamental Frequency and Inharmonicity Coefficient estimations
    % (3.3.1)
    [FF,IC] = MAT(conj(X_(:,1,u)).*Y_(:,1,u),AC_FF);
    
    for c = 1:chans
        % Algorithm for the integrated detection and cancelation of transverse and phantom partials : initialisation
        % (3.3.3)
        k = 1; l = 1;
        transverse = 0;
        phantom = 0;
        
        % Magnitude and curvature spectra for peak and bulge detections
        % (3.3.3)
        M_ = abs(E_(:,c,u));
        C_ = [NaN;diff(M_,2);NaN];
        Mp_ = find(M_(1:end-2)<M_(2:end-1)&...
            M_(2:end-1)>M_(3:end))+1;
        Cb_ = find(C_(1:end-2)>C_(2:end-1)&...
            C_(2:end-1)<C_(3:end)&C_(2:end-1)<0)+1;
        CB = length(Cb_);
        Cbi = 1; Cbj = 1;
        thre = max(find(M_(1:M/2)>N*10^(.05*thre_dB)))-1;
        
        while 1
            
            % Estimated frequencies of next transverse and phantom partial
            % (1.2.7)
            fk = k*FF*sqrt(1+IC*k^2);
            gl = l*FF*sqrt(1+.25*IC*l^2);
            
            if fk>thre | gl>thre, break; end
            
            % Find closest bulges
            [fb,Cbi]=closest(Cb_,fk+1,Cbi);
            [gb,Cbj]=closest(Cb_,gl+1,Cbj);
            
            % Check if they are also peaks
            [fIsPeak,fp] = isPeak(M_,fb);
            [gIsPeak,gp] = isPeak(M_,gb);
            
            % In case the peak has already been canceled
            flag = 1;
            if ~sum(fp==Mp_)
                k = k+1;
                flag = 0;
            end
            if ~sum(gp==Mp_)
                l = l+1;
                flag = 0;
            end
            
            if flag 
                
                % Algorithm described Section 3.3.3
                if fIsPeak & gIsPeak & fp~=gp
                    transverse = fp<gp;
                    phantom = ~transverse;
                    p = fp*transverse+gp*phantom;
                end
                if fIsPeak & gIsPeak & fp==gp
                    p = fp;
                    transverse = 1;
                    phantom = 1;
                end
                if xor(fIsPeak,gIsPeak)
                    transverse = fIsPeak;
                    phantom = gIsPeak;
                    p = fp*transverse+gp*phantom;
                end
                if ~fIsPeak & ~gIsPeak
                    transverse = 0;
                    phantom = 0;
                    p = fp;
                end
                
                if M_(p)>N*10^(.05*thre_dB)
                    % Peak measurement and cancelation
                    if fIsPeak | gIsPeak
                        
                        % CSPME measurement of 1st-order terms
                        % (2.4.3)
                        C = Y_(p,c,u)/E_(p,c,u);
                        omg = imag(log(C));
                        gam = real(log(C));
                        
                        % Cancelation of partial
                        % (2.4.3), (3.2.2), (4.3)
                        bb_ = p+b_-1;
                        nu = max(m_(u),0)-m_(u);
                        if nu, gam=0; end
                        zet = gam+1i*(omg-(p-1)*M2pi);
                        GHV_ = sum(A.*(exp(N*(OMG+zet))-exp(nu*(OMG+zet)))./(OMG+zet),2);
                        Z = E_(p,c,u)/GHV_(B2+1);
                        E_(bb_+1,c,u) = E_(bb_+1,c,u)-GHV_*Z;
                        Y_(bb_+1,c,u) = Y_(bb_+1,c,u)-GHV_*Z*C;
                        
                        % Remove p from list of peaks
                        Mp_(Mp_==p) = [];
                        
                    end
                end
                
                % Increment harmonic numbers appropriately
                if transverse & ~ phantom
                    k = k+1;
                    if fp>climbUp(M_,gb) l = l+1; end
                elseif phantom & ~transverse
                    l = l+1;
                    if gp>climbUp(M_,fb) k = k+1; end
                else
                    k = k+1;
                    l = l+1;
                end
                
                % Update spectrum if there was a situation of overlap
                if xor(transverse,phantom) & climbUp(M_,fp)==climbUp(M_,gp)
                    M_ = abs(E_(:,c,u));
                    C_ = [NaN;diff(M_,2);NaN];
                    Mp_ = find(M_(1:end-2)<M_(2:end-1)&...
                        M_(2:end-1)>M_(3:end))+1;
                    Cb_ = find(C_(1:end-2)>C_(2:end-1)&...
                        C_(2:end-1)<C_(3:end)&C_(2:end-1)<0)+1;
                    CB = length(Cb_);
                    Cbi = 1; Cbj = 1;
                end
            end
        end
    end 
end

E_ = cat(1,E_(1:M/2+1,:,:),conj(E_(flipud(2:M/2),:,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peak,Mb] = isPeak(M_,Mb);
M = length(M_);
C_ = [NaN;diff(M_,2);NaN];

peak = 0;
Mb = round(Mb);
flag = 1;
while flag
    if Mb<=1 | Mb>=M
        warning('Reached end of spectrum before finding peak.');
        return
    end
    % See if magnitude maximum, or "peak", is reached
    if M_(Mb-1)<M_(Mb) & M_(Mb)>M_(Mb+1)
        peak = 1;
        return
    end
    % See if curvature maximum, or "cave" is reached
    if C_(Mb-1)<C_(Mb) & C_(Mb)>C_(Mb+1) & C_(Mb)>0
        flag = 0;
    end
    % Climb up
    Mb = Mb+(-1)^(M_(Mb-1)>M_(Mb+1));
end

while ~(M_(Mb-1)<M_(Mb) & M_(Mb)>M_(Mb+1))
    Mb = Mb+(-1)^(M_(Mb-1)>M_(Mb+1));
    if Mb==1|Mb==M
        warning('Reached end of spectrum before finding peak.');
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Xp = climbUp(X_,Xi)

N = length(X_);
if Xi==1 | Xi==N
    warning('initial index on edge of vector')
    Xp = Xi;
    return
end

while ~(X_(Xi-1)<X_(Xi) & X_(Xi)>X_(Xi+1))
    Xi = Xi+(-1)^(X_(Xi-1)>X_(Xi+1));
    if Xi==1 | Xi==N
        warning('end of vector reached before peak was found')
        Xp = Xi;
        return
    end
end

Xp = Xi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Xi]=closest(X_,x,Xi)
% Returns value X of sorted vector X_ closest to x,
% starting search in X_ from index Xi and going up.

while Xi<length(X_) & X_(Xi)<x, Xi = Xi+1; end
if Xi>1, Xi = Xi - (X_(Xi)+X_(Xi-1)>2*x); end
X = X_(Xi);
