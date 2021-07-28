% ----------------------------------------------------------------------------------------------------------------
% FUNCTION simpl_FRS (f,ground,eig_freq,Mode_shapes,p_fact,damp_ratio_struc,xi,fc)
% ----------------------------------------------------------------------------------------------------------------
% ######################################################################
% # Simplified FRS generated as a function of basic modal parameters:
% # natural frequencies, modal shapes, participation factors.
% # Extra characteristics needed are the damping ratios of the structue
% # and of the equipment.
% #
% #
% # variable            I/O       meaning                            units
% # -------------------------------------------------------------------------------------------------------------
% # f (vector)          I       frequency                            Hz
% # ground (vector)     I       ground spectrum                      m/s^2
% # eig_freq (vector)   I       natural freq. structure              Hz
% # Mode_shapes         I       modal shapes                         m
% # p_fact              I       participation factors 
% # damp_ratio_struc    I       damping ratio structure              decimal
% # xi                  I       damping ratio equipment              decimal
% # fc                  I       upper corner freq. of ground spectr. Hz
% # FRS                 O       simplified FRS                       m/s^2
% # Author : Simone Kaldas (Fusion for Energy)
% # Version: 1.0 June 2021
% # F4E UID: F4E_D_2R9EE2
% ######################################################################


function FRS = simpl_FRS(f,ground,eig_freq,Mode_shapes,p_fact,damp_ratio_struc,xi,fc)

gamma = ground(end); %ground ZPA

n_modes = numel(eig_freq);

eigfreq_index = []; %position of the natural frequencies in the frequency vector f
for n =1:n_modes
    [~,eigfreq_index(n)] = min(abs(f-eig_freq(n)));

    Sa_n(n) = ground(eigfreq_index(n));  %Spectral acceleration at the natural frequencies
    Sa(n) = Sa_n(n)*p_fact(n)*Mode_shapes(n);
    miss_mass(n) = p_fact(n)*Mode_shapes(n);
end

Ssa =  sqrt(gamma^2*(1-sum(miss_mass)).^2 + sum(Sa.^2));

f1 = eig_freq(1);
fn = eig_freq(end);
flim = min(33,1.2*fn); 

n1 = sqrt(10/(5 + damp_ratio_struc*100)); % damping factor
n2 = sqrt(5/(xi*100));
n = n1*n2;

int1 = f > 0 & f < 0.8*f1;
int2 = f < 1.2*f1 & f >= 0.8*f1;
int3 = f <= flim & f >= 1.2*f1;
int4 = f>flim;

fc1 = 1.5*fc;
fc2 = 3*fc;

if eig_freq(1)<=fc1
    c = 5;
elseif eig_freq(1)>=fc2
    c = 1;
elseif eig_freq(1)>fc1 && eig_freq(1)<fc2
    c = 5 + (5-1)./(fc1-fc2).*(eig_freq(1)-fc1); 
end

Kt(int1) = c.*n ./ ((0.8*f1./f(int1)).^2);
Kt(int2) = c.*n;
Kt(int3) = c.*n - ((c.*n-1).*log(1.2*f1./f(int3))/log(1.2*f1/flim));
Kt(int4) = 1;

FRS = Ssa .* Kt;

[~,i] = max(ground);
fpeak = f(i(1));

hp1 = f<0.8*f1;
hp2 = f>=0.8*f1;

alpha(hp1)=1-f(hp1)./(0.8*f1);
alpha(hp2)=0;
alpha = alpha';

FRS(f<0.8*f1) = alpha(f<0.8*f1).*ground(f<0.8*f1)+(1-alpha(f<0.8*f1))*FRS(f==round(f1));

%If no isolation, impose FRS >= ground spectrum
if f1 > fpeak 
    for N = 1:length(f)
        if FRS(N) < ground(N)
            FRS(N) = ground(N);
        end
    end
end
end