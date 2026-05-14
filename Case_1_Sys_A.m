%%%%%%%%%%%%%%%%%%%%%%%%
%Sys-A: ture FP and SP
%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

N     = 70000;
Nn    = N/2;
N_eva = 7000;
T     = 100;

% narrowband reference
Fs = 2000; % sampling fre
Amp  = [ 1   1   1   1   1 ] * sqrt(4/5);  % source amplitude (variance 2)
omega = 2*pi*[ 100 150 300 400 450 ]/Fs;
c_ture = -2*cos(omega);
[wk, q]= size(omega);

% 1st (primary) path P(z) 
M_p = 48;      
cutoff = 0.4*pi;    % lowpass filter cutoff fre
s_p = fir1(M_p-1, cutoff/pi);
   
% 2nd path 
M_s = 21; 
M_change = 2;
s_0 = fir1(M_s-1,cutoff/pi);

s_1 = fir1(M_s-1-M_change,cutoff/pi);

s = zeros(M_s, N);
for n=1:Nn
    s(:, n) = ones(M_s, 1) .* s_0';
end

for n=Nn+1:N
    s(1:M_s-M_change, n) = ones(M_s-M_change, 1) .* s_1';
end

% feedback path
M_fb = 32;
s_fb_0 = fir1(M_fb-1,cutoff/pi);

s_fb_1 = fir1(M_fb-1-M_change,cutoff/pi);

s_fb = zeros(M_fb, N);
for n=1:Nn
    s_fb(:, n) = ones(M_fb, 1) .* s_fb_0';
end

for n=Nn+1:N
    s_fb(1:M_fb-M_change, n) = ones(M_fb-M_change, 1) .* s_fb_1';
end

hat_s_fb = s_fb;
hat_M_fb = M_fb;
hat_s   = s;
hat_M_s = M_s;

%%%%%%%%%%%%%%%
soft_D = 1; %Delay
%%%%%%%%%%%%%%

% Source signal
x_clean = zeros(1, N);
for i = 1:q
    x_clean = x_clean + Amp(i) * sin(omega(i)* (0:N-1) + 0);
end

%%%%
% additive noise in the narrowband source
std_v_s = sqrt(0.001);  
FLG_v_s = 1;
if FLG_v_s == 0
    v_s_org = randn(100,120000);
    save v_s_org.mat  v_s_org;
else
    load v_s_org.mat;
end

%%%
% v_p(n)
std_v_p = sqrt(0.01);  % std of additive noise in p(n)
FLG_v_p = 1;
if FLG_v_p == 0
    v_p_org = randn(100, 120000);
    save v_p_org.mat   v_p_org;
else
    load v_p_org.mat; 
end

% W_N(z)   ---- FIR type
%%%%%%%%%%%%%%%%%%%
L_N = M_p - M_s + 5;
mu_N = 0.0003;  

%%%%%%%%
%%% results for ideal system (benchmark)
%%%%%%%%%
save_e_A  = zeros(1, N);
save_e2_A = zeros(1, N);
save_G_A  = zeros(1, N);

save_NRP_A_1st = zeros(1,T);
save_NRP_A_2nd = zeros(1,T);

tic;
for t = 1:T
    
    disp(' t =====> '); disp(t);
      
    % 3 additive noises
    v_p = v_p_org(t,1:N) * std_v_p;     % additive noise in p(n)
    v_s = v_s_org(t,1:N) * std_v_s;      % additive noise in noise source x_s(n)

    x_s = x_clean + v_s;  

    %%%
    % primary noise
    p_clean = zeros(1, N);
    for n=1:N
        for j=0:M_p-1
            if n-j>0
                p_clean(n) = p_clean(n) + s_p(j+1)*x_s(n-j);
            end
        end
    end
    p   = p_clean + v_p;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Ideal system
    %%%  NANC (using ture SP and FBP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %reference
    r    = zeros(1,N);
    x    = zeros(1,N);
    hat_x = zeros(1,N);
   
    %%%%%%%%%%%%%%%%%
    % paramaters of controller
    %%%%%%%%%%%%%%%%%
    y_0 = zeros(1,N);  % output of W_N(z)
    w_N = zeros(L_N, N); % paramaters of W_N(z)
     
    %%%%%%%%%%%%%%%%%
    % paramaters of FBPM
    %%%%%%%%%%%%%%%%%
    y_f        = zeros(1, N);
    hat_y_f    = zeros(1, N);
   
    y     = zeros(1,N); %2nd source
    y_p   = zeros(1,N); % secondary path after SP
    e     = zeros(1,N); %residual noise

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% NANC
     
    for n=1:N  % time loop start    
        
       % y(n)
       if n-1>0
           y(n) = y_0(n-soft_D);
       end
                
        % y_f : true FBP output
        for m = 0:M_fb-1
            if n-m > 0
                y_f(n) = y_f(n) + s_fb(m+1,n) * y(n-m);
            end
        end
        
        % r(n) : reference signal
        r(n) = x_s(n) + y_f(n);
        
        % hat_y_f : estimate of FBP signal
        for m = 0:hat_M_fb-1
            if n-m>0
                hat_y_f(n) = hat_y_f(n) + hat_s_fb(m+1,n) * y(n-m);
            end
        end
        
        % x(n) : reference input
        x(n) = r(n) - hat_y_f(n);
        
        %  Output of W_N(z)
        for j=0:L_N-1
            if n-j > -0
                y_0(n) = y_0(n) + w_N(j+1, n) * x(n-j);
            end
        end
             
        %y_p : signal output of 2nd source passing 2nd path
        for m=0:M_s-1
            if n-m>0
                y_p(n) = y_p(n) + s(m+1,n) * y(n-m);
            end
        end
       
        %e(n) : residual noise
        e(n) = p(n) - y_p(n);      
          
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Update preparing
        %%%%%%%%%%%%%%%%%%%%%%%
        % Preparations for update
        for m=0:hat_M_s-1
            if n-m>0
                hat_x(n) = hat_x(n) + hat_s(m+1,n) * x(n-m);
            end
        end
        
        %%%%%%%%%%%%%%%%
        % Update starting
        %%%%%%%%%%%%%%%%
        % Update of W_N(z)
        for j=0:L_N-1
            if n-soft_D-j > -0
               w_N(j+1, n+1) =  w_N(j+1, n) + mu_N * e(n) * hat_x(n-soft_D-j);
            end
        end
                              
    end %  n loop end
    
    %%%%%% Remaining noise saving
    save_e_A    = save_e_A + e(1:N)/T;           %e(n)
    save_e2_A  = save_e2_A + e(1:N).*e(1:N)/T;   %e(n)^2

    save_NRP_A_1st(t) = 10*log10( var(e(Nn-N_eva:Nn)) ...
                                   / ( var( p(Nn-N_eva:Nn) ) ) );
    save_NRP_A_2nd(t) = 10*log10( var(e(N-N_eva:N)) ...
                                   / ( var( p(N-N_eva:N) ) ) );
    
end   %T loop end
 
Time_sys_A = toc;

%
save Case_1_Sys_A.mat   save_e_A    save_e2_A  ...
                          save_NRP_A_1st  save_NRP_A_2nd;

%%% The end  
