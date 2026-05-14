%%%%
%%%% FFNANC system with FBPM and SPM by Akhtar etal.
%%%%
close all
clear all

N     = 70000;
Nn    = N/2;
N_eva = 7000;
T     = 100;

% narrowband reference []
Fs = 2000; %sampling fre
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

%%%%%%%%%%%%%%%
soft_D = 1; %Delay
%%%%%%%%%%%%%%

% Source signal
x_clean = zeros(1, N);
for i = 1:q
    x_clean = x_clean + Amp(i) * sin(omega(i)* (0:N-1) + 0);
end

% Additive noise in the reference (narrowband source)
std_v_s = sqrt(0.001);  
load v_s_org.mat;

% Additive noise in p(n)
% v_p(n)
std_v_p = sqrt(0.01); 
load v_p_org.mat;

% AWGN
FLG_awgn = 1;
if FLG_awgn == 0
    v_o_org = randn(100, 120000);
    save v_o_org.mat   v_o_org;
else
    load v_o_org.mat; 
end

% OFBPM
hat_M_f = 42;
mu_f    = 0.002;    

% OSPM, AWGN scaling
std_v_o  = sqrt(1.0);
hat_M_s  = 31;
mu_s    = 0.002;   
alphh    = 0.9995;
gamma = 0.0005;
lambda = 0.998;

% W_N(z)   ---- FIR type
%%%%%%%%%%%%%%%%%%%
L_N  = M_p - M_s + 5;
mu_N = 0.0003; 

% SF: supporting filter
L_sf   = hat_M_f + 5;
mu_LOC  = 0.001;  

%%%%%%%%%
%%% results
%%%%%%%%%
save_e_B          = zeros(1, N);
save_e2_B         = zeros(1, N);
save_MSE_FBPM_B   = zeros(1, N);
save_MSE_OSPM_B   = zeros(1, N);
save_G_B          = zeros(1, N);
save_NRP_B_1st    = zeros(1,T);
save_NRP_B_2nd    = zeros(1,T);

tic;
for t = 1:T
    
    disp(' t =====> '); disp(t);
      
    % 3 additive noises
    v_o = v_o_org(t,1:N) * std_v_o;      % Auxiliary noise
    v_p = v_p_org(t,1:N) * std_v_p;             % additive noise in p(n)
    v_s = v_s_org(t,1:N) * std_v_s;           % additive noise in noise source x_s(n)

    % Source signal
    x_s = x_clean + v_s;
    
    % primary noise
    p_clean = zeros(1, N);
    for n=1:N
        for i=0:M_p-1
            if n-i>0
                p_clean(n) = p_clean(n) + s_p(i+1)*x_s(n-i);
            end
        end
    end
    
    p   = p_clean + v_p;
    
    %%%%%%%%%%%%%%%%%%%%   
    %reference
    c  = zeros(1, N);
    x  = zeros(1, N);
    hat_x  = zeros(1, N);

    y_0 = zeros(1,N);  % output of W_N(z)
    w_N = zeros(L_N, N); % paramaters of W_N(z)
 
    % OFBPM
    y_f       = zeros(1, N);
    hat_y_f   = zeros(1, N);
    hat_f  = zeros(hat_M_f, N);

    % OSPM, AWGN
    hat_s  = zeros(hat_M_s,N);
    y_s     = zeros(1,N);
    e_s     = zeros(1,N); 
    v         = zeros(1, N);
    
    %%%% for supporting filter
    w_SF = zeros(L_sf, N);
    y_h  = zeros(1,N);
    e_h  = zeros(1,N);    
    
    % paramaters of gain
    %%%%%%%%%%%%%%%%%
    Gain = zeros(1, N); % scaling factor
    Peh  = zeros(1, N);
    Pes  = zeros(1, N);
    
    y      = zeros(1,N); % 2nd source
    y_p    = zeros(1,N);   % secondary path after SP
    e      = zeros(1,N); % residual noise
   
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% hat_s & hat_f initialization
    hat_s(1:M_s, 3)  = s(1:M_s, 1) * ( 1 - 10^(-5/20) );
    hat_f(1:M_fb, 3) = s_fb(1:M_fb, 1) * ( 1 - 10^(-5/20) );
    
    for n=3:N  % time loop start

       % auxiliary noise scaled
        Gain(n) = alphh * Gain(n-1) ...
                 + gamma * max( sqrt(Peh(n-1)/(hat_f(1:hat_M_f, n)'*hat_f(1:hat_M_f, n))), ...
                                    sqrt(Pes(n-1)/ (hat_s(1:hat_M_s, n)'*hat_s(1:hat_M_s, n)))  );
        v(n) = v_o(n) *  Gain(n);       
        
        % y(n)
        y(n) = y_0(n-soft_D) + v(n);         
                
        % y_f : true FBP output
        y_f(n) = 0;
        for j = 0:M_fb-1
            if n-j > 0
                y_f(n) = y_f(n) + s_fb(j+1, n) * y(n-j);
            end
        end
        
        % c(n) : reference signal
        c(n) = x_s(n) + y_f(n);
        
        % hat_y_f : estimate of FBP signal
        hat_y_f(n) = 0;
        for m=0:hat_M_f-1
            if n-m>0
                hat_y_f(n) = hat_y_f(n) + hat_f(m+1, n) * y(n-m);
            end
        end
        
        % x(n)  
        x(n) = c(n)- hat_y_f(n);
        
        % supporting filter
        for m=0:L_sf-1
            if n-1-m > 0
               y_h(n) = y_h(n) + w_SF(m+1, n) * y_0(n-1-m);
            end
        end
        
        e_h(n) = x(n) - y_h(n);
                       
        % supporting filter update
        for j=0:L_sf-1
            if n-1-j > 0
                w_SF(j+1, n+1) = w_SF(j+1, n) + mu_LOC * e_h(n) * y_0(n-1-j);
            end
        end
        
        %  Output of W_N(z)
        for j=0:L_N-1
            if n-j > 0
                y_0(n) = y_0(n) + w_N(j+1, n) * x(n-j);
            end
        end
                             
        % y_p : signal output of 2nd source passing the SP
        y_p(n) = 0;
        for m=0:M_s-1
            if n-m>0
                y_p(n) = y_p(n) + s(m+1, n) * y(n-m);
            end
        end
       
        %e(n) : residual noise
        e(n) = p(n) - y_p(n);      
                       
        % OSPM output
        for m=0:hat_M_s-1
             if n-m>0
                   y_s(n) = y_s(n) + hat_s(m+1,n) * v(n-m);
             end
        end
        
        e_s(n) = e(n) + y_s(n);

        % Preparations for update
        for m=0:hat_M_s-1
            if n-m>0
                hat_x(n) = hat_x(n) + hat_s(m+1,n) * x(n-m);
            end
        end

        % Update of W_N(z)
        for j=0:L_N-1
            if n-soft_D-j > 0
                w_N(j+1, n+1) =  w_N(j+1, n) + mu_N * e_s(n) * hat_x(n-soft_D-j);
            end
        end
        
        % Update of hat{F}(z)
        for j=0:hat_M_f-1
            if n-j>0
                hat_f(j+1, n+1) = hat_f(j+1, n) + mu_f * e_h(n) * v(n-j);
            end
        end
        
        % Update of \hat{S}(z)
        for m=0:hat_M_s-1
            if n-m > 0
                hat_s(m+1, n+1) = hat_s(m+1, n) - mu_s * e_s(n) * v(n-m);
            end
        end
        
        Peh(n) = lambda * Peh(n-1) + (1-lambda)*e_h(n)^2;
        Pes(n) = lambda * Pes(n-1) + (1-lambda)*e_s(n)^2;
         
    end %  n loop end
    
    %%%%%% Remaining noise saving
    save_e_B   = save_e_B + e(1:N)/T;           %e(n)
    save_e2_B  = save_e2_B  + e(1:N).*e(1:N)/T;   %e(n)^2
    
    save_G_B   = save_G_B + Gain(1:N)/T;

    save_NRP_B_1st(t) = 10*log10( var(e(Nn-N_eva:Nn)) ...
                                   / ( var( p(Nn-N_eva:Nn) ) ) );
    save_NRP_B_2nd(t) = 10*log10( var(e(N-N_eva:N)) ...
                                   / ( var( p(N-N_eva:N) ) ) );

    %%%%%% MSE of FBPM estimation
    for n=1:N
        if hat_M_f > M_fb
            wk_err_fp(1:M_fb)         = hat_f(1:M_fb, n)' - s_fb(1:M_fb,n)';
            wk_err_fp(M_fb+1:hat_M_f) = hat_f(M_fb+1:hat_M_f, n)' - 0;
        else
            wk_err_fp(1:hat_M_f)   = hat_f(1:hat_M_f, n)' - s_fb(1:hat_M_f,n)';
            wk_err_fp(hat_M_f+1:M_fb) = 0 - f(hat_M_f+1:M_fb,n)';
        end

        save_MSE_FBPM_B(n)    = save_MSE_FBPM_B(n) + ( wk_err_fp * wk_err_fp') / T;
    end
    
    %%%%%% MSE of OSPM estimation
    for n=1:N
        if hat_M_s > M_s
            wk_err_sp(1:M_s)       = hat_s(1:M_s, n)' - s(1:M_s,n)';
            wk_err_sp(M_s+1:hat_M_s) = hat_s(M_s+1:hat_M_s, n)' - 0;
        else
            wk_err_sp(1:hat_M_s)   = hat_s(1:hat_M_s, n)' - s(1:hat_M_s,n)';
            wk_err_sp(hat_M_s+1:M_s) = 0 - s(hat_M_s+1:M_s,n)';
        end

        save_MSE_OSPM_B(n)    = save_MSE_OSPM_B(n) + ( wk_err_sp * wk_err_sp') / T;
    end    
       
end   %T loop end

Time_sys_B = toc;

save_MSE_FBPM_B(1:Nn) = save_MSE_FBPM_B(1:Nn)/ (s_fb(:,1)' * s_fb(:,1) );
save_MSE_FBPM_B(Nn+1:N) = save_MSE_FBPM_B(Nn+1:N)/ (s_fb(:,Nn+1)' * s_fb(:,Nn+1) );
    
save_MSE_OSPM_B(1:Nn) = save_MSE_OSPM_B(1:Nn)/ (s(:,1)' * s(:,1) );
save_MSE_OSPM_B(Nn+1:N) = save_MSE_OSPM_B(Nn+1:N)/ (s(:,Nn+1)' * s(:,Nn+1) );

save Case_1_Sys_B.mat    save_e_B          save_e2_B  ...
                              save_G_B          save_NRP_B_1st  save_NRP_B_2nd  ...
                              save_MSE_FBPM_B   save_MSE_OSPM_B;
                          
%%%%%%% THE END  %%%%%%%%
