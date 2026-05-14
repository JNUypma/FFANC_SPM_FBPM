%%%%
% %%% FFNANC system by Bai et al. 2019 
%%%%
clear all;
close all;

N     = 70000;
Nn    = N/2;
N_eva = 7000;
T     = 100;

% narrowband reference signal
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
M_change = 1;
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
soft_D = 1;   % Delay

%%%
% Source signal
x_clean = zeros(1, N);
for i = 1:q
    x_clean = x_clean + Amp(i) * sin(omega(i)* (0:N-1) + 0);
end

%%%%
% additive noise in the narrowband source
std_v_s = sqrt(0.001);  
load v_s_org.mat;

%%%
% v_p(n)
std_v_p = sqrt(0.01);  % additive noise std dev on p(n)
load v_p_org.mat; 

%%%
% AWGN
load v_o_org.mat; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Each paramaters of Sys_C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% feedback path estimation
hat_M_f = 42;
mu_f      = 0.001; %0.001;

%Secondary path estimation
hat_M_s = 31;
mu_s     = 0.001; %0.001;

% MPA Controller
mu = 0.0025* ones(1, q);  

% Notch bandpass filters : NG algorithm
rho   = 0.95 * ones(1, q);  % pole attraction para %0.925
epsil = 0.01;
betaa = 0.999;
mu_nt = 0.001 * ones(1, q);                   %0.001

% AWGN scaling factor
std_v_o  = sqrt(1.0);
alphh  = 0.999; 
betaa1 = 1-alphh;
gamma = 2;

%%%%%%%%
%%% results
%%%%%%%%%
save_e_C = zeros(1, N);
save_e2_C = zeros(1, N);

save_c_C           = zeros(q, N);
save_y_mpa_C  = zeros(q, N);
save_G_C  = zeros(1, N);

save_MSE_FBPM_C   = zeros(1, N);
save_MSE_OSPM_C   = zeros(1, N);

save_NRP_C_1st = zeros(1,T);
save_NRP_C_2nd = zeros(1,T);

tic;
for t = 1:T
    
    disp(' t =====> '); disp(t);
      
    % 3 additive noises
    v_o = v_o_org(t,1:N) * std_v_o;      % Auxiliary noise
    v_p = v_p_org(t,1:N) * std_v_p;      % additive noise in p(n)
    v_s = v_s_org(t,1:N) * std_v_s;      % additive noise in noise source x_s(n)

    % Source signal with additive noise 
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
    
    % primary noise + v_p(n)
    p   = p_clean + v_p;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% NANC with FBPM and SPM
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    %reference
    x_r     = zeros(1,N);
    hat_x_r = zeros(1,N);
    hat_x_s = zeros(q,N);  % output of BP (hat_x_r)
    hat_x_ss = zeros(q,N);  %filted-x hat_x_s
       
    %%%%%%%%%%%%%%%%
    % paramaters of notch filter
    %%%%%%%%%%%%%%%%
    c = zeros(q, N);    %notch filtere paraAmp
    g = zeros(q, N);
    G = zeros(q, N);        % lowpassed gradient power
    u_N  = zeros(q, N); %notch filter output of reference hat_x_r
    
    %%%%%%%%%%%%%%%%%
    % paramaters of MPA
    %%%%%%%%%%%%%%%%%
    y_0 = zeros(1,N); %output of MPA(sum)
    y_m =zeros(q, N); %output of MPA(one channal)
    h0 = zeros(q, N); %paramaters of MPA
    h1 = zeros(q, N);
 
    %%%%%%%%%%%%%%%%%
    % paramaters of FBPM
    %%%%%%%%%%%%%%%%%
    y_f       = zeros(1, N);
    hat_y_f   = zeros(1, N);
    hat_f  = zeros(hat_M_f, N);

    %%%%%%%%%%%%%%%%%
    % paramaters of OSPM
    %%%%%%%%%%%%%%%%%
    hat_s  = zeros(hat_M_s,N);
    y_s    = zeros(1,N);
    e_s    = zeros(1,N); 
      
    %%%%%%%%%%%%%%%%%
    % paramaters of gain
    %%%%%%%%%%%%%%%%%
    Gain = zeros(1, N); % scaling factor

    % Auxiliary noise
    v  = zeros(1, N);
    
    y     = zeros(1,N); % 2nd source
    y_p   = zeros(1,N); % secondary path after SP
    e     = zeros(1,N); % residual noise
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% NANC
     
    for n=2:N  % time loop start

        % auxiliary noise scaled
        if n-1>0
            Gain(n) = alphh * Gain(n-1) + betaa1 *(abs(e(n-1))^gamma);
        end
        v(n) = v_o(n) *  Gain(n);

        % y(n)
        if n-1>0
            y(n) = y_0(n-soft_D) + v(n); %%% plus
        end
                
        % y_f : true FBP output
        for j = 0:M_fb-1
            if n-j > 0
                y_f(n) = y_f(n) + s_fb(j+1,n) * y(n-j);
            end
        end
        
        % x_r(n) : reference signal
        x_r(n) = x_s(n) + y_f(n);
        
        % hat_y_f : estimate of FBP signal
        for i=0:hat_M_f-1
            if n-i>0
                hat_y_f(n) = hat_y_f(n) + hat_f(i+1, n) * y(n-i);
            end
        end
        
        %hat_x_r(n) : input of Notch filter
        hat_x_r(n) = x_r(n)- hat_y_f(n);
        
        %%%%%%%%%%%%%%%%%%%%
        % Notch filter
        %%%%%%%%%%%%%%%%%%%%
        if n-2 >0
            for i=1:q
                if i > 1 %(i > 1)
                    u_N(i, n) = -rho(i)*c(i, n)*u_N(i, n-1)-rho(i)^2*u_N(i, n-2) ...
                                +u_N(i-1, n) + c(i, n) * u_N(i-1, n-1) + u_N(i-1, n-2);

                    % hat_x_s : BP output (contain the sinusoid)
                    hat_x_s(i, n) = u_N(i-1,n) - u_N(i, n);
                else %(i = 1)
                    u_N(i, n) = -rho(i)*c(i, n)*u_N(i, n-1)-rho(i)^2*u_N(i, n-2) ...
                                  + hat_x_r(n) + c(i, n) * hat_x_r(n-1) + hat_x_r(n-2);

                    % hat_x_s
                    hat_x_s(i, n) = hat_x_r(n) - u_N(i,n);
                end
            end
        end

        % ANFB update: normalized gradient algorithm (NG)
        if n-1>0
            for i = 1:q
                if i>1
                    g(i, n) =  u_N(i-1, n-1) - rho(i) * u_N(i, n-1);
                else
                    g(i, n) =  hat_x_r(n-1) - rho(i) * u_N(i, n-1);
                end
                G(i, n) = alphh * G(i, n-1) + (1-alphh) * g(i, n)^2;
                c(i, n+1) = c(i, n) - mu_nt(i) * u_N(i, n)* g(i, n) / (epsil + G(i, n));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%
        %   Output of MPA
        %%%%%%%%%%%%%%%%%%%
        if n-1>0
            for i = 1:q
                % MPA outputs (one channal)
                y_m(i, n) = h0(i, n) * hat_x_s(i, n) + h1(i, n) * hat_x_s(i, n-1);
            end
        end

        % MPA output (sum)
        y_0(n) = sum(y_m(:, n));

        %y_p : signal output of 2nd source passing 2nd path
        for m=0:M_s-1
            if n-m>0
                y_p(n) = y_p(n) + s(m+1,n) * y(n-m);
            end
        end
       
        %e(n) : residual noise
        e(n) = p(n) - y_p(n);    

        %OSPM
        %y_s: out put of hat_s
        for m=0:hat_M_s-1
            if n-m>0
                y_s(n) = y_s(n) + hat_s(m+1,n) * v(n-m);
            end
        end

        e_s(n) = e(n) + y_s(n); %%%% plus
               
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% Update preparing
        %%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%
        % preparing for MPA update
        %%%%%%%%%%%%%
        % hat_x_ss(i, n) : filted-x of hat_x_s
        for i=1:q
            for m=0:hat_M_s-1
                if n-m>0
                    hat_x_ss(i, n) = hat_x_ss(i, n) + hat_s(m+1,n) * hat_x_s(i, n-m);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%
        % Update starting
        %%%%%%%%%%%%%%%%
        %%% update of  h0 , h1
        if n-1- soft_D>0
            for i = 1:q
                h0(i, n+1) = h0(i, n) + mu(i) * e(n) * hat_x_ss(i, n - soft_D);
                h1(i, n+1) = h1(i, n) + mu(i) * e(n) * hat_x_ss(i, n-1 - soft_D);
            end
        end
              
        % updating the FBPM
        for i=0:hat_M_f-1
            if n-i>0
                hat_f(i+1, n+1) =  hat_f(i+1, n) + mu_f * u_N(q, n) * v(n-i);
            end
        end
        
        % updating the SPM
        for i=0:hat_M_s-1
            if n-i > 0
                hat_s(i+1, n+1) = hat_s(i+1, n) - mu_s * e_s(n) * v(n-i); %%minus
            end
        end
                
    end %  n loop end
    
    %%%%%% Remaining noise saving
    save_e_C   = save_e_C + e(1:N)/T;           %e(n)
    save_e2_C  = save_e2_C + e(1:N).*e(1:N)/T;   %e(n)^2

    save_G_C   = save_G_C +  Gain(1:N)/T;

    for i=1:q
        save_c_C(i, 1:N)     = save_c_C(i, 1:N) + c(i, 1:N)/T;
        save_y_mpa_C(i, 1:N) = save_y_mpa_C(i, 1:N) + y_m(i, 1:N)/T;
    end

    save_NRP_C_1st(t) = 10*log10( var(e(Nn-N_eva:Nn)) ...
                                   / ( var( p(Nn-N_eva:Nn) ) ) );
    save_NRP_C_2nd(t) = 10*log10( var(e(N-N_eva:N)) ...
                                   / ( var( p(N-N_eva:N) ) ) );


    %%%%%% MSE of FBP estimation
    for n = 1:N
        if hat_M_f > M_fb
            wk_err_fp(1:M_fb)       = hat_f(1:M_fb, n)' - s_fb(1:M_fb,n)';
            wk_err_fp(M_fb+1:hat_M_f) = hat_f(M_fb+1:hat_M_f, n)' - 0;
        else
            wk_err_fp(1:hat_M_f)   = hat_f(1:hat_M_f, n)' - s_fb(1:hat_M_f,n)';
            wk_err_fp(hat_M_f+1:M_fb) = 0 - f(hat_M_f+1:M_fb,n);
        end

        save_MSE_FBPM_C(n)    = save_MSE_FBPM_C(n) + ( wk_err_fp * wk_err_fp') / T;      
    end

    %%%%%% MSE of SPM estimation
    for n = 1:N
        if hat_M_s > M_s
            wk_err_sp(1:M_s)         = hat_s(1:M_s, n)' - s(1:M_s,n)';
            wk_err_sp(M_s+1:hat_M_s) = hat_s(M_s+1:hat_M_s, n)' - 0;
        else
            wk_err_sp(1:hat_M_s)   = hat_s(1:hat_M_s, n)' - s(1:hat_M_s,n)';
            wk_err_sp(hat_M_s+1:M_s) = 0 - s(hat_M_s+1:M_s,n)';
        end

        save_MSE_OSPM_C(n)    = save_MSE_OSPM_C(n) + ( wk_err_sp * wk_err_sp') / T;
    end

end   %T loop end

Time_sys_C = toc;

save_MSE_FBPM_C(1:Nn) = save_MSE_FBPM_C(1:Nn)/ (s_fb(:,1)' * s_fb(:,1) );
save_MSE_FBPM_C(Nn+1:N) = save_MSE_FBPM_C(Nn+1:N)/ (s_fb(:,Nn+1)' * s_fb(:,Nn+1) );

save_MSE_OSPM_C(1:Nn) = save_MSE_OSPM_C(1:Nn)/ (s(:,1)' * s(:,1) );
save_MSE_OSPM_C(Nn+1:N) = save_MSE_OSPM_C(Nn+1:N)/ (s(:,Nn+1)' * s(:,Nn+1) );


%saving
save Case_1_Sys_C.mat  ...
                              save_e_C          save_e2_C  ...
                              save_G_C          save_NRP_C_1st    save_NRP_C_2nd  ...
                              save_MSE_FBPM_C   save_MSE_OSPM_C;

%%%%%%% the end  %%%%%%%%


