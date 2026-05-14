%%%%% New system with OSPM and OFBPM
close all
clear all

N     = 70000;
Nn    = N/2;
N_eva = 7000;
T     = 100;

% narrowband reference []
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

%%%%%%%%%%%%%%%
soft_D = 1;    % Delay

% Source signal
x_clean = zeros(1, N);
for i = 1:q
    x_clean = x_clean + Amp(i) * sin(omega(i) * (0:N-1) + 0);
end

% Additive noise in the reference (narrowband source)
std_v_s = sqrt(0.001); 
load v_s_org.mat;

% Additive noise in p(n)
% v_p(n)
std_v_p = sqrt(0.01); 
load v_p_org.mat;

% AWGN for OFBPM and OSPM
load v_o_org.mat;

% OFBPM
hat_M_f = 42;
mu_f    = 0.0025;   

% OSPM, AWGN scalingfd
std_v_o  = sqrt(1.0);
hat_M_s = 31;
mu_s   = 0.003; 
alphh  = 0.9996;
betaa  = 0.0018;
gamma  = 2;

% W_N(z)   ---- FIR type
%%%%%%%%%%%%%%%%%%%
L_N = M_p - M_s +5;
mu_N = 0.0003;  

% Q_1(z): LPF 1
L_sf   = hat_M_f+5;        
mu_LOC = 0.003;    

% Q_2(z): LPF 2
L_H  = L_N+5;          
mu_GLO = 10 * mu_N;  

%%%%%%%%%
%%% results
%%%%%%%%%
save_p          = zeros(1, N);
save_e_D          = zeros(1, N);
save_e2_D         = zeros(1, N);
save_MSE_FBPM_D   = zeros(1, N);
save_MSE_OSPM_D   = zeros(1, N);
save_G_D          = zeros(1, N);
save_NRP_D_1st = zeros(1,T);
save_NRP_D_2nd = zeros(1,T);

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
    x_r      = zeros(1, N);
    hat_x_r  = zeros(1, N);
    hat_hat_x_r = zeros(1, N);

    y_0 = zeros(1, N);  % output of W_N(z)
    w_N = zeros(L_N, N); % paramaters of W_N(z)
 
    % OFBPM
    y_f      = zeros(1, N);
    hat_y_f  = zeros(1, N);
    hat_f    = zeros(hat_M_f, N);

    % OSPM, AWGN
    hat_s   = zeros(hat_M_s, N);
    y_s     = zeros(1,N);
    e_s     = zeros(1,N); 
    v       = zeros(1, N);
    
    % SF1
    h_loc = zeros(L_sf, N);
    y_loc = zeros(1,N);
    e_loc = zeros(1,N);            
    
    % Global SF (GSF)
    h_glo = zeros(L_H, N);
    e_glo = zeros(1, N);               
    y_glo = zeros(1, N);         
    
    % paramaters of gain
    %%%%%%%%%%%%%%%%%
    Gain = zeros(1, N); % scaling factor
    
    y     = zeros(1,N); %2nd source
    y_p   = zeros(1,N);%secondary path after SP
    e     = zeros(1,N); %residual noise
                
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%% NANC
     
    for n=3:N  % time loop start
        
        % auxiliary noise scaled
        Gain(n) = alphh * Gain(n-1)  + betaa*abs(y_glo(n-1))^gamma;
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
        
        % x_r(n) : reference signal
        x_r(n) = x_s(n) + y_f(n);
        
        % hat_y_f : estimate of FBP signal
        hat_y_f(n) = 0;
        for m=0:hat_M_f-1
            if n-m>0
                hat_y_f(n) = hat_y_f(n) + hat_f(m+1, n) * y(n-m);
            end
        end
        
        % hat_x_r(n) : input of Notch filter
        hat_x_r(n) = x_r(n)- hat_y_f(n);
        
        %  Local SF      
        for j=0:L_sf-1
            if n-soft_D-j > 0
                y_loc(n) = y_loc(n) + h_loc(j+1, n) * y_0(n-soft_D-j);
            end
        end
        
        % LSF error
        e_loc(n) = hat_x_r(n) - y_loc(n);   
        
        % LSF update 
        for j=0:L_sf-1
            if n-soft_D-j > 0
                h_loc(j+1, n+1)  =  h_loc(j+1, n) ...
                                  + mu_LOC * e_loc(n) * y_0(n-soft_D-j);
            end
        end
        
        %  Output of W_N(z)
        for j=0:L_N-1
            if n-j > 0
                y_0(n) = y_0(n) + w_N(j+1, n) * hat_x_r(n-j);
            end
        end
                     
        % y_p : signal output of 2nd source passing the SP
        y_p(n) = 0;
        for m=0:M_s-1
            if n-m>0
                y_p(n) = y_p(n) + s(m+1,n) * y(n-m);
            end
        end
       
        %e(n) : residual noise
        e(n) = p(n) - y_p(n);      
        
        % GSF output
        for j=0:L_H-1
            if n-j > 0
                y_glo(n) = y_glo(n) + h_glo(j+1, n) * hat_x_r(n-j);   % remaining narrowband component
            end
        end
        
        % GSF error
        e_glo(n) = e(n) - y_glo(n);    % GSF error: remaining broadband component
                       
        % OSPM output, error
        for m=0:hat_M_s-1
             if n-m>0
                   y_s(n) = y_s(n) + hat_s(m+1,n) * v(n-m);
             end
        end
        
        e_s(n) = e_glo(n) + y_s(n);

        % Preparations for update
        for m=0:hat_M_s-1
            if n-m>0
                hat_hat_x_r(n) = hat_hat_x_r(n) + hat_s(m+1,n) * hat_x_r(n-m);
            end
        end

        % Update of W_N(z)
        for j=0:L_N-1
            if n-soft_D-j > -0
               w_N(j+1, n+1) =  w_N(j+1, n) ...
                   + mu_N * y_glo(n) * hat_hat_x_r(n-soft_D-j);
            end
        end
        
        % Update of hat{F}(z)
        for j=0:hat_M_f-1
            if n-j>0
                hat_f(j+1, n+1) = ...
                    hat_f(j+1, n) + mu_f * e_loc(n) * v(n-j);
            end
        end
        
        % Update of \hat{S}(z)
        for m=0:hat_M_s-1
            if n-m > 0
                hat_s(m+1, n+1) = hat_s(m+1, n) - mu_s * e_s(n) * v(n-m);
            end
        end
        
        % G_SF update
        for j=0:L_H-1
            if n-j > 0
                h_glo(j+1, n+1)  =  h_glo(j+1, n) + mu_GLO * e_s(n) * hat_x_r(n-j);
            end
        end
                                  
    end %  n loop end
    
    save_p   = save_p + p(1:N)/T;           %e(n)
    
    %%%%%% Remaining noise saving
    save_e_D   = save_e_D + e(1:N)/T;           %e(n)
    save_e2_D  =  save_e2_D  + e(1:N).*e(1:N)/T;   %e(n)^2
    
    save_G_D   = save_G_D +  Gain(1:N)/T;

    save_NRP_D_1st(t) = 10*log10( var(e(Nn-N_eva:Nn)) ...
                                   / ( var( p(Nn-N_eva:Nn) ) ) );
    save_NRP_D_2nd(t) = 10*log10( var(e(N-N_eva:N)) ...
                                   / ( var( p(N-N_eva:N) ) ) );

    %%%%%% MSE of FBPM estimation
    for n=1:N
        if hat_M_f > M_fb
            wk_err_fp(1:M_fb)       = hat_f(1:M_fb, n)' - s_fb(1:M_fb,n)';
            wk_err_fp(M_fb+1:hat_M_f) = hat_f(M_fb+1:hat_M_f, n)' - 0;
        else
            wk_err_fp(1:hat_M_f)   = hat_f(1:hat_M_f, n)' - s_fb(1:M_fb,n)';
            wk_err_fp(hat_M_f+1:M_fb) = 0 - f(hat_M_f+1:M_fb,n)';
        end

        save_MSE_FBPM_D(n)    = save_MSE_FBPM_D(n) + ( wk_err_fp * wk_err_fp') / T;
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

        save_MSE_OSPM_D(n)    = save_MSE_OSPM_D(n) + ( wk_err_sp * wk_err_sp') / T;
    end   
    
end   %T loop end

Time_sys_D = toc;

save_MSE_FBPM_D(1:Nn) = save_MSE_FBPM_D(1:Nn)/ (s_fb(:,1)' * s_fb(:,1) );
save_MSE_FBPM_D(Nn+1:N) = save_MSE_FBPM_D(Nn+1:N)/ (s_fb(:,Nn+1)' * s_fb(:,Nn+1) );
    
save_MSE_OSPM_D(1:Nn) = save_MSE_OSPM_D(1:Nn)/ (s(:,1)' * s(:,1) );
save_MSE_OSPM_D(Nn+1:N) = save_MSE_OSPM_D(Nn+1:N)/ (s(:,Nn+1)' * s(:,Nn+1) );
    
save  Case_1_Sys_D.mat    save_e_D          save_e2_D  ...
                              save_G_D          save_NRP_D_1st  save_NRP_D_2nd  ...
                              save_MSE_FBPM_D   save_MSE_OSPM_D;

load  Case_1_Sys_A.mat;
load  Case_1_Sys_B.mat;
load  Case_1_Sys_C.mat;

N1  =  20;
set(0, 'defaultAxesFontName', 'Arial' );
set(0, 'defaultAxesFontSize', 16);

figure(1);
plot(1:N1:N,10*log10(save_e2_A(1:N1:N)),'-k', ...
        1:N1:N,10*log10(save_e2_B(1:N1:N)),'-b', ...
        1:N1:N,10*log10(save_e2_C(1:N1:N)),'-m',  ...
        1:N1:N,10*log10(save_e2_D(1:N1:N)),'-r', 'LineWidth',1);
axis( [1    N     -25   5] );
legend('Sys-A','Sys-B','Sys-C','Sys-D');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
   'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('E[e^2(n)] [dB]', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

figure(2);
plot(1:N,save_G_B(1:N),'-b', ...
        1:N,save_G_C(1:N),'--m', ...
        1:N,save_G_D(1:N),'.-r');                
legend('Sys-B','Sys-C','Sys-D');
axis( [1    N     -0.2   1.6] );
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
   'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('G_s(n)', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

figure(3);
subplot(2,1,1);
plot(1:N1:N,10*log10(save_MSE_FBPM_B(1:N1:N)),'-b', ...
        1:N1:N,10*log10(save_MSE_FBPM_C(1:N1:N)),'--m', ...
        1:N1:N,10*log10(save_MSE_FBPM_D(1:N1:N)),'.-r','LineWidth',1);
axis( [1    N     -40   10] );
legend('Sys-B','Sys-C','Sys-D');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
   'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('J_f(n) [dB]', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

subplot(2,1,2);
plot(1:N1:N,10*log10(save_MSE_OSPM_B(1:N1:N)),'-b',...
        1:N1:N,10*log10(save_MSE_OSPM_C(1:N1:N)),'--m', ...
        1:N1:N,10*log10(save_MSE_OSPM_D(1:N1:N)),'.-r',  'LineWidth',1);
axis( [1    N     -40   10] );
legend('Sys-B','Sys-C','Sys-D');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
   'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('J_s(n) [dB]', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

figure(4)
plot(1:T,  save_NRP_A_2nd(1:T),'*-k', ...
        1:T,  save_NRP_B_2nd(1:T),'x-b', ...
        1:T,  save_NRP_C_2nd(1:T),'o-m', ...
        1:T,  save_NRP_D_2nd(1:T),'s-r', 'LineWidth',1);
legend('Sys-A','Sys-B','Sys-C','Sys-D');
axis( [1    T     -24   -2] );
title('2nd half')
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Independent runs', ...
   'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('NRP [dB]', ...
    'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

disp('1st half NRPs ( A,  B, C,  D)  ==>');
disp( [ mean(save_NRP_A_1st(1:T))  mean(save_NRP_B_1st(1:T))  ...
            mean(save_NRP_C_1st(1:T))     mean(save_NRP_D_1st(1:T))  ] );
         
disp('2nd half NRPs ( A,  B, C,  D)  ==>');
disp( [ mean(save_NRP_A_2nd(1:T))  mean(save_NRP_B_2nd(1:T))  ...
            mean(save_NRP_C_2nd(1:T))     mean(save_NRP_D_2nd(1:T)) ] );

%%%%%%%%% THE END  %%%%%%%%%%%%%%%



