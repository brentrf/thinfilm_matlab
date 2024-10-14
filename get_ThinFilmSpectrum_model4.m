%%%%%%%%%%%
%  "Coating Model 3"  
%
%  Purpose is to model the transmission spectrum for a multilayer coating 
%  for light at a given incidence angle...
%
%  This is for a multilayer coating that has 
%  ** arbitrary layers ** (each has different thickness and n) 
%  e.g  "ABCDEF..." structure
%
%   NEW(3):  Added S-polarization for incident beam (earlier versions
%   considered only P-polarized incident beam)
%
%   NEW(4):  Adding dispersion:  vary index with wavelength
%   NEW(4):  Adding losses:   imaginary refractive index (kappas)
%
%
%%%%%%%%%%

function [Rp,Rs,Tp,Ts,T_tot,T_np,Abs_tot,Abs_np, Mtest] = get_ThinFilmSpectrum_model4(theta0_deg, theta_pol_deg, T, mat_nk, n0, ns, lam)

    %INPUT PARAMETERS
    T;                              %Layer thicknesses  (nm) 
    mat_n       = mat_nk(:,:,1);  	%real refractive index (layer, lambda) 
    mat_kappas  = mat_nk(:,:,2);	%imag refractive index (layer, lambda) 
    %theta0_deg = 30;               %angle of incidence (degrees)
    %n0 = 1;                %index of refraction before first layer   (examples:  n(fused-SiO2)=1.46,  n()=   )
    %ns = 1.0;              %index of refraction after last layer   (examples:  n(TiO2-amorph)=2.3,  n(TiO2-anatase)=2.4, n(TiO2-rutile)=2.6 )
    %theta_pol_deg = 00;  %polarization angle (0 = p-polarization; 90=s-polarization)

    %CALCULATIONS
    eps0=8.854*10^-12;      %C/(V*m)  = 1/(4*pi*8.987551788e9) C^2/N-m2
    mu0=4*pi*10^-7;         %N/A^2
    Z0 = (eps0/mu0)^-0.5;   %Ohm
    dlam = median(lam(2:end)-lam(1:end-1));        
    P = length(T); %number of layers
    if size(mat_n,1)~=P    error('Number of index inputs does not match number of layers');    end
    theta0 = pi*theta0_deg/180;         %angle of incidence
    thetaS = asin(n0/ns*sin(theta0));   %angle of transmittance
       %assume *no loss* for incident and substrate 
    Y0 = n0*sqrt(eps0/mu0)*cos(theta0); %parameter for p-polarization  
    Ys = ns*sqrt(eps0/mu0)*cos(thetaS); %parameter for p-polarization
    X0 = n0*sqrt(eps0/mu0)/cos(theta0); %parameter for s-polarization
    Xs = ns*sqrt(eps0/mu0)/cos(thetaS); %parameter for s-polarization
    polvec = [cos(theta_pol_deg*pi/180),sin(theta_pol_deg*pi/180)];
       %[p-polarization;s-polarization] %norm pol vector to unit length ... %polvec=polvec/sqrt(dot(polvec,polvec));


    for k=1:length(lam)
        lambda = lam(k);
        nreal  = mat_n(:,k);       %vertical vector of (real) refract.index at this wavelength
        kappas = mat_kappas(:,k);  %vertical vector of (imag) refract.index at this wavelength
        n      = nreal + i*kappas; %vertical vector of complex refract.index  @ this wavelength
        
        M=eye(2); S=eye(2);  
        for kp=1:P  %loop through stack and generate product matrix, "M"
            %   theta_kp =  asin(n0/n(kp)*sin(theta0));          %angle of incidence at bottom of kp_th layer
            cos_theta_kp = sqrt(1 - (n0/n(kp)*sin(theta0))^2 );    %cos(angle of incidence at bottom of kp_th layer)
            phi_kp = -2*pi*n(kp)*T(kp)*cos_theta_kp/lambda;  %phase delay of one-pass through "kp_th" layer of material

            %this is the matrix for P-polarization
                Ykp = n(kp)*sqrt(eps0/mu0)*cos_theta_kp;  
                Mkp = [cos(phi_kp), i*sin(phi_kp)/Ykp; i*sin(phi_kp)*Ykp, cos(phi_kp)];        
                M = M*Mkp;

            %this is the matrix for S-polarization
                Xkp = n(kp)*sqrt(eps0/mu0)/cos_theta_kp;  
                Skp = [cos(phi_kp), i*sin(phi_kp)/Xkp; i*sin(phi_kp)*Xkp, cos(phi_kp)];
                S = S*Skp;
        end

        if k==1
            Mtest = M; %total layer stack M  @ lam(1)
        end    
        
        r_p(k) = ( Y0*M(1,1) + Y0*Ys*M(1,2) - M(2,1) - Ys*M(2,2) ) / ( Y0*M(1,1) + Y0*Ys*M(1,2) + M(2,1) + Ys*M(2,2) );
        t_p(k) = 2*Y0 / ( Y0*M(1,1) + Y0*Ys*M(1,2) + M(2,1) + Ys*M(2,2) );
            %the following are by analogy to the p-polarization equations...
        r_s(k) = ( X0*S(1,1) + X0*Xs*S(1,2) - S(2,1) - Xs*S(2,2) ) / ( X0*S(1,1) + X0*Xs*S(1,2) + S(2,1) + Xs*S(2,2) );
        t_s(k) = 2*X0 / ( X0*S(1,1) + X0*Xs*S(1,2) + S(2,1) + Xs*S(2,2) );
            %the following are by my derivation ... 
        %r_s(k) = ( -X0*S(1,1) - X0*Xs*S(1,2) + S(2,1) + Xs*S(2,2) ) / ( X0*S(1,1) + X0*Xs*S(1,2) + S(2,1) + Xs*S(2,2) );
        %t_s(k) = 2*Xs*(n0/ns) / ( X0*S(1,1) + X0*Xs*S(1,2) + S(2,1) + Xs*S(2,2) );

    end


    %calculate reflectivity/transmissivity by polarization component (s and p polarization)
    Rp = (abs(r_p).^2); Rs = (abs(r_s).^2); 
    Tp = (abs(t_p).^2); Ts = (abs(t_s).^2);  %note: abs(t_s).^2 = t_s.*conj(t_s)
    R_tot =  (polvec(1)^2)*Rp + (polvec(2)^2)*Rs;  %note Rs and Ts are SQUARED -- Watts_out/Watts_in
    T_tot =  (polvec(1)^2)*Tp + (polvec(2)^2)*Ts  ;  %note Rs and Ts are SQUARED -- Watts_out/Watts_in
    Abs_p = 1 - (Rp+Tp);  Abs_s = 1-(Rs+Ts);
    Abs_tot = 1 - (R_tot+T_tot);

    % Now calculate the result for a statistical distribution of polarizations
    % (non-polarized light)
    R_np =  0.5*Rp + 0.5*Rs;  %note Rs and Ts are SQUARED -- Watts_out/Watts_in
    T_np =  0.5*Tp + 0.5*Ts  ;  %note Rs and Ts are SQUARED -- Watts_out/Watts_in
    Abs_np = 1 - (R_np+T_np);                   
return;


%    %Plot BOTH S & P Transmission Spectrum
%    figure; plot(lam,[Tp',Ts',T_tot',T_np']); xlabel('wavelength [nm]'); ylabel('T := |Et/Ei|^2'); grid on; xlim([min(lam),1700]);
%    legend('P-polarized','S-polarized',strcat('Total(pol=',num2str(theta_pol_deg),'deg)'),'Total(NONpolarized Light)'); 
%    title({strcat('angle of incidence = ',num2str(theta0_deg),' degrees'); ...
%        strcat('polarization angle =  ',num2str(theta_pol_deg),' degrees
%        (p-pol:=0deg):')}); 

%    %Plot BOTH S & P Reflectivity Spectrum
%    figure; plot(lam,[Rp',Rs',R_tot',R_np']); xlabel('wavelength [nm]'); ylabel('T := |Et/Ei|^2'); grid on; xlim([min(lam),1700]);
%    legend('P-polarized','S-polarized',strcat('Total(pol=',num2str(theta_pol_deg),'deg)'),'Total(NONpolarized Light)'); 
%    title({strcat('angle of incidence = ',num2str(theta0_deg),' degrees'); ...
%        strcat('polarization angle =  ',num2str(theta_pol_deg),' degrees (p-pol:=0deg):')}); 

%    %for diagnostics... 
%    figure; plot3(real(t_s),imag(t_s),lam); xlim([-1,1]);ylim([-1,1]);
%    figure; plot(Rs+Ts); figure; plot(Rp+Tp);

