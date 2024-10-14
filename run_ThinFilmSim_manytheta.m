%%%%%%%%%%%
%  "run "Coating Model 3" for many input angles  
%%%%%%%%%%

function run_coating_model_3_manytheta()

    clear; close all;
    fdir = '/Users/brentfisher/Documents/MATLAB/XDC/LightTrapMetal_ThinFilm/ThinFilmSim_Models';
    ddir = '/Users/brentfisher/Documents/MATLAB/XDC/LightTrapMetal_ThinFilm/ThinFilmSim_Models';
    ddir_materials = [ddir,'/refractive_index_data'];
    cd (ddir)
    
%% set up    
 
    %INPUT PARAMETERS:   Air & Substrate
    theta0_deg_vec      = [0:10:85];
%     theta0_deg          = 0;     %degrees
    theta_pol_deg       = 0;  %polarization angle (0 = p-polarization; 90=s-polarization)
    dlam                = 1.0; 
    lam                 = [300:dlam:1000]; %nm ... define wavelength space.
    N_upsample          = 20;  %# samples to increase & smooth  (used only if thick layer) 
    
   
    %MATERIAL PARAMETERS  - Above & Substrate
    n0      = 1.0;  %index of refraction before first layer   (examples:  n(fused-SiO2)=1.46,  n()=   )
    ns      = 1.0;  %index of refraction after last layer     (examples:  n(TiO2-amorph)=2.3,  n(TiO2-anatase)=2.4, n(TiO2-rutile)=2.6 )\
    kappa_s = 0;  %kappa of after last layer  (substrate) 
    d_subst_nm = 10000000; %thickness of substrate in nm

    
    %MATERIAL PARAMETERS  - LAYERS (from XLS)   
%     TAB_NAME        = 'XDC_LTM0 (from RefPaper)';
    TAB_NAME        = 'Validation1';%'XDC_LTM';
    TAB_NAME        = 'sandbox';
    TAB_NAME        = 'XDC_LTM';
    d_layers                = xlsread('INPUT_coating_prescription.xlsx',TAB_NAME,'b2:b100'); %nm   layer thicknesses
    n                       = xlsread('INPUT_coating_prescription.xlsx',TAB_NAME,'c2:c100'); %nm   refractive index (real part)
    kappas                  = xlsread('INPUT_coating_prescription.xlsx',TAB_NAME,'d2:d100'); %nm   refractive index (imaginary part)
    [dummy,material_names]  = xlsread('INPUT_coating_prescription.xlsx',TAB_NAME,'e2:e100'); %nm   refractive index (imaginary part)
    
%%     
    if max(d_layers)>1000 %if layers > 10um ==> THICK FILM
        BOOL_APPLYSMOOTHING = 1;       
        lam_orig = lam;
        dlam = dlam / N_upsample;
        lam  = [lam_orig(1)-dlam*N_upsample^2: dlam : lam_orig(end)+dlam*N_upsample^2]; %nm ... define wavelength space.
    else
        BOOL_APPLYSMOOTHING=0;
    end
    
    %set up refractive index v. lambda
    if length(material_names)>0  %use material data (if available
        if length(d_layers) ~= length(material_names)
            error('Layer Input counts do not match:  # layer thicknesses different than number of materials.');
        end   
        disp('Using Material Properties....')        
        for kl = 1:length(d_layers)
            disp(['                         ',material_names{kl},'.txt']);
            NK_material{kl} = dlmread([ddir_materials,'/',material_names{kl},'.txt'],'\t',1,0); 
            mat_n(kl,:)     = interp1(NK_material{kl}(:,1),NK_material{kl}(:,2),lam);
            mat_kappa(kl,:) = interp1(NK_material{kl}(:,1),NK_material{kl}(:,3),lam);
            matstack_names{kl,1}    = [material_names{kl},' : ',num2str(d_layers(kl)),'nm' ];
        end
        isel = find(lam>=550,1,'first');NK_material
        n = mat_n(:,isel); %pull vector of single n  per layer (@550nm)
        
    else  %use constant n,k  (no wavelength dependence)
        if  length(d_layers)>0            
            disp('********   USING CONSTANT N, K for each material ********')
            if length(d_layers) ~= length(n) & length(n) ~= length(kappas) 
                error('Layer Input counts do not match:  # layer thicknesses different than number of [n, k] values.');       
            end
            mat_n(kl,:)     = n(kl)*ones(size(lam));
            mat_kappa(kl,:) = kappas(kl)*ones(size(lam));
            matstack_names{kl,1} = ['n=',num2str(n(kl)),' kap=',num2str(kappas(kl)),' : ',num2str(d_layers(kl)),'nm' ];
        else
            mat_n=[]; mat_kappa=[]; matstack_names{1,1} = 'SingleInterface';
        end
    end
    mat_nk(:,:,1) = mat_n;
    mat_nk(:,:,2) = mat_kappa;
    
    
    
%% Loop over Angles

    for k=1:length(theta0_deg_vec)  %loop over angle of incidence         
        theta0_deg = theta0_deg_vec(k);  
        disp(['Angle: ',num2str(theta0_deg_vec(k)),'deg']);
        legnames{k} = [num2str(theta0_deg_vec(k)),'deg'];
                  
        % "jreftan_tr" = T,R,A Calculator (includes absorption loss)
        d_thick = [1000000;d_layers;d_subst_nm];
        for kk = 1:length(lam)
            %d = layer thickness vector, nm
            %n = layer complex refractive index vector
            if length(mat_n)>0
               n_complex_list = [n0;mat_n(:,kk);ns]+ i*[0;mat_kappa(:,kk);kappa_s];
            else
               n_complex_list = [n0;ns]+ i*[0;kappa_s];
            end
            l = lam(kk);                                    %l = free space wavelength, nm
            t0 = deg2rad(theta0_deg);                       %t0= angle of incidence, radians             
            
            polarization = 0; % 0 for TE (s-polarized)            
            [r_s(k,kk),t_s(k,kk),Rs_mat(k,kk),Ts_mat(k,kk),As_mat(k,kk),Mtest2a] ...
                    = jreftran_rt(l,d_thick,n_complex_list,t0,polarization);
            
            polarization = 1; % otherwise for TM (p-polarized)
            [r_p(k,kk),t_p(k,kk),Rp_mat(k,kk),Tp_mat(k,kk),Ap_mat(k,kk),Mtest2b] ...
                    = jreftran_rt(l,d_thick,n_complex_list,t0,polarization);  
        end
        
        %Apply Smoothing
        smoothtype = 'moving'; %'rloess';
        if BOOL_APPLYSMOOTHING %if layers > 10um ==> THICK FILM
            disp('Smoothing Applied');
                tmp.r_p(1,:)    = smooth(r_p(k,:)',N_upsample^2,smoothtype)';
                tmp.t_p(1,:)    = smooth(t_p(k,:)',N_upsample^2,smoothtype)';                
                tmp.Rp_mat(1,:) = smooth(Rp_mat(k,:)',N_upsample^2,smoothtype)';                
                tmp.Tp_mat(1,:) = smooth(Tp_mat(k,:)',N_upsample^2,smoothtype)';                
                tmp.Ap_mat(1,:) = smooth(Ap_mat(k,:)',N_upsample^2,smoothtype)';
                r_p(k,:)    = tmp.r_p(1,:);
                t_p(k,:)    = tmp.t_p(1,:);       
                Rp_mat(k,:) = tmp.Rp_mat(1,:);                
                Tp_mat(k,:) = tmp.Tp_mat(1,:);        
                Ap_mat(k,:) = tmp.Ap_mat(1,:);
                
                tmp.r_s(1,:)    = smooth(r_s(k,:)',N_upsample^2,smoothtype)';
                tmp.t_s(1,:)    = smooth(t_s(k,:)',N_upsample^2,smoothtype)';                
                tmp.Rs_mat(1,:) = smooth(Rs_mat(k,:)',N_upsample^2,smoothtype)';                
                tmp.Ts_mat(1,:) = smooth(Ts_mat(k,:)',N_upsample^2,smoothtype)';                
                tmp.As_mat(1,:) = smooth(As_mat(k,:)',N_upsample^2,smoothtype)';                      
                r_s(k,:)    = tmp.r_s(1,:);
                t_p(k,:)    = tmp.t_s(1,:);       
                Rs_mat(k,:) = tmp.Rs_mat(1,:);                
                Ts_mat(k,:) = tmp.Ts_mat(1,:);        
                As_mat(k,:) = tmp.As_mat(1,:); 
                clear tmp;
        end

        %Original Calculator: brent's 2007 code 
        %   (  "get_ThinFilmSpectrum_model3" / "get_ThinFilmSpectrum_model3"...)
            %         %[Rp,Rs,Tp,Ts,T_tot,T_np,Abs_tot,Abs_np] = get_ThinFilmSpectrum_model3(theta0_deg, theta_pol_deg, d_layers, n,      n0, ns, lam);
            %         [Rp,Rs,Tp,Ts,T_tot,T_np,Abs_tot,Abs_np, Mtest] = get_ThinFilmSpectrum_model4(theta0_deg, 90, d_layers, mat_nk, n0, ns, lam);
            %         Tp_mat(k,:)  = Tp;
            %         Ts_mat(k,:)  = Ts;
            %         Tnp_mat(k,:) = T_np;
            %         Rp_mat(k,:)  = Rp;
            %         Rs_mat(k,:)  = Rs;
            %         Rnp_mat(k,:) = 1-(T_np+Abs_np);        
        
    end
    legnames{1} = ['AOI=',legnames{1}];
    
    
    
    clear tmp;
    if BOOL_APPLYSMOOTHING
        lam = downsample(lam,N_upsample);
        isel = find(lam_orig(1) <= lam & lam <= lam_orig(end)); %keep only original lambda samples
        lam = lam(isel);
        for k = 1:length(theta0_deg_vec)
            tmp.r_s(k,:)    =  downsample(r_s(k,:)',N_upsample)';
            tmp.t_s(k,:)    =  downsample(t_s(k,:)',N_upsample)';
            tmp.Rs_mat(k,:)    =  downsample(Rs_mat(k,:)',N_upsample)';
            tmp.Ts_mat(k,:)    =  downsample(Ts_mat(k,:)',N_upsample)';
            tmp.As_mat(k,:)    =  downsample(As_mat(k,:)',N_upsample)';
            tmp.r_p(k,:)    =  downsample(r_p(k,:)',N_upsample)';
            tmp.t_p(k,:)    =  downsample(t_p(k,:)',N_upsample)';
            tmp.Rp_mat(k,:)    =  downsample(Rp_mat(k,:)',N_upsample)';
            tmp.Tp_mat(k,:)    =  downsample(Tp_mat(k,:)',N_upsample)';
            tmp.Ap_mat(k,:)    =  downsample(Ap_mat(k,:)',N_upsample)';               
        end
        r_s = tmp.r_s(:,isel);  t_s = tmp.t_s(:,isel);  
        Rs_mat = tmp.Rs_mat(:,isel); Ts_mat = tmp.Ts_mat(:,isel); As_mat = tmp.As_mat(:,isel);
        r_p = tmp.r_p(:,isel);  t_p = tmp.t_p(:,isel);  
        Rp_mat = tmp.Rp_mat(:,isel); Tp_mat = tmp.Ts_mat(:,isel); Ap_mat = tmp.Ap_mat(:,isel);
        %clear tmp;
        
    end
    
    %FINAL CALCS
    T_np_mat =  0.5*Tp_mat + 0.5*Ts_mat;
    R_np_mat =  0.5*Rp_mat + 0.5*Rs_mat;
    A_np_mat =  1 - (T_np_mat + R_np_mat);
        
    T_polarized_at_angle = cosd(theta_pol_deg)*Tp_mat + sind(theta_pol_deg)*Ts_mat;
    R_polarized_at_angle = cosd(theta_pol_deg)*Rp_mat + sind(theta_pol_deg)*Rs_mat;
    %========================================================================================     
    

    %FIGURE...
    figure; plot(lam,R_np_mat);  xlabel('wavelength [nm]');
    ylim([0,1]); legend(legnames);  ylabel('Reflectivity (non-polarized light)');
    title('Reflectivity Spectra of Nonpolarized Light at many Angles of Incidence for this stack:');
    text(550, 0.95, matstack_names,'VerticalAlignment','top') 
 
    %DATA SAVE....
    save(['SpectrumResults_',date,'.mat'],'theta0_deg_vec','lam','d_layers','n','n0','ns','material_names','mat_n','mat_kappa',...
            'Tp_mat','Ts_mat','T_np_mat','Rp_mat','Rs_mat','R_np_mat','Ap_mat','As_mat','A_np_mat','T_polarized_at_angle','R_polarized_at_angle');

    %Export in Zemax "TABLE" format (for coating.dat)    
    %WAVE <wavelength 1 in micrometers> Rs Rp Ts Tp Ars Arp Ats Atp
    dlam_zmx = 20;lam_zmx = [lam(1):dlam_zmx:lam(end)];
    isel = interp1(lam,1:length(lam),lam_zmx); %indices of selected lambda samples
    table_matrix = [];
    for k=1:length(theta0_deg_vec)  
         table_matrix_thisangle = [lam_zmx/1000; Rs_mat(k,isel); Rp_mat(k,isel); Ts_mat(k,isel); Tp_mat(k,isel); ]';
         table_matrix = [table_matrix;theta0_deg_vec(k)*ones(3,5);table_matrix_thisangle];
    end
    csvwrite('temp_table_matrix.csv',table_matrix);
   
    
    %--------------------------------------------------------------------------------------------------
  %  return;
    %comment out this return to keep running code below...
    %--------------------------------------------------------------------------------------------------
    
   
    fdate = date;  %***UPDATE*****
    
    load(['SpectrumResults_',fdate,'.mat']);

    
    %Plot 2D Reflectivity Spectrum v. Angle (NonPolarized Light)
    temp(:,:,1) = ones(size(R_np_mat));
    temp(:,:,2) = R_np_mat;
    temp2 = min(temp,3);
    figure; surf(lam,theta0_deg_vec,R_np_mat); xlabel('wavelength [nm]'); ylabel('incidence angle (deg)'); zlabel('R := |Eref/Ei|^2'); 
    shading interp; xlim([min(lam),1000]); 
    title('Reflectance of Non Polarized Light'); 
    
    
    %Calculate Reflectivity per *** ASTM E971 ***
    A = xlsread([ddir,'/Solar Spectra.xlsx'],'Sheet1','a2:b4000');
    dlam = median(lam(2:end)-lam(1:end-1));
    Ysolar = interp1(A(:,1),A(:,2),lam); %[W/m2/nm]
    
    A = csvread([ddir,'/PhotopicSpectrum.Csv']);
    Yphotopic = interp1(A(:,1),A(:,2),lam); %[a.u.]
    isel_nan = find(isnan(Yphotopic));
    Yphotopic(isel_nan) = 0;
    
    for k=1:size(R_np_mat,1) %loop over AOI
        SolarANDPhotopicWeightFactors(k,:)     = ones(1,length(lam)).*Ysolar.*Yphotopic/sum(Ysolar.*Yphotopic);
        R_np_mat_SolarANDPhotopicWeighted(k,:) = R_np_mat(k,:)      .*Ysolar.*Yphotopic/sum(Ysolar.*Yphotopic);
        %integrate over wavelengths
        R_np_E971(k,1) = sum( R_np_mat(k,:).*Ysolar.*Yphotopic/sum(Ysolar.*Yphotopic) ); 
    end
   
    figure; surf(lam,theta0_deg_vec,SolarANDPhotopicWeightFactors); xlabel('wavelength [nm]'); ylabel('incidence angle (deg)');
    zlabel('R := |Eref/Ei|^2  **per nm**');     shading interp; xlim([min(lam),1000]); 
    title('Weighting Factors = AM1.5/sum(AM1.5) x  Photopic/sum(Photopic)'); 
        %double check:  weighting factors sum to 1: 
        disp('double check:  weighting factors sum to 1:');
        sum(SolarANDPhotopicWeightFactors,2)    
    figure; surf(lam,theta0_deg_vec,R_np_mat_SolarANDPhotopicWeighted); xlabel('wavelength [nm]'); ylabel('incidence angle (deg)');
    zlabel('R := |Eref/Ei|^2  **per nm**');     shading interp; xlim([min(lam),1000]); 
    title('Reflectance of Non Polarized Light: Weighted by Solar Spectrum AND Photopic Response of Eye'); 
    
    %Reflectivity v. Angle
    figure; plot(theta0_deg_vec,R_np_E971); xlabel('theta incidence (deg)'); ylabel('E971-Weighted Reflectivity nonpolarized Light'); grid on;
    title({'E971-Weighted Reflectivity (Weighting = AM1.5 x Photopic)';' For Light Trap Metal Stack:'});
    text(10,max(R_np_E971)*0.9, matstack_names,'VerticalAlignment','top') 
 
% 
%     %Integrate over Wavelength to get Rel.Intensity = f(theta_incidence)
%     %where rel.intensity := Watts_Transmitted / Watts_incident , integrated
%     %over all wavelengths.
%     %The Watts_Incident is given by solar insolation data: 
%     A = xlsread('Solar Spectra.xls','Sheet1','a2:b4000');
%     dlam = median(lam(2:end)-lam(1:end-1))
%     Y0=0; Y = zeros(size(Tnp_mat,1),1); 
%     for k=1:length(lam)
%         indx_hi = find( lam(k) <= A(:,1), 1, 'first' );
%         indx_lo = find( lam(k) > A(:,1), 1, 'last' ); %should be: indx_lo = indx_hi-1;
% 
%         %interpolate (linearly) if lam(k) is between datapoints of spectrum data file...
%         if length(indx_lo) == 1 & length(indx_hi) == 1
%             value = ( (lam(k)-A(indx_lo,1))*A(indx_hi,2) + (A(indx_hi,1)-lam(k))*A(indx_lo,2) ) / (A(indx_hi,1)-A(indx_lo,1));
%         elseif length(indx_lo) == 1 & length(indx_hi) == 0
%             value = A(indx_lo,2);
%         elseif length(indx_lo) == 0 & length(indx_hi) == 1
%             value = A(indx_hi,2);
%         else
%             indx = [find( A(:,1) > lam(k) , 1, 'first') , find( A(:,1) < lam(k) , 1, 'last')];
%             warning(strcat('lambda = ',num2str(lam(k)),' is outside the sampling domain for the spectrum!! (',num2str(A(1,1)),'nm...',num2str(A(end,1)),'nm)'));
%         end
% 
%         Y=Y+Tnp_mat(:,k)*(value*dlam);  %integrate this lambda contribution
%         Y0=Y0+(value*dlam); %integrate this lambda contribution
%     end
% 
%     %Plot Transmission v. Incidence Angle (for *NONPOLARIZED*, *AM1.5 SPECTRUM* light!!)
%     T_np_AM15 = Y/Y0;
%     figure; plot(theta0_deg_vec,T_np_AM15); xlabel('theta incidence (deg)'); ylabel('Transmissivity of AR coating'); grid on;
%     title({'Field of View of 100% efficient flat-panel Solar Cell';' with AR coating under *NonPolarized* *AM1.5* Insolation'});



end %function "run_coating_Model"







