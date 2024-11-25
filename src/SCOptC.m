%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [models_out,N,D,results,foptions_out] = SCOptC(wl,theta,models,foptions)
    
    %% UltimaRI
    
    %% Files
    % Variables inside the data files must be called:

    %  Size of Reflectance & Transmittance data must be (  length(wl_exp) ,  length(theta_exp) )
    %
    
    models_out = models;
    foptions_out = foptions;


    %% Refractive Index Models
    % Different refractive index models can be used for each unknown layer
    % The following models can be used:
    %
    % 'Fh-N' Forouhi Bloomer model  ( Eg , n0 , fi , Ei , Gi ) with N oscillators (max N=4):
    %
    % 'Ch-n' real Cauchy model      ( A1 , A2 , A3 )
    %
    % 'Ch-nk' complex Cauchy model  ( A1 , A2 , A3 , A4 , A5 , A6 )
    %
    % 'cnst' constant refractive index ( n0 )
    %
    % '.mat' import data from file (only for known layers)
    %       Variables inside data files containing RI data must be:
    %       'wl_exp'  for wavelength
    %       'n'       for real part of the RI
    %       'k'       for imaginary part of the RI
    %

    lcoher = foptions.lcoher;
    n_layers = length(models);
    %N = zeros(length(wl),n_layers);
    %D = zeros(n_layers,1);


    
    auxind = 0;
    for k2=1:length(models)
        switch models{k2}.type
            case "lin-grad"
                auxind = auxind+1;
                models{k2}.index = auxind;
                n1 = models{k2}.n1;
                n2 = models{k2}.n2;
                nlayers = models{k2}.nlayers;
                nvec = linspace(n1,n2,nlayers);

                for jj=1:nlayers
                    N(:,models{k2}.index+jj-1) = ones(length(wl),1)*nvec(jj);
                    D(models{k2}.index+jj-1) = models{k2}.D/nlayers/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index+jj-1) = 0;
                    else
                        active(models{k2}.index+jj-1) = 1;
                    end
                end
                auxind = auxind+nlayers-1;

             case "DBR"
                auxind = auxind+1;
                models{k2}.index = auxind;
                n1 = models{k2}.n1;
                n2 = models{k2}.n2;
                D1 = models{k2}.D1;
                D2 = models{k2}.D2;
                nperiod = models{k2}.nperiod;
                nlayers = nperiod*2+1;

                for jj=1:2:(nlayers-1)
                    N(:,models{k2}.index+jj-1) = ones(length(wl),1)*n1;
                     N(:,models{k2}.index+jj-1+1) = ones(length(wl),1)*n2;
                    D(models{k2}.index+jj-1) = D1/1000;
                    D(models{k2}.index+jj-1+1) = D2/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index+jj-1) = 0;
                        active(models{k2}.index+jj-1+1) = 0;
                    else
                        active(models{k2}.index+jj-1) = 1;
                        active(models{k2}.index+jj-1+1) = 1;
                    end
                end
                N(:,models{k2}.index+nlayers-1) = ones(length(wl),1)*n1;
                D(models{k2}.index+nlayers-1) = D1/1000;
                if models{k2}.active == "false"
                        active(models{k2}.index+nlayers-1) = 0;
                else 
                        active(models{k2}.index+nlayers-1) = 1;
                end
                auxind = auxind+nlayers-1;

           
            case "Fh-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                Eg = models{k2}.Eg;
                n0 = models{k2}.n0;
                fi = models{k2}.fi;
                Ei = models{k2}.Ei;
                Gi = models{k2}.Gi;
                N(:,models{k2}.index) = f_nk_ForouhiBloomer(wl,Eg, n0, fi, Ei, Gi);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "Lnz-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                E0 = models{k2}.E0;
                n0 = models{k2}.n0;
                Ep = models{k2}.Ep;
                g = models{k2}.g;
                N(:,models{k2}.index) = f_nk_lorentz(wl, n0, E0, Ep, g);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "eFh-N"  % 5 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                Eg = models{k2}.Eg;
                n0 = models{k2}.n0;
                fi = models{k2}.fi;
                Ei = models{k2}.Ei;
                Gi = models{k2}.Gi;
                n_pvk = f_nk_ForouhiBloomer(wl,Eg, n0, fi, Ei, Gi);
                ff_pvk = models{k2}.ff_pvk;
                n_mat = models{k2}.n_mat;
                ff_mat = models{k2}.ff_mat;
                N(:,models{k2}.index) = f_nk_EMA(n_mat,n_pvk,1.0, ff_mat, ff_pvk ,2);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "Ch-n"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                A = models{k2}.A;
                N(:,models{k2}.index) = Cauchy_n(wl,A);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "Ch-nk"   % 3 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                A = models{k2}.A;
                N(:,models{k2}.index) = Cauchy_nk(wl,A);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "cnst" % 1 + 1 par
                auxind = auxind+1;
                models{k2}.index = auxind;
                ncnst = models{k2}.n;
                N(:,models{k2}.index) = ones(length(wl),1)*ncnst;
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end
            case "file"
                auxind = auxind+1;
                models{k2}.index = auxind;
                nkdata = load(models{k2}.filename);
                N(:,models{k2}.index) = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);
                if k2~=1 && k2~=n_layers
                    D(models{k2}.index) = models{k2}.D/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index) = 0;
                    else
                        active(models{k2}.index) = 1;
                    end
                end

             case "DBRf"
                auxind = auxind+1;
                models{k2}.index = auxind;

                nkdata = load(models{k2}.filename1);
                n1 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                nkdata = load(models{k2}.filename2);
                n2 = interp1(nkdata.wl,nkdata.n,wl)+1i*interp1(nkdata.wl,nkdata.k,wl);

                D1 = models{k2}.D1;
                D2 = models{k2}.D2;
                nperiod = models{k2}.nperiod;
                nlayers = nperiod*2+1;

                for jj=1:2:(nlayers-1)
                    N(:,models{k2}.index+jj-1) = n1;
                     N(:,models{k2}.index+jj-1+1) = n2;
                    D(models{k2}.index+jj-1) = D1/1000;
                    D(models{k2}.index+jj-1+1) = D2/1000;
                    if models{k2}.active == "false"
                        active(models{k2}.index+jj-1) = 0;
                        active(models{k2}.index+jj-1+1) = 0;
                    else
                        active(models{k2}.index+jj-1) = 1;
                        active(models{k2}.index+jj-1+1) = 1;
                    end
                end
                N(:,models{k2}.index+nlayers-1) = n1;
                D(models{k2}.index+nlayers-1) = D1/1000;
                if models{k2}.active == "false"
                        active(models{k2}.index+nlayers-1) = 0;
                else 
                        active(models{k2}.index+nlayers-1) = 1;
                end

                auxind = auxind+nlayers-1;
                
            
        end
    end



    D = D(2:end);
    active = active(2:end);

    D = D'*1000;

%% Plot

    if foptions.backwards == "true"
        N = fliplr(N);
        D = flipud(D);
        active = fliplr(active);
    end
    active = active>0;
    
    z = -200:foptions.zstep:0;
    if isempty(foptions.zstep) == false
        h = [0 ; cumsum(D)];
        for ww = 1 : length(D)
            if D(ww)<lcoher
                z = [z , h(ww):foptions.zstep:h(ww+1)];
            else
                z = [z , linspace(h(ww),h(ww+1),20)];
            end
        end
        z = [z , z(end)+1:z(end)+200];
    end

    [Rt_S,Rt_P,Tt_S,Tt_P,E2_s,E2_p,Pabs,Abs] = f_plot_RTc(N, D, lcoher, wl, theta.values,z,active,foptions);
    results.Rt_S = Rt_S;
    results.Rt_P = Rt_P;
    results.Tt_S = Tt_S;
    results.Tt_P = Tt_P;
    results.Pabs = Pabs;
    results.Abs = Abs;
    results.z = z;
    results.active = active;

    %% Calculate Jsc
    load AM15_spec.mat AM15
    load LED_3000K.mat wl_LED VarName540
    load LED_6500K.mat VarName1412
    qe        = 1.602176565e-19;
    h         = 6.62606957e-34;
    c         = 299792458;


    Irr_n = AM15(:,4)./trapz(AM15(:,1), AM15(:,4))*1000; % 1000 W/m2
    Iam15 = interp1(AM15(:,1), Irr_n, wl);
    Nam15 = Iam15.*wl*1e-9/h/c;
    Nam15 = Nam15.';
    Iam15 = Iam15.';

    Irr_led1 = VarName540./trapz(wl_LED, VarName540)*4; % 4 W/m2
    Iled1 = interp1(wl_LED, Irr_led1, wl);
    N3000K = Iled1.*wl*1e-9/h/c;
    N3000K = N3000K.';

    Irr_led2 = VarName1412./trapz(wl_LED, VarName1412)*4; % 4 W/m2
    Iled2= interp1(wl_LED, Irr_led2, wl);
    N6500K = Iled2.*wl*1e-9/h/c;
    N6500K = N6500K.';

    IPCE     = squeeze(Abs(:,:,active));
    IPCE_SQ =(wl<(1239.8/foptions.SQgap)).*ones(1,length(wl));
    Jsc_SQ = trapz(wl, qe.*IPCE_SQ.'.*Nam15)/1e4*1e3;
    results.IPCE_SQ = IPCE_SQ;
    results.Jsc_SQ = Jsc_SQ;
    Jsc      = trapz(wl, qe.*IPCE.*Nam15)/1e4*1e3;
    results.Jsc_total = sum(Jsc);
    results.IPCE = IPCE;
    results.Jsc = Jsc;
    IPCE_loss     = squeeze(Abs(:,:,~active));
    Jsc_loss      = trapz(wl, qe.*IPCE_loss.*Nam15)/1e4*1e3;
    results.Jsc_loss_total = sum(Jsc_loss);
    results.IPCE_loss = IPCE_loss;
    results.Jsc_loss = Jsc_loss;

%     Jsc_3000K      = trapz(wl, qe.*IPCE.*N3000K)/1e4*1e3;
%     results.Jsc_3000K = Jsc_3000K;
%     Jsc_loss_3000K      = trapz(wl, qe.*IPCE_loss.*N3000K)/1e4*1e3;
%     results.Jsc_loss_3000K = Jsc_loss_3000K;
% 
%     Jsc_6500K      = trapz(wl, qe.*IPCE.*N6500K)/1e4*1e3;
%     results.Jsc_6500K = Jsc_6500K;
%     Jsc_loss_6500K      = trapz(wl, qe.*IPCE_loss.*N6500K)/1e4*1e3;
%     results.Jsc_loss_6500K = Jsc_loss_6500K;

    load CIE_photopic.mat CIE1931 wl_CIE

    CIE = interp1(wl_CIE,CIE1931,wl);
    Tt = 0.5*(Tt_P+Tt_S);

    Sun_Tt = Iam15.*Tt;

    aux_wl = wl;
    aux_Tt = Sun_Tt;

    % Build CRI data
    if wl(1)<380 || wl(end)>780
        aux_wl = (380:780).';
        aux_Tt = interp1(wl,Sun_Tt,aux_wl);
        aux_Iam15 = interp1(wl,Iam15,aux_wl);
    end

    

    [~,~,~,~,~,~,Ra,~] =pspectro([aux_wl ,aux_Tt]);

    [Ra_Sun,~] =getcri1995([aux_wl ,aux_Tt],[aux_wl ,aux_Iam15],aux_wl.');

    results.CRI = Ra/100;
    results.CRI_Sun = Ra_Sun/100;

    AVT = trapz(wl,CIE.'.*Iam15.*Tt,1)./trapz(wl,CIE.'.*Iam15); % Average visible transmission
    AT = trapz(wl,Iam15.*Tt,1)./trapz(wl,Iam15); % Average visible transmission
    results.AVT = AVT;
    results.AT = AT;

%     AVT_3000K = trapz(wl,CIE.'.*N3000K.*Tt,1)./trapz(wl,CIE.'.*N3000K); % Average visible transmission
%     AT_3000K = trapz(wl,N3000K.*Tt,1)./trapz(wl,N3000K); % Average visible transmission
%     results.AVT_3000K = AVT_3000K;
%     results.AT_3000K = AT_3000K;
% 
%     AVT_6500K = trapz(wl,CIE.'.*N6500K.*Tt,1)./trapz(wl,CIE.'.*N6500K); % Average visible transmission
%     AT_6500K = trapz(wl,N6500K.*Tt,1)./trapz(wl,N6500K); % Average visible transmission
%     results.AVT_6500K = AVT_6500K;
%     results.AT_6500K = AT_6500K;

    results.E2 = 0.5*(E2_s+E2_p);


end


