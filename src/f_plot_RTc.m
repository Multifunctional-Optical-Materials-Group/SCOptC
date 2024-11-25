%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   This work is licensed under the Creative Commons Reconocimiento-NoComercial-CompartirIgual 4.0 Internacional License.
%   To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Rt_S,Rt_P,Tt_S,Tt_P,E2_s,E2_p,Pabs,Pabs_s,Pabs_p,Abs,Abs_s,Abs_p] = f_plot_RTc(N, D, lcoher, wl, theta,z,active,foptions)
    
   

    Tt_S=zeros(length(wl),length(theta));
    Tt_P=zeros(length(wl),length(theta));
    Rt_S=zeros(length(wl),length(theta));
    Rt_P=zeros(length(wl),length(theta));


    E2_s = zeros(length(wl),length(theta),length(z));
    E2_p = zeros(length(wl),length(theta),length(z));
    Pabs_s = zeros(length(wl),length(theta),length(z));
    Pabs_p = zeros(length(wl),length(theta),length(z));

    Abs = zeros(length(wl),length(theta),length(D));
    Abs_s = zeros(length(wl),length(theta),length(D));
    Abs_p = zeros(length(wl),length(theta),length(D));


    for k1=1:length(theta)
        for k2=1:length(wl)
            
            [Rt_S(k2,k1), Rt_P(k2,k1), Tt_S(k2,k1), Tt_P(k2,k1), E2_s(k2,k1,:), E2_p(k2,k1,:), Pabs_s(k2,k1,:), Pabs_p(k2,k1,:)] = RTF_Abeles_F(N(k2,:), D',[], wl(k2),theta(k1)*pi/180,z,lcoher,30);
            
        end

    end

    Rt = 0.5*(Rt_S + Rt_P);
    Tt = 0.5*(Tt_S + Tt_P);
    Pabs =0.5*(Pabs_s+Pabs_p);

    h = [0 ; cumsum(D)];
    for w = 1:length(D)
        if D(w)<lcoher
            [~,z1] = min(abs(h(w)-z));
            [~,z2] = min(abs(h(w+1)-z));
            Abs(:,:,w) = trapz(z(z1:z2),Pabs(:,:,z1:z2),3);
            Abs_s(:,:,w) = trapz(z(z1:z2),Pabs_s(:,:,z1:z2),3);
            Abs_p(:,:,w) = trapz(z(z1:z2),Pabs_p(:,:,z1:z2),3);
        end
    end

    if foptions.plot == "true"
        figure;
        clf
        for jj=1:length(theta)
            subplot(1,3,1)
            plot(wl,Tt(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
            xlabel("\lambda (nm)")
            ylabel("T")
            xlim([min(wl),max(wl)])
    
            subplot(1,3,2)
            plot(wl,Rt(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
            xlabel("\lambda (nm)")
            ylabel("R")
            xlim([min(wl),max(wl)])
    
            subplot(1,3,3)
    
            plot(wl,1-Tt(:,jj)-Rt(:,jj),'LineWidth',1.5,'color','k','LineStyle','-')
            hold on
            for ww = 1:length(active)
                if active(ww) == true
                    plot(wl,squeeze(Abs(:,jj,ww)),'LineWidth',1.5,'color','r','LineStyle','--')
                else
                    plot(wl,squeeze(Abs(:,jj,ww)),'LineWidth',1.5,'color','b','LineStyle','--')
                end
            end
            plot(wl,sum(Abs(:,:,active),3),'LineWidth',1.5,'color','r','LineStyle','-')
            plot(wl,sum(Abs(:,:,~active),3),'LineWidth',1.5,'color','b','LineStyle','-')

            xlabel("\lambda (nm)")
            ylabel("A")
            xlim([min(wl),max(wl)])
    
        end
        
        set(gcf,'units','points','position',[150,150,900,300])
    end
    
    
end