% ========================================================================
% Fronteras de Estabilidad para convertidor CD-CD Boost
% Autor : Erick Moreno Negrete
% Version : 1.3
% ========================================================================

clc
format longG
set(groot,'defaultAxesTickLabelInterpreter','latex');


% ========================================================================
%                   Parametros del Sistema
% ========================================================================

L = 2.7648e-3;
C = 1.666667e-6;
R = 144;
E = 48;
d = 0.6;
dp = 1-d;
Vo = E/(1-d); % Equilibrio de Voltaje
IL = E/(R*((1-d)^2)); % Equilibrio de Corriente

% ========================================================================
%                     Modelo de pequeña señal
% ========================================================================

A = [   0     -dp/L;
     dp/C  -1/(C*R)];

F = [Vo/L;
     -IL/C];
 
C1 = [1 0];
C2 = [0 1];

D = 0;

% ========================================================================
%          Funciones de transferencia C1=corriente; C2=voltaje
% ========================================================================

[num,den] = ss2tf(A,F,C2,D);

numC = [0 Vo/L (2*dp*IL)/(L*C)];
denC = [1 1/(R*C) (dp^2)/(L*C)];

numV = [0 -IL/C (Vo*dp)/(L*C)];
denV = [1 1/(R*C) (dp^2)/(L*C)];

G = tf(num,den) 

Gi = tf(numC,denC)
Gv = tf(numV,denV)

% ========================================================================
%     Coeficientes de la funcion de transferencia de corriente Gi
% ========================================================================

n1 = numV(2);
n0 = numV(3);
d2 = denV(1);
d1 = denV(2);
d0 = denV(3);

% ========================================================================
%       Caracterizacion de las fronteras de estabilidad Gc
% ========================================================================

sigma_inicial = 0; 
sigma_final = 2000;
inc = 200;
sgma = sigma_inicial:inc:sigma_final;



 for iter = 1:1:length(sgma)
     
     sigma = sgma(iter);

    w = linspace(0,100000,100000);
    sb = -sigma+1i.*w;
    
    %Polinomios en s=sigma
    Nds0 = n1*(sigma) + n0;
    Dds0 = d2*((sigma)^2) + d1*(sigma) + d0;
    
    %Polinomios en s=sigma
    Nds = n1*sb + n0;
    Dds = d2*sb.^2 + d1*sb + d0;
    
    figure(1);
    hold on
    %axis([-0.1e-4 6e-5 -0.1 .5]);
    xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
    ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
    title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');



    if sigma==0
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        Kp_r = linspace(-5e-3,0.3e-3,1000);
        Ki_r = sigma*Kp_r + sigma*(Dds/Nds);
        Ki_r2 = linspace(-0.2,.6,1000);
        figure(1)
        plot(vKps,vKis,'k--')
        plot(Kp_r,Ki_r,'k--')
        plot(Ki_r,Ki_r2,'k--')
    end


    if sigma==200
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#7400b8';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)

    end
 
    
    if sigma==400
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#6930c3';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end
    
    if sigma==600
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#5e60ce';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end
    
    if sigma==800
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#5390d9';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==1000
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#4ea8de';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==1200
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#48bfe3';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==1400
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#56cfe1';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==1600
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#64dfdf';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==1800
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#72efdd';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
    
    if sigma==2000
       
        vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds));
        vKis = (w+((sigma^2)./w)).*imag(Dds./Nds);
        
        str = '#80ffdb';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        
        figure(1)
        plot(vKps,vKis,'Color', color)
    end 
   
 
end

% ========================================================================
%       Caracterizacion de las fronteras de estabilidad Gv
% ========================================================================

% sigma_inicial2 = 0; 
% sigma_final2 = 10000;
% inc2 = 1000;
% sgma2 = sigma_inicial2:inc2:sigma_final2;
% 
% 
% for iter = 1:1:length(sgma2)
%     
%     sigma = sgma2(iter);
%     w = linspace(0,1000000,100);
%     sb = -sigma+1i.*w;
%     
%     Polinomios en s=sigma
%     Nds0 = n1v*(-sigma) + n0v;
%     Dds0 = d2*((-sigma)^2) + d1*(-sigma) + d0;
%     
%     Polinomios en s=sigma
%     Nds = n1v*sb + n0v;
%     Dds = d2*(sb.^2) + d1*sb + d0;
%     
%     figure(2);
%     hold on
%     axis([-0.1e-4 6e-5 -0.1 .5]);
%     xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
%     ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
%     title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');
% 
%         vKps = -real(Dds./Nds)-((sigma./w).*imag(Dds./Nds))
%         vKis = (w+((sigma^2)./w)).*imag(Dds./Nds)
%         
%         Kp_r = linspace(-5e-3,0.3e-3,1000);
%         Ki_r = sigma*Kp_r + sigma*(Dds/Nds);
%         Ki_r2 = linspace(-0.2,.6,1000);
%         figure(2)
%         
%         plot(vKps,vKis,'k')
%         hold on
%         plot(Kp_r,Ki_r,'k')
%         plot(Ki_r,Ki_r2,'k--')
% 
% end
 

% ========================================================================
%       Caracterizacion de las fronteras de estabilidad Gv
% ========================================================================


% 
% while (1)
%     figure(1)
%     [kp,ki]=ginput(1)
%     
%     c3 = d2;
%     c2 = (kp*n1 + d1);
%     c1 = (ki*n1 + kp*n0 + d0);
%     c0 = ki*n0;
%     
%     Pcl=[c3 c2 c1 c0];
%     roots(Pcl)
%     
%     Gc=pid(kp,ki);
%     sys=series(Gc,G);
%     H=1;
%     OL=G*Gc;
%     Mc=feedback(sys,H);
%     
%     figure(2)
%     step(Mc)
% end
