% ========================================================================
% Fronteras de Estabilidad para convertidor CD-CD Boost
% Autor : Erick Moreno Negrete
% Version : 1.4
% ========================================================================

clc
format longG



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



    w = linspace(0,10000,1000);
    sb = 1i.*w;
    
    s1 = -500+(1i.*w);
    s1n = -500;
    
    s2 = -1000+(1i.*w);
    s2n = -1000;
    
    s3 = -1300+(1i.*w);
    s3n = -1300;
    
    %Polinomios en s=jw; fronteras de estabilidad
    Nds0 = n1*sb + n0;
    Dds0 = d2*sb.^2 + d1*sb + d0;
    
    %Polinomios en s=-500+jw; 
    Nds1 = n1*s1 + n0;
    Dds1 = d2*s1.^2 + d1*s1 + d0;
    
    %Polinomios en s=-500; 
    Nds1n = n1*s1n + n0;
    Dds1n = d2*s1n.^2 + d1*s1n + d0;
    
    %Polinomios en s=-1000+jw; 
    Nds2 = n1*s2 + n0;
    Dds2 = d2*s2.^2 + d1*s2 + d0;
    
    %Polinomios en s=-1000; 
    Nds2n = n1*s2n + n0;
    Dds2n = d2*s2n.^2 + d1*s2n + d0;
    
    %Polinomios en s=-2000+jw; 
    Nds3 = n1*s3 + n0;
    Dds3 = d2*s3.^2 + d1*s3 + d0;
    
    %Polinomios en s=-2000; 
    Nds3n = n1*s3n + n0;
    Dds3n = d2*s3n.^2 + d1*s3n + d0;
    
    figure(1);
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(gca,'GridLineStyle','--')
    
    

       
        vKp0 = -real(Dds0./Nds0);
        vKi0 = w.*imag(Dds0./Nds0);
        rKp0 = linspace(-4e-3,5e-3,1000);
        rKi0 = linspace(0,0,1000);
        rKi0y = linspace(-1,13,1000);
        
        
        plot(vKp0,vKi0,'k--')
        hold on
        plot(rKp0,rKi0,'k--','HandleVisibility','off')
        plot(rKi0,rKi0y,'k--','HandleVisibility','off')
        
        % sigma = 500
        vKp1 = -real(Dds1./Nds1);
        vKi1 = w.*imag(Dds1./Nds1);
        
        %recta 
        rKp1 = linspace(-4e-3,5e-3,1000);
        rKi1 = (500*rKp1)+500*(Dds1n/Nds1n);
    
        plot(vKp1,vKi1,'k-.')
        plot(rKp1,rKi1,'k-.','HandleVisibility','off')
        
        % sigma = 1000
        vKp2 = -real(Dds2./Nds2);
        vKi2 = w.*imag(Dds2./Nds2);
        
        %recta 
        rKp2 = linspace(-4e-3,5e-3,1000);
        rKi2 = (1000*rKp2)+1000*(Dds2n/Nds2n);
    
        plot(vKp2,vKi2,'k:')
        plot(rKp2,rKi2,'k:','HandleVisibility','off')
        axis([-4e-3 5e-3 -1 13]);
        xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
        ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
        title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');

        legend('$\sigma=0$','$\sigma<500$','$\sigma<1000$','$\sigma<1300$','interpreter', 'latex')
        % sigma = 2000
        vKp3 = -real(Dds3./Nds3);
        vKi3 = w.*imag(Dds3./Nds3);
        
        %recta 
        rKp3 = linspace(-4e-3,5e-3,1000);
        rKi3 = (1000*rKp3)+1000*(Dds3n/Nds3n);
    
%         plot(vKp3,vKi3,'g')
%         plot(rKp3,rKi3,'g')

 

% ========================================================================
%       Caracterizacion de las fronteras de estabilidad Gv
% ========================================================================


% % 
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
