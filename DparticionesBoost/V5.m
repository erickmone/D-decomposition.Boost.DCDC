% ========================================================================
% Fronteras de Estabilidad para convertidor CD-CD Boost
% Autor : Erick Moreno Negrete
% Version : 1.5 
% Lazo interno de corriente
% ========================================================================

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

[num,den] = ss2tf(A,F,C1,D);
G = tf(num,den) 

numC = [0 Vo/L (2*dp*IL)/(L*C)];
numV = [0 -IL/C (Vo*dp)/(L*C)];
den = [1 1/(R*C) (dp^2)/(L*C)];

G0 = tf(numV,den)
G1 = tf(numC,den)
G2 = tf(numV,numC)

% ========================================================================
%     Coeficientes de la funcion de transferencia de corriente G1
% ========================================================================

n11 = numC(2);
n10 = numC(3);
d12 = den(1);
d11 = den(2);
d10 = den(3);

% ========================================================================
%     Coeficientes de la funcion de transferencia de voltaje G0
% ========================================================================

n01 = numV(2);
n00 = numV(3);
d02 = den(1);
d01 = den(2);
d00 = den(3);

% ========================================================================
%     Coeficientes de la funcion de transferencia de corriente G2
% ========================================================================

n21 = numV(2);
n20 = numV(3);
d21 = numC(2);
d20 = numC(3);

% ========================================================================
%       Caracterizacion de las fronteras de estabilidad G1
% ========================================================================

% Fronteras de sigma
%%%%%%%

% w = linspace(0,12000,10000);
% sb = 1i.*w;
%     
% s1 = -100+(1i.*w);
% s1n = -100;
%     
% s2 = -500+(1i.*w);
% s2n = -500;
%     
% s3 = -1000+(1i.*w);
% s3n = -1000;
% 
% s4 = -2000+(1i.*w);
% s4n = -2000;
% 
% s5 = -3000+(1i.*w);
% s5n = -3000;
% 
% s6 = -6000+(1i.*w);
% s6n = -6000;
% 
% %Polinomios en s=jw; fronteras de estabilidad
% Nds0 = n11*sb + n10;
% Dds0 = d12*sb.^2 + d11*sb + d10;
% 
% %Polinomios en s=-100+jw; 
% Nds1 = n11*s1 + n10;
% Dds1 = d12*s1.^2 + d11*s1 + d10;
% 
% %Polinomios en s=-100; 
% Nds1n = n11*s1n + n10;
% Dds1n = d12*s1n.^2 + d11*s1n + d10;
% 
% %Polinomios en s=-500+jw; 
% Nds2 = n11*s2 + n10;
% Dds2 = d12*s2.^2 + d11*s2 + d10;
%     
% %Polinomios en s=-500; 
% Nds2n = n11*s2n + n10;
% Dds2n = d12*s2n.^2 + d11*s2n + d10;
% 
% %Polinomios en s=-1000+jw; 
% Nds3 = n11*s3 + n10;
% Dds3 = d12*s3.^2 + d11*s3 + d10;
%     
% %Polinomios en s=-1000; 
% Nds3n = n11*s3n + n10;
% Dds3n = d12*s3n.^2 + d11*s3n + d10;
% 
% %Polinomios en s=-2000+jw; 
% Nds4 = n11*s4 + n10;
% Dds4 = d12*s4.^2 + d11*s4 + d10;
%     
% %Polinomios en s=-2000; 
% Nds4n = n11*s4n + n10;
% Dds4n = d12*s4n.^2 + d11*s4n + d10;
% 
% %Polinomios en s=-3000+jw; 
% Nds5 = n11*s5 + n10;
% Dds5 = d12*s5.^2 + d11*s5 + d10;
%     
% %Polinomios en s=-3000; 
% Nds5n = n11*s5n + n10;
% Dds5n = d12*s5n.^2 + d11*s5n + d10;
% 
% %Polinomios en s=-4000+jw; 
% Nds6 = n11*s6 + n10;
% Dds6 = d12*s6.^2 + d11*s6 + d10;
%     
% %Polinomios en s=-4000; 
% Nds6n = n11*s6n + n10;
% Dds6n = d12*s6n.^2 + d11*s6n + d10;
%     
% vKp0 = -real(Dds0./Nds0);
% vKi0 = w.*imag(Dds0./Nds0);
% rKp0 = linspace(-0.1,0.25,1000);
% rKi0 = linspace(0,0,1000);
% rKi0y = linspace(-1,1800,1000);
% 
% figure(1);
% set(groot,'defaultAxesTickLabelInterpreter','latex');
% set(gca,'GridLineStyle','--')
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#1a759f';
% color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% % sigma = 0
% plot(vKp0,vKi0,'Color',color1)
% hold on
% plot(rKp0,rKi0,'b--','HandleVisibility','off')
% plot(rKi0,rKi0y,'b--','HandleVisibility','off')
% 
% % sigma = 100
% vKp1 = -real(Dds1./Nds1);
% vKi1 = w.*imag(Dds1./Nds1);
% 
% %recta 
% rKp1 = linspace(-0.1,0.25,1000);
% rKi1 = (50*rKp1)+50*(Dds1n/Nds1n);
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#168aad';
% color2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp1,vKi1,'Color',color2)
% plot(rKp1,rKi1,'Color',color2,'HandleVisibility','off')
% 
% % sigma = 500
% vKp2 = -real(Dds2./Nds2);
% vKi2 = w.*imag(Dds2./Nds2);
%         
% %recta 
% rKp2 = linspace(-0.1,0.25,1000);
% rKi2 = (1000*rKp2)+1000*(Dds2n/Nds2n);
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#34a0a4';
% color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp2,vKi2,'Color',color3)
% plot(rKp2,rKi2,'Color',color3,'HandleVisibility','off')
% 
% % sigma = 1000
% vKp3 = -real(Dds3./Nds3);
% vKi3 = w.*imag(Dds3./Nds3);
%         
% %recta 
% rKp3 = linspace(-0.1,0.25,1000);
% rKi3 = (1000*rKp3)+1000*(Dds3n/Nds3n);
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#52b69a';
% color4 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp3,vKi3,'Color',color4)
% plot(rKp3,rKi3,'Color',color4,'HandleVisibility','off')
% 
% % sigma = 2000
% vKp4 = -real(Dds4./Nds4);
% vKi4 = w.*imag(Dds4./Nds4);
%         
% %recta 
% rKp4 = linspace(-0.1,0.25,1000);
% rKi4 = (1000*rKp4)+1000*(Dds4n/Nds4n);
% 
% str = '#76c893';
% color5 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp4,vKi4,'Color',color5)
% plot(rKp4,rKi4,'Color',color5,'HandleVisibility','off')
% 
% % sigma = 3000
% vKp5 = -real(Dds5./Nds5);
% vKi5 = w.*imag(Dds5./Nds5);
%         
% %recta 
% rKp5 = linspace(-0.1,0.25,1000);
% rKi5 = (1000*rKp5)+1000*(Dds5n/Nds5n);
% 
% str = '#99d98c';
% color6 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp5,vKi5,'Color',color6)
% plot(rKp5,rKi5,'Color',color6,'HandleVisibility','off')
% 
% % sigma = 6000
% vKp6 = -real(Dds6./Nds6);
% vKi6 = w.*imag(Dds6./Nds6);
%         
% %recta 
% rKp6 = linspace(-0.1,0.25,1000);
% rKi6 = (1000*rKp6)+1000*(Dds6n/Nds6n);
% 
% str = '#b5e48c';
% color7 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(vKp6,vKi6,'Color',color7)
% plot(rKp6,rKi6,'Color',color7,'HandleVisibility','off')
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#f72585';
% color8 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% % punto de seleccion
% plot(0.213453038674033,1681.51315789474,'p','Color',color8)
% 
% legend('$\sigma=0$','$\sigma<100$','$\sigma<500$','$\sigma<1000$',...
%     '$\sigma<2000$','$\sigma<3000$','$\sigma<6000$','$(k_p^c,k_i^c)$', 'interpreter', 'latex','Location', 'northwest')
% 
% axis([-0.02 0.25 -1 1800]);
% xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
% ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
% title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');
% 
% % while (1)
% %     figure(1)
% %     [kp,ki]=ginput(1)
% 
% kp = 0.213453038674033;
% ki = 1681.51315789474;
%     
%     c3 = d12;
%     c2 = (kp*n11 + d11);
%     c1 = (ki*n11 + kp*n10 + d10);
%     c0 = ki*n10;
%     
%     n1 = (kp*Vo)/L;
%     n2 = ((2*kp*dp*IL)+(ki*Vo*C))/(L*C);
%     n3 = (2*ki*dp*IL)/(L*C);
%     
%     Pcl=[c3 c2 c1 c0];
%     roots(Pcl)
%     
%     Zcl=[n1 n2 n3];
%     roots(Zcl)
%     
%     Cc=pid(kp,ki);
%     sys=series(Cc,G);
%     H=1;
%     OL=G*Cc;
%     Mc=feedback(sys,H);
%     
%     figure(2)
%     step(Mc)
%     
%     % funcion en lazo cerrado analitica
%     nLC = [(kp*Vo)/L ((2*dp*IL*kp)+(Vo*C*ki))/(L*C) (2*dp*IL*ki)/(L*C)];
%     dLC = [1 ((kp*Vo*R*C)+L)/(R*L*C) ...
%            ((2*dp*IL*kp)+(dp^2)+(Vo*C*ki))/(L*C) ...
%            (2*dp*IL*ki)/(L*C)];
%        
%     GCL = tf(nLC,dLC);
%     
%     figure(3)
%     pzplot(GCL)
%     hold on
% 
% % end

