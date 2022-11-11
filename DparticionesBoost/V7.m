% ========================================================================
% Fronteras de Estabilidad para convertidor CD-CD Boost
% Autor : Erick Moreno Negrete
% Version : 1.6
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
%                    Funciones de transferencia 
% ========================================================================

numC = [0 Vo/L (2*dp*IL)/(L*C)];
numV = [0 -IL/C (Vo*dp)/(L*C)];
den =  [1 1/(R*C) (dp^2)/(L*C)];

G0 = tf(numV,den);  % Voltaje/d
G1 = tf(numC,den);  % Corriente/d
G2 = tf(numV,numC); % Voltaje/Corriente

% elementos de diseño
a1 = -numV(2);
a2 = numV(3);
a3 = numC(2);
a4 = numC(3);

b1 = den(2);
b2 = den(3);

% ========================================================================
%                       Funcion lazo abierto G3
% ========================================================================

% ========================================================================
%                             Opcion #1

kp = 0.202513477088949;
ki = 703.894701986755;
%                               Raices
%            -5368.21202904105 +      9266.54604394599i
%            -5368.21202904105 -      9266.54604394599i
%            -2219.88921834814 +                     0i
%
%                               Ceros
%                         -8333.331666667
%                         -3475.79189348266
% ========================================================================


% ========================================================================
%                             Opcion #2

%kp = 0.213453038674033;
%ki = 1681.51315789474;
%                               Raices
%           -4637.91993256466 +      11173.8225981438i
%           -4637.91993256466 -      11173.8225981438i
%           -4155.28077176463 +                     0i
%
%                               Ceros
%                         -8333.33166666699
%                         -7877.67261754752
% ========================================================================





nLC = [-a1*kp (a2*kp)-(a1*ki) a2*ki];
dLC = [1 b1+(a3*kp) b2+(a4*kp)+(a3*ki) a4*ki];
   
G3ol = tf(nLC,dLC);

% elementos de diseño
n1 = nLC(1);
n2 = nLC(2);
n3 = nLC(3);

d1 = dLC(1);
d2 = dLC(2);
d3 = dLC(3);
d4 = dLC(4);

% ========================================================================
%     Caracterizacion de las fronteras de estabilidad lazo externo
% ========================================================================

% Fronteras de sigma

w = linspace(0,17000,10000);
sb = 1i.*w;

s1 = -1000+(1i.*w);
s1n = -1000;
    
s2 = -2000+(1i.*w);
s2n = -2000;
    
s3 = -3000+(1i.*w);
s3n = -3000;

%Polinomios en s=jw; fronteras de estabilidad
Nds0 = n1*sb.^2 + n2*sb + n3;
Dds0 = d1*sb.^3 + d2*sb.^2 + d3*sb + d4;

%Polinomios en s=-100+jw; 
Nds1 = n1*s1.^2 + n2*s1 + n3;
Dds1 = d1*s1.^3 + d2*s1.^2 + d3*s1 + d4;

%Polinomios en s=-100; 
Nds1n = n1*s1n.^2 + n2*s1n + n3;
Dds1n = d1*s1n.^3 + d2*s1n.^2 + d3*s1n + d4;

%Polinomios en s=-500+jw; 
Nds2 = n1*s2.^2 + n2*s2 + n3;
Dds2 = d1*s2.^3 + d2*s2.^2 + d3*s2 + d4;
    
%Polinomios en s=-500; 
Nds2n = n1*s2n.^2 + n2*s2n + n3;
Dds2n = d1*s2n.^3 + d2*s2n.^2 + d3*s2n + d4;

%Polinomios en s=-1000+jw; 
Nds3 = n1*s3.^2 + n2*s3 + n3;
Dds3 = d1*s3.^3 + d2*s3.^2 + d3*s3 + d4;
    
%Polinomios en s=-1000; 
Nds3n = n1*s3n.^2 + n2*s3n + n3;
Dds3n = d1*s3n.^3 + d2*s3n.^2 + d3*s3n + d4;


vKp0 = -real(Dds0./Nds0);
vKi0 = w.*imag(Dds0./Nds0);
rKp0 = linspace(-0.05,0.03,1000);
rKi0 = linspace(0,0,1000);
rKi0y = linspace(-1,320,1000);

figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(gca,'GridLineStyle','--')

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#2f3e46';
color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

% sigma = 0
plot(vKp0,vKi0,'Color',color1)
hold on
plot(rKp0,rKi0,'b--','HandleVisibility','off')
plot(rKi0,rKi0y,'b--','HandleVisibility','off')

% sigma = 100
vKp1 = -real(Dds1./Nds1);
vKi1 = w.*imag(Dds1./Nds1);

%recta 
rKp1 = linspace(-0.05,0.03,1000);
rKi1 = (1000*rKp1)+1000*(Dds1n/Nds1n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#354f52';
color2 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp1,vKi1,'--','Color',color2)
plot(rKp1,rKi1,'--','Color',color2,'HandleVisibility','off')

% sigma = 500
vKp2 = -real(Dds2./Nds2);
vKi2 = w.*imag(Dds2./Nds2);
        
%recta 
rKp2 = linspace(-0.05,0.03,1000);
rKi2 = (2000*rKp2)+2000*(Dds2n/Nds2n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#52796f';
color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp2,vKi2,':','Color',color3)
plot(rKp2,rKi2,':','Color',color3,'HandleVisibility','off')

% sigma = 1000
vKp3 = -real(Dds3./Nds3);
vKi3 = w.*imag(Dds3./Nds3);
        
%recta 
rKp3 = linspace(-0.05,0.03,1000);
rKi3 = (3000*rKp3)+3000*(Dds3n/Nds3n);

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#84a98c';
color3 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(vKp3,vKi3,'-.','Color',color3)
plot(rKp3,rKi3,'-.','Color',color3,'HandleVisibility','off')

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#f72585';
color4 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(0.5e-3,50,'p','Color',color4)

legend('$\sigma=0$','$\sigma<100$','$\sigma<500$','$\sigma<1000$',...
    '$(k_p^v,k_i^v)$', 'interpreter', 'latex','Location', 'Best')
axis([-0.04 0.04 0 250]);
xlabel('$$k_p$$','FontSize', 20 , 'interpreter', 'latex');
ylabel('$$k_i$$','FontSize', 20 ,  'interpreter', 'latex');
title('$\mathcal{D}-$particiones Convertidor Boost','FontSize', 12 ,'interpreter', 'latex');



while (1)
    figure(1)
    [kpv,kiv]=ginput(1)
    
    c4 = 1;
    c3 = (b1 + a3*kp - a1*kpv*kp);
    c2 = (b2 + a4*kp + a3*ki + a2*kpv*kp - a1*kpv*ki - a1*kiv*kp);
    c1 = (a4*ki + a2*kpv*ki + a2*kiv*kp - a1*kiv*ki);
    c0 = a2*kiv*ki;
    
    Pcl=[c4 c3 c2 c1 c0];
    roots(Pcl)
    
    z1 = -a1*kpv*kp;
    z2 = a2*kpv*kp - a1*kpv*ki - a1*kiv*kp;
    z3 = a2*kpv*ki + a2*kiv*kp - a1*kiv*ki;
    z4 = a2*kiv*ki;
    
    Zcl = [z1 z2 z3 z4];
    roots(Zcl)
    
    
%     
%     Cc=pid(kpv,kiv);
%     sys=series(Cc,G);
%     H=1;
%     OL=G*Cc
%     Mc=feedback(sys,H);
%     
%     figure(2)
%     step(Mc)
%     
%     funcion en lazo cerrado analitica
%     nLC = [(kp*Vo)/L ((2*dp*IL*kp)+(Vo*C*ki))/(L*C) (2*dp*IL*ki)/(L*C)];
%     dLC = [1 ((kp*Vo*R*C)+L)/(R*L*C) ...
%            ((2*dp*IL*kp)+(dp^2)+(Vo*C*ki))/(L*C) ...
%            (2*dp*IL*ki)/(L*C)];
%        
%     GCL = tf(nLC,dLC);

end


