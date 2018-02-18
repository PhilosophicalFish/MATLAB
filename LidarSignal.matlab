%Lidar Signal Simulation%

%Creator: Maurice E. D. Roots | meroots015@gmail.com 
%Updated: Feb. 11th, 2018 [02/11/2018]

%This script requires Mie.m and Mie_abcd.m be in the same folder as this
%one. Otherwise the "Mie Function" section can not be run.

clear all; %House Keeping%
%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%

% ************* molecular extinctiont ************* %
str='C:\Users\meroo\Desktop\Research\Atmo\Mie';

cd(str);

[FileName,PathName] = uigetfile({'*.txt';'*.*'},'Molecular Extinction Data File Selector');
full_name = [PathName,FileName];
    
       ME = importdata(full_name); 
       me355 = ME(:,1);
       me532 = ME(:,2);
       me1064 = ME(:,3);
       height = ME(:,4);
       
% ************* molecular backscatter ************* %       
cd(str);

[FileName,PathName] = uigetfile({'*.txt';'*.*'},'Molecular Backscattering Data File Selector');
full_name = [PathName,FileName];

       MB = importdata(full_name); 
       mb355 = MB(:,1);
       mb532 = MB(:,2);
       mb1064 = MB(:,3);
       height1 = MB(:,4);
%%%%%%%%%%%%%%%%%%%%% VALUES %%%%%%%%%%%%%%%%%%%%%

rm = [0.0285, 0.0285, 0.0285];    %Size Distribution Mode Radius
    %  [Clean-Continental, Urban, Martime]

sd = [0.350, 0.350, 0.400];    %Size Distribution Standard Deviation
    %  [Clean-Continental, Urban, Martime]

z = [0.0075:0.0075:9];    %Altitude (km)

r = [0.01:0.01:10];    %Aerosol Size (um)
    
C = (3000*(1e-6))./(exp((z./1.2)));   %Aerosol Density Number

%format long   %to double cheack indexies of refraction
maritime = [complex(1.51, 4e-7); complex(1.5, 1e-8); complex(1.47, 2e-4)];
cleanContinental = [complex(1.53, 1.5e-3); complex(1.53, 1.5e-3); complex(1.52, 1.5e-3)];
urban = [complex(1.75, 0.47); complex(1.75, 0.44); complex(1.75, 0.44)];

M = [cleanContinental, urban, maritime];

%%%%%%%%%%%%%%%%% Mie Function %%%%%%%%%%%%%%%%

for i = 1:1000
    x1(i) = (2*pi*r(i)/(1e6))/(355/(1e9));     %355nm wavelength
    result1 = Mie(M(1), x1(i));     %results of mie with 100 points
    qext1(i) = result1(4);     %Efficiencies for extinction
    qsca1(i) = result1(5);     %scattering
    qabs1(i) = result1(6);     %absorbtion
    qb1(i) = result1(7);     %backscattering
    asy1(i) = result1(8);     %asymetry parameter (asy = <cos(theta)>)
    qratio1(i) = result1(9);     %(qratio = (backscatering qb)/(scattering qsca)
end

%%%%%%%%%%%%% Aerosol Data %%%%%%%%%%%%%

 for i=1:1200
     for j=1:1000
         sae1(j)=qext1(j)*C(i)/((r(j)/(1e6))*log(10)*(2*pi)^0.5)*exp(-1*(log((r(j)/(1e6))/(rm(1)/(1e6))))^2/(2*sd(1)^2))*0.0075;
         sab1(j)=qb1(j)*C(i)/((r(j)/(1e6))*log(10)*(2*pi)^0.5)*exp(-1*(log((r(j)/(1e6))/(rm(1)/(1e6))))^2/(2*sd(1)^2))*0.0075;
     end
     
     ae1(i)=sum(sae1);
     ab1(i)=sum(sab1);
 end

 ae1(200:1200) = 0;
 ab1(200:1200) = 0;
 
%%%%%%%%%%%%% Lidar Signal %%%%%%%%%%%%%

lidars1 =(ab1(1:800)+mb355').*exp(-2*cumsum(ae1(1:800)+me355')*0.0075);

%%%%%%%%%%%%% Gradient Method %%%%%%%%%%%%%

fz = lidars1.*((z(1:800)).^2);    %f(z) = P(z)z^2

gm = gradient(fz);

%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%
figure;    %New Figure Window
plot(log(lidars1), z(1:800));
title('Lidar Signal Simulation 355nm Clean-Continental');    %Title of Graph
ylabel('Altitude');    %X-axis Label
xlabel('Lidar Signal (log)');    %Y-axis Label

figure
plot(gm, z(1:800));
pbl = find(min(gm) == gm);
xmin = gm(pbl);
ymin = z(pbl);
strmin = ['PBL Height = ',num2str(ymin), 'km'];
text(xmin,ymin+1,strmin,'HorizontalAlignment','center');
%figure;plot(r,qb1);
%figure;plot(z(1:1000), sae1);
%figure;plot(z(1:1000), sab1);