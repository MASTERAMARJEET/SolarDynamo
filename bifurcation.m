close all;
clear all;

T0 = 2;
T1 = 0.5;
%%%%%%%%%%%%%%%%%%%%%%
lags = [T1 T0];

tspan = [0 1000];
sol = dde23(@ddefun, lags, @history, tspan);

figure(1)
syms t
Bphi = sol.y(1,:);
Bphimax = max(Bphi);
size(Bphimax)
deriv = diff(Bphi);
plot(sol.x(1,:),sol.y(1,:))

figure(2)
omega_by_L = -0.34;
tau = 15;
Ndmean = -13;
alphamean = Ndmean/(omega_by_L*(tau^2));
percentage = 0.30;
global alpha
alpha = alphamean*(1+percentage+(2*percentage*rand(1,1)));
Nd = alpha*omega_by_L*(tau^2);
size(Nd)
%scatter(abs(Nd),Bphimax);

function s = history(t)
  Bmin = 1;
  Bmax = 7;
  s = ones(2,1);
  s(1,1) = (Bmax + Bmin)/2;
  s(2,1) = (Bmax + Bmin)/2;
end

function dydt = ddefun(t,y,Z)
  omega_by_L = -0.34;
  tau = 15;
  Ndmean = -13;
  alphamean = Ndmean/(omega_by_L*(tau^2));
  percentage = 0.30;
  epsilon_lim = sqrt(3+percentage^2) * 0.17 * alphamean;
  ylag1 = Z(:,1);
  ylag2 = Z(:,2);
  B = y(1);
  global alpha
  %disp(t)

  dydt = [(omega_by_L*ylag2(2))-(B/tau); 
          (alpha*quenching(ylag1(1))*ylag1(1))-(y(2)/tau)-epsilon_lim+2*epsilon_lim*rand(1,1)];
end

function f = quenching(B)
    Bmin = 1;
    Bmax = 7;
    f = (1 + erf(B^2 - Bmin^2)) * (1 - erf(B^2 - Bmax^2)) * 0.25;
end