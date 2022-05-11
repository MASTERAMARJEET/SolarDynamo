close all;
clear all;

tau = 15;
Bmin = 1;
Bmax = 7;
T0 = 2;
T1 = 0.5;
alpha = 0.17;
omega_by_L = -2*alpha;
Nd = alpha*omega_by_L*(tau^2);
epsilon_lim = sqrt(3) * 0.17 * alpha;

%%%%%%%%%%%%%%%%%%%%%%
lags = [T1 T0];

tspan = [0 5000];
sol = dde23(@ddefun, lags, @history, tspan);

figure(1)
plot(sol.x,sol.y)
xlabel('Time t');
ylabel('Solution y');

figure(2)
syms t
Bphi = sol.y(1,:);
deriv = diff(Bphi);
plot(Bphi(1,2:end), deriv);
title(["Bphi vs dBphi/dt phase plot (with noise) with dynamo number ",num2str(Nd)]);
%title(["Bphi vs dBphi/dt phase plot (without noise) with dynamo number ",num2str(Nd)]);
xlabel("Bphi");
ylabel("dBphi/dt");
hold on
plot(Bphi(1,2), deriv(1,1),'r*');
plot(Bphi(1,end), deriv(1,end),'g*');


function s = history(t)
  Bmin = 1;
  Bmax = 7;
  s = ones(2,1);
  s(1,1) = (Bmax + Bmin)/2;
  s(2,1) = (Bmax + Bmin)/2;
end

function dydt = ddefun(t,y,Z)
  alpha = 0.17;
  omega_by_L = -2*alpha;
  tau = 15;
  epsilon_lim = sqrt(3) * 0.17 * alpha;
  ylag1 = Z(:,1);
  ylag2 = Z(:,2);
  B = y(1);
  %dydt = [(omega_by_L*ylag2(2))-(B/tau); 
  %        (alpha*quenching(ylag1(1))*ylag1(1))-(y(2)/tau)]; %%for without
  %        noise

  dydt = [(omega_by_L*ylag2(2))-(B/tau); 
          (alpha*quenching(ylag1(1))*ylag1(1))-(y(2)/tau)-epsilon_lim+2*epsilon_lim*rand(1,1)];
end

function f = quenching(B)
    Bmin = 1;
    Bmax = 7;
    f = (1 + erf(B^2 - Bmin^2)) * (1 - erf(B^2 - Bmax^2)) * 0.25;
end