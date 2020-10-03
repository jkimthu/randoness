%%  Plasticity in growth rate across metabolic state, M



%  Goal: plot growth rate in high and low nutrient for M1 (low-adapted),
%        M2 (high-adapted) and intermediate M


%  Strategy:
%
%  Part 1:
%     1. import data from M1_low, M2_high and Mi_intermediate
%     2. plot!


%  Last edit: jen, 2020 Oct 2
%  commit: plot plasticity in growth rate across metabolic state


% Okie, go go let's go!


clc
clear

% 1. import data from M1_low, M2_high and Mi_intermediate
%    coordinates: (mu_high, mu_low) 


Mi = [1.49, 0.82;
      1.70, 0.66;
      1.88, 0.70;
      2.01, 0.57;
      2.12, 0.46;
      2.06, 0.53];

  
% 2. plot
color_m1 = rgb('Indigo');
color_m2 = rgb('FireBrick');
color_mi = rgb('DarkGray');

figure(1)
plot(1.36, 1.05,'o','MarkerFaceColor',color_m1,'MarkerEdgeColor',color_m1,'MarkerSize',10) % M1
hold on
plot(2.83, 0.11,'o','MarkerFaceColor',color_m2,'MarkerEdgeColor',color_m2,'MarkerSize',10) % M2
hold on
plot(Mi(:,1),Mi(:,2),'o','MarkerFaceColor',color_mi,'MarkerEdgeColor',color_mi,'MarkerSize',10) % Mi

axis([1,3,0,1.2])
ylabel('low growth rate (1/h)')
xlabel('high growth rate (1/h)')
title('plasticity in growth rate by metabolic state')
legend('M_1','M_2','M_i')

