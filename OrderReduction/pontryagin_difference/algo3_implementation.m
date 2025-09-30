 
 % load('pontry_example_reducable_example1.mat')
 % load('pontry_example_nonreducable_example1.mat')
 % load('pontry_example_nonreducable_example2.mat')

set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 13)
set(0,'defaultTextFontSize' , 13)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')
fig = figure('Position',[100 100 400 315]);
set(groot,'defaultLegendInterpreter','latex');

 G = [1, 0.5, -0.5; 
      0, 0.9, 0.9];
c=zeros(2,1);
Z1=zono(G,c);

G = [0.2  0.1 -0.1  0.15 -0.15;  
      0.1 -0.2  0.15 0.1  -0.1];

c=zeros(2,1);
Z2=zono(G,c);

% Zi = new_pontry_approach(Z1,Z2);
Zi = pontryDiff(Z1,Z2);

 [v,~]=plot(Zi,'r',0.5);

 disp(Zi);
 tic
 [isZono, generators] = isZonotope_ds6(v);

 result_zono=zono(generators,[0;0]);
toc
 disp(result_zono);

 hold on;

 plot(Z1,'b','0.3');
 plot(Z2,'g',0.3);
 legend("$\mathcal{Z}_3=\mathcal{Z}_3'$",'$\mathcal{Z}_1$','$\mathcal{Z}_2$', 'Location','northeast','Interpreter','latex');
yticks(-2:1:2) 
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
set(gcf, 'Color', 'w');
