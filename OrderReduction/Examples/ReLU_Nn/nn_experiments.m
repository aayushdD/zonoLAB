
leaves = getLeaves(NN);
nLeaves = size(leaves,2);
hold on
equal=[];
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times') 
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 16)
set(0,'defaultTextFontSize' , 16)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

for i=1:nLeaves
    Zi = conZono(NN.Gc,NN.c+NN.Gb*leaves(:,i),NN.Ac,NN.b-NN.Ab*leaves(:,i));
    Zi=projection(Zi,[1,2]);
    % Zj=lp_bounds_rem_function(Zi);
    Zj=reduce(Zi);
    % Zj=Iterative_bounds_rem_function(Zi);

    %To verify exactness, contain_check_v3 must be applied to both Zj and
    %Zi as the including set and the included set. If both contain each
    %other, equality can be considered.

    equal(end+1,:)=contain_check_v3(Zj,Zi,1e-3,0);

    if i==833
         plot(Zi,'y',1);
         drawnow
         continue
    end
    plot(Zi,'b',0.3);
    hold on;
    plot(Zj,'r',0.3);
    drawnow
end

