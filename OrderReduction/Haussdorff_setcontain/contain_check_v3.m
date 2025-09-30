function [verdict] = contain_check_v3(A, B, d, tol)

new_B = B + zono(d*eye(B.n),zeros(B.n,1));
[Va,~] = plot(A,plotOptions('Display','off'));
[Vb,~] = plot(new_B,plotOptions('Display','off'));

[H,f,~,~] = vrep_to_Hrep_mpt3(Vb);

N = size(Va,1);
insideMask = (H*Va'-repmat(f,1,N))<=tol;
verdict= all(insideMask(:));
end


