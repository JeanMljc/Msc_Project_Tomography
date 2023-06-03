function[f]=ART(R,g,N)
    %programme assez long 2-3 min 
    f=(1/N^2).*ones(N^2,1);
    k=1;
    lambda=0.1; %lambda optimal 0.1 mais image pas belle pour Ni=m
    nP=ceil(N*sqrt(2))+3;
    m=nP*180; %m=6600=6.10^3 si on veut les restes Ni>m
    Ni=m; %attention si Ni<m la reconstruction n'est pas complète => distorsion

    while k<Ni
        f_new=[];
        for j=(1:m)
            if (norm(R(j,:))>0) && (k==mod(j,m))
                corr=lambda*((g(j)-(R(j,:)*f(:,end)))*(R(j,:)'))/(norm(R(j,:))^2);
                f_new=f(:,end)+corr;
            end
        end
        f=[f,f_new];
        k=k+1;
    end