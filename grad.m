function[x]=grad(R_carre,g,N)
    x=[zeros((N^2),1)];
    r=g-R_carre*x(:,1);
    while ((r')*r)>100
        gam=((r')*r)/((r')*R_carre*r);
        x_next=x(:,end)+gam.*r;
        
        %actualisation et stockage
        r=g-R_carre*x_next;
        x=[x,x_next];
    end
end

