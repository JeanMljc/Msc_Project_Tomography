function[x]=gradC(R_carre,g,N)
    x=[zeros((N^2),1)];
    d_pre=g-R_carre*x(:,1);
    b_pre=-d_pre(:,1);
    N=1000;

    %%while ((r_pre')*r_pre)>100
    for k=(1:N)
        alpha=(-(b_pre')*d_pre)/((d_pre')*R_carre*d_pre);
        %maj
        x_next=x(:,k)+alpha.*d_pre;
        b_next=R_carre*x_next-g;
        beta=(-(b_next')*R_carre*d_pre)/((d_pre')*R_carre*d_pre);
        d_next=b_next+beta.*d_pre;
    
        % if ((r_next')*r_next)>((r_pre')*r_pre)
        %   break
        % end
   
        %maj variable
        b_pre=b_next;
        d_pre=d_next;
        x=[x,x_next];
        %disp(norm(r_pre))
    end
end

