function[R]=ker_radon(N)
    theta=(0:179);
    N_phi=length(theta);

    nP=ceil(N*sqrt(2))+3; %+3 ou +4
    xp=(-(nP/2):1:(nP/2));
    R=zeros(nP*N_phi,N*N);

    for t=theta
        for p=xp
            for v=xp
                x=(N/2)+ceil(p*cos(((pi/180)*t)+pi/2)-v*sin(((pi/180)*t)+pi/2))-1;
                y=(N/2)+ceil(p*sin(((pi/180)*t)+pi/2)+v*cos(((pi/180)*t)+pi/2));
                if(x>0 & y>0 & x<=N & y<=N)
                    presproj=(nP*t)+p+(nP/2)+1;
                    pix=(y-1)*N+x;
                    R(presproj,pix)=1;
                end
            end
        end
    end
end
