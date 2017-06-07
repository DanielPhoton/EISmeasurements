function [Voc Jsc FF Eff scanrate] = SCparameters(V,T,J)
    Vrate = gradient(V);
    Vrate(isnan(Vrate(:,1)),:)=[];
    medianVrate = median(abs(Vrate));
    Vrate = gradient(V);
    a=abs(gradient(V)./gradient(T));
    a(isnan(a(:,1)),:)=[];
    scanrate=median(a);
    Vmax=max(V); Vmin=min(V);
    Voc= [ 0; 0 ]; Jsc=[0;0]; Mpp=[0;0]; FF=[0;0]; Eff=[0;0];
    ju=0;vu=0;jd=0;vd=0;ffu=0;ffd=0;
    power=V.*J;
    for i = 2:length(V)
        if Vrate(i) > 0.5*medianVrate
            if V(i) > 0 && V(i-1) <= 0
                ju=ju+1;
                Jsc(1,ju)=J(i);% add to Jsc up
                juind(ju)=i;
            elseif J(i) >0 && J(i-1) <=0
                vu=vu+1;
                Voc(1,vu)=V(i);% add to Voc up
                vuind(vu)=i;
            end
        elseif Vrate(i) < -0.5*medianVrate
            if V(i) < 0 && V(i-1) >=0
                jd=jd+1;
                Jsc(2,jd)=J(i);% add to Jsc down
                jdind(jd)=i;
            elseif J(i) <0 && J(i-1) >=0
                vd=vd+1;
                Voc(2,vd)=V(i);% add Voc up
                vdind(vd)=i;
            end
        else
            % do nothing
        end
    end
    try
        for i= 1:length(vuind)
            if juind(i) > vuind(i) && vdind(i) < jdind(i)
                Mpp(1,i)= min(power(juind(i)- juind(1)+1:vuind(i)));
                Mpp(2,i)= min(power(vdind(i):jdind(i)));
            elseif vdind(i) > jdind(i) && juind(i) < vuind(i)
                Mpp(1,i)= min(power(juind(i):vuind(i)));
                Mpp(2,i)= min(power(vdind(i)-vdind(1)+1:jdind(i)));
            else
                Mpp(1,i)= min(power(juind(i):vuind(i)));
                Mpp(2,i)= min(power(vdind(i):jdind(i)));
            end
        end
        FF=Mpp./(Voc.*Jsc);
        Eff = Voc.*(-Jsc).*FF; % correct for area and Amps to mAmps converstion
    catch ME
        disp('Error! setting FF and Eff = 0')
        FF= zeros(2,length(vuind));
        Eff = zeros(2,length(vuind));
    end
end