clear all
clear all

    clc


CostFunction=@(x) Cost(x);
nvar=30;
D=nvar;
xmin=-100;
xmax=-xmin;

dx=xmax-xmin;

vmax=0.5*dx;

    npop=60;
maxit=100000/npop;




for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for i=1:npop
         velocity(i,:)=zeros(1,nvar);
            delta(i)=unifrnd(0,2*pi);
           position(i,:)=xmin+(xmax-xmin)*rand(1,nvar);
         cost(i)=CostFunction(position(i,:));
            
           pbest(i,:)=position(i,:);
            pbestcost(i)=cost(i);
            
            if pbestcost(i)<gbestcost(it)
                gbest=pbest(i,1:nvar);
                gbestcost(it)=pbestcost(i);
            end
        end
    else

        gbestcost(it)=gbestcost(it-1);

        for i=1:npop

 
aa=2*(sin(delta(i)));
bb=2*(cos(delta(i)));
ee=abs(cos(delta(i)))^aa;
tt=abs(sin(delta(i)))^bb;
                         

 KK=(ee)*(pbest(i,:)-position(i,:)) +(tt)*(gbest-position(i,:));

velocity(i,:)=((ee/(i))*velocity(i,:))+KK;

           
           velocity(i,:)=min(max(velocity(i,:),-vmax),vmax);

            position(i,:)=position(i,:)+velocity(i,:);


           position(i,:)=min(max(position(i,:),xmin),xmax);
            
            cost(i)=CostFunction(position(i,:));
            
  delta(i)=delta(i)+(abs(aa+bb)*(2*pi));


     vmax=(abs(cos(delta(i)))^2)*dx;
            if cost(i)<pbestcost(i)
                pbest(i,:)=position(i,:);
                pbestcost(i)=cost(i);

                if pbestcost(i)<gbestcost(it)
                    gbest=pbest(i,:);
                    gbestcost(it)=pbestcost(i);
                end

                end
        end
    end
    
    disp(['Iteration ' num2str(it) ':   Best Cost = ' num2str(gbestcost(it))]);
    
  
end
hold on;
figure(1)
 it=1:50:maxit;

        semilogy((it),(log(gbestcost(it))),'-ob','linewidth',1.4,'Markersize',3);
        
        xlabel('\fontsize{12}\bf Iteration');
        ylabel('\fontsize{12}\bf Best value');
        legend('\fontsize{10}\bf PPSO',1);
        Cost_Rsult=gbestcost(end);

