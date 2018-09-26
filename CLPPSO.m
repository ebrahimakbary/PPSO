clear all
clear all

 clc

disp('CLPPSO');

CostFunction=@(x) Cost(x);
nvar=30;
D=nvar;
xmin=-100;
xmax=-xmin;

dx=xmax-xmin;

vmax=0.5*dx;
npop=40;
maxit=100000/npop;


        ai = zeros(npop, nvar);
        f_pbest = 1 : npop;
        f_pbest = repmat(f_pbest', 1, nvar);
          t = 0 : 1 / (npop - 1) : 1;
        t = 5 .* t;
        Pc = 0.0 + (0.5 - 0.0) .* (exp(t) - exp(t(1))) ./ (exp(t(npop)) - exp(t(1)));
        m = 0 .* ones(npop, 1);
for it=1:maxit
    if it==1
        gbestcost(1)=inf;

        for i=1:npop
         velocity(i,:)=zeros(1,nvar);
            delta(i)=unifrnd(0,2*pi);
           position(i,:)=xmin+(xmax-xmin)*rand(1,nvar);
         cost(i)=CostFunction(position(i,:));%%,M ,o ,A, a, b, A
            MM(i)=7;
            FF(i)=0;
           pbest(i,:)=position(i,:);
            pbestcost(i)=cost(i);
            
            if pbestcost(i)<gbestcost(it)
                gbest=pbest(i,1:nvar);
                gbestcost(it)=pbestcost(i);
            end
        end
                for k = 1 : npop

            ar = randperm(nvar);
            ai(k, ar(1 : m(k))) = 1;
            fi1 = ceil(npop * rand(1, nvar));
            fi2 = ceil(npop * rand(1, nvar));
            fi = (pbestcost(fi1) < pbestcost(fi2)) .* fi1 + (pbestcost(fi1) >= pbestcost(fi2)) .* fi2;
            bi = ceil(rand(1, nvar) - 1 + Pc(k));

            if bi == zeros(1, nvar),
                rc = randperm(nvar);
                bi(rc(1)) = 1;
            end

            f_pbest(k, :) = bi .* fi + (1 - bi) .* f_pbest(k, :);

        end
    else

        gbestcost(it)=gbestcost(it-1);

        for i=1:npop
if FF(i)>=2
    FF(i)=0;
    ai(i, :) = zeros(1, nvar);
                    f_pbest(i, :) = i .* ones(1,nvar);
                    ar = randperm(nvar);
                    ai(i, ar(1 : m(i))) = 1;
                    fi1 = ceil(npop * rand(1, nvar));
                    fi2 = ceil(npop * rand(1, nvar));

fi = (pbestcost(fi1) < pbestcost(fi2)) .* fi1 + (pbestcost(fi1) >= pbestcost(fi2)) .* fi2;
                    bi = ceil(rand(1, nvar) - 1 + Pc(i));

                    if bi == zeros(1, nvar),
                        rc = randperm(nvar);
                        bi(rc(1)) = 1;
                    end

                    f_pbest(i, :) = bi .* fi + (1 - bi) .* f_pbest(i, :);
end

aa=2*(sin(delta(i)));
bb=2*(cos(delta(i)));
ee=abs(cos(delta(i)))^aa;
tt=abs(sin(delta(i)))^bb;

for j1=1:nvar
   velocity(i,j1)=((ee/(npop))*velocity(i,j1))+(ee)*(pbest(f_pbest(i, j1),j1)-position(i,j1));
end
         
           velocity(i,:)=min(max(velocity(i,:),-vmax),vmax);
       

            position(i,:)=position(i,:)+velocity(i,:);

            
           position(i,:)=min(max(position(i,:),xmin),xmax);
            
            cost(i)=CostFunction(position(i,:));
%             
  delta(i)=delta(i)+(abs(aa+bb)*(2*pi));


  
     vmax=(abs(cos(delta(i)))^2)*dx;
     
if cost(i)<pbestcost(i)
 FF(i)=0;
else
FF(i)=FF(i)+1;
end
            if cost(i)<pbestcost(i)
                pbest(i,:)=position(i,:);
                pbestcost(i)=cost(i);
            else

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

        semilogy((it),(log(gbestcost(it))),'-or','linewidth',1.4,'Markersize',3);
        
        xlabel('\fontsize{12}\bf Iteration');
        ylabel('\fontsize{12}\bf Best value');
        legend('\fontsize{10}\bf PPSO',1);
       Cost_Rsult=gbestcost(end)
%% toc