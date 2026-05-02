
clc
clear all

c = [1 2 3 0 0];
b = [20; 30];
A = [1 2 0 1 0; 3 0 4 0 1];

m = size(A,1);
n = size(A,2);

bv_index = n-m+1:n;
Y = [A b];


for s=1:50    
    cb = c(bv_index);
    Xb = Y(:,end);
    
    zj = cb*Xb;
    zjcj = cb*Y(:,1:n) - c;
    Table=[zjcj zj;Y]

    if(zjcj>=0)
        disp(Xb)
        disp(zj)
        break
    else 
        [a EV]=min(zjcj);
        if(Y(:,EV)<=0)
            disp('Unbounded')
            break
        else
            for j=1:m
                if(Y(j,EV)>0)
                    ratio(j)=Xb(j)/Y(j,EV);
                else
                    ratio(j)=inf;
                end
            end
        end
        [k, LV]=min(ratio);
        bv_index(LV)=EV;
    end
    pivot=Y(LV,EV);
    %pivot row is row/pivot element ie divide by 4 
    Y(LV,:)=Y(LV,:)/pivot;
    for i=1:m
        if(i~=LV)
            Y(i,:)=Y(i,:)-Y(i,EV)*Y(LV,:);
        end
    end
end

           
% big m
clc;
clear;

M = 1e5;

C = [4 3 0 -M -M];

A = [1 2 -1 1 0;
     3 1  0 0 1];

b = [6; 9];

[m,n] = size(A);

Y = [A b];

bv_index = [4 5];

for s = 1:50
    
    Xb = Y(:,end);
    Cb = C(bv_index);
    
    Z = Cb * Xb;
    
    ZjCj = Cb * Y(:,1:n) - C;
    
    Table = [ZjCj Z; Y];
    disp(Table);
    
    if all(ZjCj >= 0)
        sol = zeros(1,n);
        sol(bv_index) = Xb;
        disp(sol);
        disp(Z);
        break;
        
    else
        [~,Ev] = min(ZjCj);
        
        if all(Y(:,Ev) <= 0)
            disp('unbounded');
            break;
        else
            for j = 1:m
                if Y(j,Ev) > 0
                    ratio(j) = Xb(j)/Y(j,Ev);
                else
                    ratio(j) = inf;
                end
            end
        end
        
        [~,Lv] = min(ratio);
        bv_index(Lv) = Ev;
        
        pivot = Y(Lv,Ev);
        Y(Lv,:) = Y(Lv,:)/pivot;
        
        for i = 1:m
            if i ~= Lv
                Y(i,:) = Y(i,:) - Y(i,Ev)*Y(Lv,:);
            end
        end
    end
end


    %least cost method
clc;
clear;
cost=[11 13 17 14;
      16 18 14 10;
      21 24 13 10];
supply=[10 5 9];
demand=[8 7 15 4];
[m,n]=size(cost);
S=sum(supply);
D=sum(demand);


if(S==D)
    disp("Problem is balanced")
elseif(S<D)
    cost(end+1,:)=zeros(1,n);
    supply(end+1)=D-S;
else 
    cost(:,end+1)=zeros(m,1);
    demand(end+1)=S-D;
end
disp("Balanced TP")
Table=[cost supply'; demand sum(supply)]
[m,n]=size(cost);
X=zeros(m,n);
Icost=cost;
%allocation
while(any(supply)~=0||any(demand)~=0)
    min_cost=min(cost(:));
    [r c]=find(cost==min_cost);
    y=min(supply(r),demand(c));
    [aloc,index]=max(y);
    rr=r(index);
    cc=c(index);
    X(rr,cc)=aloc;
    supply(rr) = supply(rr) - aloc;
    demand(cc) = demand(cc) - aloc;
    cost(rr, cc) = inf; 
end
cost_ec=X.*Icost;
final_cost=sum(cost_ec(:))
disp(X)
disp(final_cost)


format short 
clear all    
clc          

% steepest graident 
syms x1 x2
f1 = x1-x2+2*x1^2+2*x1*x2+x2^2;
fx = inline(f1); 
fobj = @(x) fx(x(:,1),x(:,2)); 


grad = gradient(f1);
G = inline(grad); 
gradx = @(x) G(x(:,1),x(:,2)); 


H1 = hessian(f1);
Hx = inline(H1); 

x0 = [1 1]; 
maxiter = 4; 
tol = 10^(-3); 
iter = 0; 
X = []; 

while norm(gradx(x0)) > tol && iter < maxiter
    X = [X; x0]; 
    S = -gradx(x0); 
    H = Hx(x0); 
    lambda = S'*S ./ (S'*H*S); 
    Xnew = x0 + lambda.*S'; 
    x0 = Xnew; 
    iter = iter + 1; 
end

%% Phase 5- Print the solution
fprintf('Optimal Solution X = %f, %f\n',x0(1),x0(2))
fprintf('Optimal Value f(x) = %f \n',fobj(x0))
