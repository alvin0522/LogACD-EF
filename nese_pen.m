%% Modified at 4/15/2017
%1. NLEQN with penalty
%2. Assume true model is (1,0), (2,0) and (1,1); fit a (10,0) model.
%3. Test for June 3rd BAC using (10,0). It works for all cases.  
%4. It should be g_theta - penalty = 0.
function F=nese3(para)
 
fid = fopen('data_info1.csv');
out = textscan(fid,'%f','Delimiter',',');
fclose(fid);
sim=out{1};

fid = fopen('initial1.csv');
out = textscan(fid,'%f','Delimiter',',');
fclose(fid);
init=out{1};

simN=length(sim);
paraN=length(para);
p=sim(simN-1);
q=sim(simN);
sim=sim(1:(simN-2));
omg = para(1);
alp = para(2:(1+p));
beta = para((p+2):paraN-1);
if q == 0 
    beta = 0;
else 
    beta = para((p+2):paraN-1);
end;
%z = para(paraN);
n = length(sim);
pdim = p+q;
psih = ones(1,n);
gtheta = zeros(paraN,1);

fid1 = fopen('moments1.csv');
%%fid1 = fopen('C:\Users\jzou\Dropbox\Durations April 2015\Code(R and Matlab) log ACD1\dist.csv');
out1 = textscan(fid1,'%f','Delimiter',',');
fclose(fid1);
dist=out1{1};
mue=dist(1);
vare=dist(2);
skewe=dist(3);
kurte=dist(4);
        
z = para(paraN);

 theta0 = abs(init);
 theta = theta0.^2./abs(theta0);
%theta = abs(para);
% sgn = sign(para);
fid = fopen('lambda1.csv');
out2 = textscan(fid,'%f','Delimiter',',');
fclose(fid);
lm = out2{1};
lambda = lm(1);
a = lm(2);
penalty = zeros((pdim+1),1);

  for i= 2: (length(theta))
   term1 = theta(i) <= lambda;
   term2 = max((a*lambda - theta(i)),0)/((a-1)*lambda);
   penalty(i) = (term1 + (theta(i) > lambda)*term2)*lambda/theta0(i)*para(i);
   end

%  for i= 2: (length(theta)-1)
%   term1 = theta(i) <= lambda;
%   term2 = max((a*lambda - theta(i)),0)/((a-1)*lambda);
%   penalty(i) = (term1 + (theta(i) > lambda)*term2)*lambda*sign(para(i));
%   end
            
mxpq = max(p,q);
t0 = mxpq + 1;
psih(1:mxpq)=omg;  % omega 
for t=t0:n
    alpde=log(flip(sim(t-p:t-1)));
    betade=flip(psih(t-q:t-1));
    psih(t) = omg+sum(transpose(alp)*alpde)+sum(transpose(beta)*transpose(betade));
    %psih(t) = omg+sum(alp*alpde)+sum(beta*transpose(betade));
    %% mu(t), sigsq(t), gamma(t), kappa(t)	
    mu = exp(psih(t));
    sigsq = vare*exp(2*psih(t))/mue^2;
    gama = skewe*exp(3*psih(t))/mue^3;
    kappa = kurte*exp(4*psih(t))/mue^4;
    
    %% Define derivatives of psih(t) wrt theta: pdim*1 vector  						
    %% First derivative	

    derpsi = transpose([1,transpose(alpde),betade]);
    
    % First Derivatives of mu(t),sigsq(t), gamma(t) and kappa(t)												
    dermu = exp(psih(t))*derpsi;							
    dersigsq = 2*vare*exp(2*psih(t))*derpsi/mue^2;
    %dergamma = 3*skewe*exp(3*psih(t))*derpsi;
    %derkappa = 4*kurte*exp(4*psih(t))*derpsi;
    
    %% Compute m(t) and M(t)							
    m = sim(t)-mu;       							
    qm = m^2-sigsq;
    
    %% Compute Quadratic variations of m(t) and M(t) and 
    %% covariance of (m(t), M(t))							
    vm = sigsq;
    vqm = kappa - vm^2;
    vmqm = gama;
    
    %% Define rho^2(t) 															
    termr = 1-((vmqm^2)/(vm*vqm));								
    rho = 1/termr;
    
    %% Define eta 								
    eta = vmqm/(vm*vqm);	
    
    %% Define Derivatives of m(t) and M(t)     					
    derm = -dermu;
    derqm = 2*m*derm - dersigsq;
    
    %% Define vectors astr and bstr						
    astr = rho*(-dermu/vm + dersigsq*eta);							
    bstr = rho*(dermu*eta - dersigsq/vqm);
    
    %% Constrain Transformation %%
    fz = z - 1 + sum(alp) + sum(beta);
    
    %% compute thehat[t]
    gtheta(1:(pdim+1)) = astr*m + bstr*qm + gtheta(1:(pdim+1)) - penalty;
    gtheta(pdim+2) = fz;
end

F=gtheta;
    
    
    
    
 
 



