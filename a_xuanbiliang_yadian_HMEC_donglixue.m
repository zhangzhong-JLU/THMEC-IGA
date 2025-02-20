tic 

addpath ../fem_util/
addpath ../C_files/  
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-util/


clear variables
close all
global sdof edof   Al Pe Pm Dt
global  matmtx matmtx1  matmtx2  matmtx3 matmtx4  matmtx5 matmtx6 matmtx7
global p q

suanlidata%L=0.02m,h=0.002m

refineCount   = 5; %
noGPs         = 5; % # of Gauss points surface integral along one direction
noGPs1        = 2; % # of GPs for line integral
E0  =7.6e10;  % Young's modulus
nu0 = 0.31;  % Poisson ratio
stressState ='PLANE_STRAIN';  
force = 1;   % disp control, force=1: force control
ubar = 0.4;
F=-100; 
if (refineCount) 
    hRefinement2d 
end

generateIGA2DMesh
plotMesh(controlPts,weights,uKnot,vKnot,p,q,50,'r-','try.eps');



%%-----------------------------------------------------------------------
% 
sxc11=99.201e9; %%C11  
sxc12=54.016e9; %%C12  
sxc13=50.778e9; %%C13  
sxc33=86.856e9; %%C33  
sxc44=21.1e9; %%C44  

sxe31=-7.209; %%e31  
sxe33=15.118; %%e33  
sxe15=12.322; %%e15  

sxt11=1.53e-9; %s11  
sxt33=1.5e-9; %s33  

ll1=9e-7;
ll3=9e-7;
pp3=2.5e-5;

eta1 = 0;
eta2 = 1.1e-4;
eta3 = 1.1e-4;

qqt=1.0;  
qqm=0.01;
AlphaS=0.0;
BetaS=0.0;

thick=1.0;% 
pou=7600;% 

% % 
lengthy=0.002;
lemta=0.0;

mmatn=element;
 
wjyc =[sxc11-sxc12*sxc12/sxc11, sxc13-sxc13*sxc12/sxc11,0;sxc13-sxc13*sxc12/sxc11,sxc33-sxc13*sxc13/sxc11,0; 0, 0, sxc44];%   
wjye =[0 sxe31-sxc12/sxc11*sxe31;0 sxe33-sxc13/sxc11*sxe31;sxe15 0]; %      
wjyth=[ sxt11 0;0  sxt33+sxe31*sxe31/sxc11];%   

wjyc = wjyc*(1+AlphaS*qqt+BetaS*qqm); 


matmtx=wjyc; 
matmtx1=wjye; 
matmtx2=wjyth;  

        matmtx4=[ll1-ll1*sxc13/sxc33; ll3-ll3*sxc13/sxc33; 0]; 
        matmtx5=[0;  pp3+pp3*sxe31/sxc11];

        matmtx6=[eta1-sxc12*eta2/sxc11 ; eta3-eta3*sxc13/sxc33; 0]; 
        matmtx7=[0;  eta1+eta2*sxe31/sxc11];
 
noCtrPts       = noPtsX * noPtsY;%      
noDofs         = noCtrPts * 2;%     

 

leftNodes = find(controlPts(:,1)==0)'; 
rightNodes = find(controlPts(:,1)==0.3)'; 
bottomNodes = find(controlPts(:,2)==0)';  
topNodes = find(controlPts(:,2)==0.02)';  

fixedNodes  = leftNodes;  
forcedNodes = topNodes;      

rightPoints    = controlPts(forcedNodes,:);  
rightEdgeMesh  = zeros(noElemsV,q+1);  

% 
for i=1:noElemsV
    rightEdgeMesh(i,:) = forcedNodes(i:i+q);
end

K = sparse(3*noCtrPts,3*noCtrPts);            
f  =zeros(3*noCtrPts,1);                       
fe = zeros(noDofs,1);                       
M = sparse(noDofs,noDofs);                   

jacob   = zeros(2,2);
Nxi     = zeros(1,p+1);
Neta    = zeros(1,q+1);
dNdxi   = zeros(1,p+1);
dNdeta  = zeros(1,q+1);

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(size(fixedNodes));

udofs=fixedNodes;            
vdofs=fixedNodes+noCtrPts;      
diandofs=fixedNodes+noCtrPts*2;    

fixededge = [udofs vdofs];       

if (force==0)
    uFixed = [uFixed ubar*ones(1,length(forcedNodes))];
    udofs = [udofs forcedNodes];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature


% Loop over elements (knot spans)

sdof=noCtrPts*2; %     
K2=sparse(sdof/2.0,sdof/2.0);   
K3=sparse(sdof,sdof/2.0);      
KK = sparse(noDofs,noDofs);    

         Ket=sparse(sdof,1);  
         Kvt=sparse(sdof/2.0,1); 

         Kem=sparse(sdof,1);  
         Kvm=sparse(sdof/2.0,1);  

for e=1:noElems
   idu    = index(e,1);
   idv    = index(e,2);
   xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
   etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
   connU  = elConnU(idu,:);
   connV  = elConnV(idv,:);
   
   noFnsU = length(connU);
   noFnsV = length(connV);
   
   sctr   = element(e,:);          %  element scatter vector

   % % 10.7 
   yy = mean(controlPts(sctr, 2));
  % 
   sts=exp(lemta*(yy/lengthy));
  %    
   sts1=exp(lemta*(yy/lengthy));
  %    
   sts2=exp(lemta*(yy/lengthy));
   % %
   
   sctrB  = [sctr sctr+noCtrPts]; %  vector that scatters a B matrix
   nn     = length(sctr); 
   B      = zeros(3,2*nn);
   
    edof=2*nn;  %
 

   % loop over Gauss points 
   
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      wm      = W(gp);

      % compute coords in parameter space
      Xi      = parent2ParametricSpace(xiE,pt(1));
      Eta     = parent2ParametricSpace(etaE,pt(2)); 
      J2      = jacobianPaPaMapping(xiE,etaE);
      
      dRdxi   = [];
      dRdeta  = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFnsU
       [Ni,dNi]  = NURBSbasis (connU(in),p,Xi,uKnot,weights);
       Nxi(in)    = Ni;
       dNdxi(in)  = dNi;
      end
      
      for in=1:noFnsV
       [Ni,dNi]  = NURBSbasis (connV(in),q,Eta,vKnot,weights);
       Neta(in)   = Ni;
       dNdeta(in) = dNi;
      end

      for j=1:noFnsV
          for i=1:noFnsU
              dRdxi  = [dRdxi  dNdxi(i) * Neta(j)];
              dRdeta = [dRdeta Nxi(i)   * dNdeta(j)];
          end
      end

   % derivate of R=Nxi*Neta w.r.t xi and eta
      % this is derivative of shape functions in FEM
      % the following code is for B-spline only!!!
   
%%%%%%%%%%%%%%%%%%%%%%    for     %%%%%%%%%%%%%%%%%%%%%
       ind=0;
      for j=1:noFnsV
          for i=1:noFnsU
              ind=ind+1;
              R(ind)  =Nxi(i) * Neta(j);              
          end
      end
      RR(1,1:nn)=R;
      RR(2,nn+1:nn*2)=R;

      pts = controlPts(sctr,:);
      
      % Jacobian matrix
      
      jacob(1,1) = dRdxi  * pts(:,1);
      jacob(1,2) = dRdeta * pts(:,1);
      jacob(2,1) = dRdxi  * pts(:,2);
      jacob(2,2) = dRdeta * pts(:,2);
      
      J1         = det(jacob);
      
      % Jacobian inverse and spatial derivatives
            
      dRdx       = [dRdxi' dRdeta']/jacob;%1x 2y
      
      % B matrix
      %        _                                      _
      %        |  N_1,x  N_2,x  ...      0      0  ... |
      %  B  =  |      0      0  ... N_1,y  N_2,y  ... |
      %        |  N_1,y  N_2,y  ... N_1,x  N_2,x  ... |
      %        -                                      -
      
      B(1,1:nn)       = dRdx(:,1)';%x
      B(2,nn+1:2*nn)  = dRdx(:,2)';%y
      B(3,1:nn)       = dRdx(:,2)';
      B(3,nn+1:2*nn)  = dRdx(:,1)';
      %%%         
      Be(1,1:nn) = dRdx(1:nn,1); % %%%%%%1 st row of B
      Be(2,1:nn) = dRdx(1:nn,2); %%%%%% 2 nd row of B      
      Bm(1,1:nn) = dRdx(1:nn,1); % %%%%%%1 st row of B
      Bm(2,1:nn) = dRdx(1:nn,2); %%%%%% 2 nd row of B
  
SS11=sxc11*sts;
SS12=sxc12*sts;
SS13=sxc13*sts;
SS33=sxc33*sts;
SS32=sxc13*sts;
SS44=sxc44*sts; 
%  
ee31=sxe31*sts1;
ee33=sxe33*sts1;
ee15=sxe15*sts1;
%  
eej11=sxt11*sts2;
eej33=sxt33*sts2; 

        matmtx=[SS11-SS12*SS12/SS11, SS13-SS13*SS12/SS11,0;SS13-SS13*SS12/SS11,SS33-SS13*SS13/SS11,0; 0, 0, SS44];
        matmtx1=[0 ee31-SS12/SS11*ee31;0 ee33-SS13/SS11*ee31;ee15 0]; 
        matmtx2=[ eej11 0;0  eej33+ee31.^2/SS11];      
      
      
      
    K2(sctr,sctr) = K2(sctr,sctr) - thick*Be'*matmtx2*Be*J1*J2*wt;     
    K3(sctrB,sctr) = K3(sctrB,sctr) + thick*B'*matmtx1*Be*J1*J2*wt;     
    KK(sctrB,sctrB) = KK(sctrB,sctrB) + thick*B' * matmtx* B * J1 * J2 * wt;   

    M(sctrB,sctrB) = M(sctrB,sctrB) + thick*RR' * pou * RR * J1 * J2 * wt;    

       Ket(sctrB,1) =Ket(sctrB,1)-thick*B'*matmtx*matmtx4*qqt*J1*J2*wt;  
       Kvt(sctr,1)=Kvt(sctr,1)+thick*Be'*matmtx5*qqt*J1*J2*wt; 

       Kem(sctrB,1) =Kem(sctrB,1)-thick*B'*matmtx*matmtx6*qqm*J1*J2*wm; 
       Kvm(sctr,1) =Kvm(sctr,1)+thick*Be'*matmtx7*qqm*J1*J2*wm;  
    end
end


%me+
K=[KK,K3;K3',K2];
ff1=[Ket; Kvt]+[Kem; Kvm];

f(2*noCtrPts)=0;
ffa=-ff1+f;

KK = full(KK);
K2 =-full(K2);
K3 = full(K3);
M  = full(M);

Ket= full(Ket);
Kvt= full(Kvt);
Kem= full(Kem);
Kvm= full(Kvm);
 
bcwt=mean(diag(K));
f=ffa-K(:,udofs)*uFixed';  % modify the  force vector
f=ffa-K(:,vdofs)*vFixed';
f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;

xleft=[udofs vdofs];
xleftdian = bottomNodes+noCtrPts*2 ;

K(xleft,:)=0;
K(:,xleft)=0;
K(xleft,xleft)=speye(length(xleft));
ffa(xleft)=0;

K(xleftdian,:)=0;
K(:,xleftdian)=0;
K(xleftdian,xleftdian)=speye(length(xleftdian));
ffa(xleftdian)=0;


V=0;
for i = xleftdian
    K(i,i) = 1e12*K(i,i);
    ffa(i)   = K(i,i)*(-V);
end



U=K\ffa;
toc
%%%%%%%%%%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % K2  %    
% % K3  %      
% % KK  %  

% %-----    -----

xleft = [udofs vdofs];
xleftdian = bottomNodes;

KK(xleft,:)=[];
KK(:,xleft)=[];

K2(xleftdian,:)=[];
K2(:,xleftdian)=[];

K3(xleft,:)=[];
K3(:,xleftdian)=[];

M(xleft,:)=[];
M(:,xleft)=[];

% KK = K(1:2*noCtrPts,1:2*noCtrPts);
% K2 = K(2*noCtrPts+1:3*noCtrPts, 2*noCtrPts+1:3*noCtrPts);
% K3 = K(1:2*noCtrPts, 2*noCtrPts+1:3*noCtrPts);

Keq=KK+K3*inv(K2)*K3';%      
% Keq = full(Keq);
% KK = full(KK);

% freedofs = setdiff(1:2*noCtrPts,fixededge);
% Keq = Keq(freedofs,freedofs);
% M = M(freedofs,freedofs);

KKK=Keq;
MMM=M;
 
%    
ttime=0.1;
P1=100;%  
wut=2*pi/10;
% wut=571.43*pi;
time=0:ttime:20;
F=P1*sin(wut*time);%  

dtUU=zeros(length(MMM),length(time));
dispUU=zeros(length(MMM),length(time));
vUU=zeros(length(MMM),length(time));
vvUU=zeros(length(MMM),length(time));

alpa=0.25;
detaa=0.5;
% alpa=0.7;
% detaa=0.7;

c0=1.0/(alpa*ttime^2);
c1=detaa/(alpa*ttime);
c2=1.0/(alpa*ttime);
c3=1.0/(2*alpa)-1.0;
c4=detaa/alpa-1.0;
c5=ttime/2.0*(detaa/alpa-2.0);
c6=ttime*(1-detaa);
c7=detaa*ttime;

LKK=full(KKK)+c0*full(MMM);%     
ForceA=zeros(length(MMM),length(time));
ForceA(noCtrPts*2,:)=F';%        
ForceA(xleft, :) = [];
fht = ffa(1:2*noCtrPts);
fht(xleft) = [];
fpe = ffa(2*noCtrPts+1:3*noCtrPts);
fpe(xleftdian) = [];
ffh = fht+K3/(K2)*(fpe);
for i = 1 : size(ForceA, 2)
    ForceA(:, i) = ForceA(:, i)+ffh;
end

for i=1:length(time)-1
    Forcess=ForceA(:,i)+full(MMM)*(c0*dispUU(:,i)+c2*vUU(:,i)+c3*vvUU(:,i));%t+dt     
    dispUU(:,i+1)=LKK\Forcess; %t+dt   
    vvUU(:,i+1)=c0*(dispUU(:,i+1)-dispUU(:,i))-c2*vUU(:,i)-c3*vvUU(:,i);%t+dt    
    vUU(:,i+1)=vUU(:,i)+c6*vvUU(:,i)+c7*vvUU(:,i+1);%t+dt          
end
 
Kvt_fake = Kvt;
Kvt_fake(xleftdian) = [];
Kvm_fake = Kvm;
Kvm_fake(xleftdian) = [];
fanUU=inv(K2)*(K3'*dispUU+Kvt_fake+Kvm_fake);

outputNodes = 1122;

Awe=dispUU(outputNodes*2-1,:)';%A  x
AwA=dispUU(outputNodes*2,:)';%A  z
Adian=fanUU(outputNodes,:)';%A  
Aresult=[Awe Adian];

% figure;
% plot(time(2:end),Awe(2:end));
% % % figure;
% % % plot(time(2:end),AwA(2:end));
% figure;
% plot(time(2:end),Adian(2:end));

figure;
plot(time(1:end),Awe(1:end));
figure;
plot(time(1:end),AwA(1:end));
figure;
plot(time(2:end),Adian(2:end));

% figure;
% plot(time,Aci);
% % %%%%%%%%%%%%%%%%%%%%

outputNodes1 = 1089:1122;
totalTime = 20; %    
outputTime = 5; %       
OutputTime1 = round(outputTime / totalTime * (length(time)));
dispTimeWeiyi = dispUU(:, OutputTime1); %       
dispTimePhi = fanUU(:, OutputTime1); %       

UX = dispTimeWeiyi(outputNodes1);
UY = dispTimeWeiyi(1122 + outputNodes1);
PHI = dispTimePhi(outputNodes1);
