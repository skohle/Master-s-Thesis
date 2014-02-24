function [ P,p ] = DEIM( POD, x1 , x2, mode )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Sophia Kohle
%
% Computing the interpolation indices for the POD-DEIM method
%
% INPUT: POD basis
%
% OUTPUT: DEIM interpoalation indices p
%         Projection matrix P
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,r]=size(POD);
[rho,pindex(1,:)]=max(sign(POD(:,1)).*POD(:,1));
U=POD(:,1);
e_p=zeros(m,1);
e_p(pindex)=1;
P=e_p;
p=pindex;
if mode==0;
    figure(110);
    plot(x1(pindex(1,:)),x2(pindex(1,:)),'*')
    text(x1(pindex(1,:)),x2(pindex(1,:)), [' ' num2str(1) '.'])
end
for j=2:r
    c=(P'*U)\P'*POD(:,j);
    r=POD(:,j)-U*c;
    [rho1,pindex(j,:)]=max(sign(r).*r);
    rho=[rho;sign(rho1)*rho1];
    U=[U,POD(:,j)];
    e_p=zeros(m,1);
    e_p(pindex(j,:))=1;
    P=[P e_p];
    p=[p; pindex(j,:)];
    
    if mode==0
        hold on
        plot(x1(pindex(j,:)),x2(pindex(j,:)),'*')
        title('DEIM interpolation points')
        text(x1(pindex(j,:)),x2(pindex(j,:)), [' ' num2str(j) '.'])
        hold off
    end
end  

P=sparse(P);

