function [xi,epsilon_dx,epsilon_gradx,iter]=CGNE_Val(CCD_img_noise,FCS,Gr_x,Gr_y, param_CGNE)

[R,C,W,S]=size(FCS);

sqrt_mu_x=sqrt(param_CGNE.mu_x);
sqrt_mu_y=sqrt(param_CGNE.mu_y);
sqrt_mu_lambda=sqrt(param_CGNE.mu_lambda); 


%*** Init

if(param_CGNE.x0==0)
    xi=zeros(R,C,W);
elseif (size(param_CGNE.x0,3)>1)
    xi=param_CGNE.x0;
else
    xi=param_CGNE.x0(:,:,ones(W,1));  
end

%%
%*** Considering the noise

if(size(param_CGNE.Sig_noise)==1)
    FCS=FCS/param_CGNE.Sig_noise;
    CCD_img_noise=CCD_img_noise/param_CGNE.Sig_noise;
    sqrt_mu_x=sqrt(mu_x)/(param_CGNE.Sig_noise);
    sqrt_mu_y=sqrt(mu_y)/(param_CGNE.Sig_noise);
    sqrt_mu_lambda=sqrt(mu_lambda)/(param_CGNE.Sig_noise); 
    
elseif(size(param_CGNE.Sig_noise)>=2)
    for s=1:S
        
    FCS(:,:,:,s)=squeeze(FCS(:,:,:,s))./repmat(squeeze(param_CGNE.Sig_noise(:,:,s)), [1 1 W]);
    CCD_img_noise(:,:,s)=CCD_img_noise(:,:,s)./squeeze(param_CGNE.Sig_noise(:,:,s));
    
    end
end
%%%

%%
% *** Computing the corresponding Epsilon( dx(n)=norm(oi(:))/norm(xi(:)) )
% ***                  and        Epsilon( Grdi/Res)
%      of the param_CGNE.epsilonx given (input argument param_CGNE.epsilonx)

% 
EpGrdX=exp(-9.5)*param_CGNE.epsilonx.^(.52); % Epsilon( Grdi/Res)
EpdX=exp(- 6.62)*param_CGNE.epsilonx.^(.51); % Epsilon( dx(n)=norm(oi(:))/norm(xi(:))

%%%

ri_A=zeros(R,C,S);
pi_A=zeros(R,C,W);
for s=1:S
    ri_A(:,:,s)=CCD_img_noise(:,:,s)-sum(FCS(:,:,:,s).*xi,3);
end
ri_x=(-sqrt_mu_x*cat(2,Gr_x(:,1:C-1,ones(W,1)).*diff(xi,1,2),zeros(R,1,W)));
ri_y=(-sqrt_mu_y*cat(1,Gr_y(1:R-1,:,ones(W,1)).*diff(xi,1,1),zeros(1,C,W))); 
ri_z=(-sqrt_mu_lambda*cat(3,diff(xi,1,3),zeros(R,C,1)));% diff(xi,1,3), diff(xi,1,2) and diff(xi,1,1) <==> difference along lambda, along x and along y

for s=1:S
    pi_A=pi_A+FCS(:,:,:,s).*( squeeze(ri_A(:,:,s,ones(W,1))) );
end
pi_x=-sqrt_mu_x*diff(cat(2,zeros(R,1,W),(Gr_x(:,:,ones(W,1)).*ri_x)),1,2);
pi_y=-sqrt_mu_y*diff(cat(1,zeros(1,C,W),(Gr_y(:,:,ones(W,1)).*ri_y)),1,1);
pi_z=-sqrt_mu_lambda*diff(cat(3,zeros(R,C,1),ri_z),1,3);
pii=pi_A+pi_x+pi_y+pi_z;
clear pi_A pi_x pi_y pi_z;
zi_A=pii;
%**************************************************************************

wi_A=zeros(R,C,S);
count=0;
epsilon_dx=100*EpdX;
epsilon_gradx=100*EpGrdX;
while ((epsilon_dx >=EpdX || epsilon_gradx >=EpGrdX) ) && count < param_CGNE.stop_crit

    

    for s=1:S
        wi_A(:,:,s)=sum(FCS(:,:,:,s).*pii,3);
    end
    wi_x=sqrt_mu_x*cat(2,Gr_x(:,1:C-1,ones(W,1)).*diff(pii,1,2),zeros(R,1,W));
    wi_y=sqrt_mu_y*cat(1,Gr_y(1:R-1,:,ones(W,1)).*diff(pii,1,1),zeros(1,C,W)); 
    wi_z=sqrt_mu_lambda*cat(3,diff(pii,1,3),zeros(R,C,1));% diff(xi,1,3), diff(xi,1,2) and diff(xi,1,1) <==> difference along lambda, along x and along y
    norm2_zi=sum(sum(sum(zi_A.^2)));
    alphai=norm2_zi/(sum(sum(sum(wi_A.^2))) +sum(sum(sum(wi_x.^2)))+sum(sum(sum(wi_y.^2)))+sum(sum(sum(wi_z.^2))));
    
    oi=xi;
    xi= xi + alphai*pii;
    oi=oi-xi; % This variable will be used to compute the succesive differences : x_i-x_i+1
    
    ri_A=ri_A - alphai*wi_A;
    ri_x=ri_x - alphai*wi_x;
    ri_y=ri_y - alphai*wi_y;
    ri_z=ri_z - alphai*wi_z;
    
    zi_A=zi_A*0;
    for s=1:S
        zi_A=zi_A+FCS(:,:,:,s).*( squeeze(ri_A(:,:,s,ones(W,1))) );        
    end
    zi_x=-sqrt_mu_x*diff(cat(2,zeros(R,1,W),(Gr_x(:,:,ones(W,1)).*ri_x)),1,2);
    zi_y=-sqrt_mu_y*diff(cat(1,zeros(1,C,W),(Gr_y(:,:,ones(W,1)).*ri_y)),1,1);
    zi_z=-sqrt_mu_lambda*diff(cat(3,zeros(R,C,1),ri_z),1,3);
    zi_A=zi_A+zi_x+zi_y+zi_z;
    
    betai= sum(sum(sum(zi_A.^2)))/norm2_zi;    
    pii=zi_A + betai*pii;
    
%% Epsilon dX and GrdX
    Res=sum(sum(sum(ri_A.^2)))+sum(sum(sum(ri_x.^2)))+sum(sum(sum(ri_y.^2)))+sum(sum(sum(ri_z.^2)));
    Grdi=sqrt(norm2_zi);
 
    epsilon_dx=norm(oi(:))/norm(xi(:));
    epsilon_gradx=Grdi/Res;
%% 
count=count+1
epsilon_dx;
%imagesc(squeeze(xi(:,:,30)))
end

iter=count;
%Xi=xi;


end
