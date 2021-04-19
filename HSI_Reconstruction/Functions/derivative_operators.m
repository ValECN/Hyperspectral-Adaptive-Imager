function [D_x, D_y, D_lambda]= derivative_operators(Gr_x,Gr_y,derivative_order_xy,derivative_order_lambda,Hyper_cube_size)


R=Hyper_cube_size(1);
C=Hyper_cube_size(2);
W=Hyper_cube_size(3);


if(derivative_order_xy==2)
    
size_Gr_x=size(Gr_x);
size_Gr_y=size(Gr_y);

if(size_Gr_x(2)~=(C-2) || size_Gr_y(1)~=(R-2) )
    
    error('Dimension mismatches: Gr_x & Gr_y should be of size [ R X (C-2) ] & [ (R-2) X C ] respectively');    
end    


%*** Construction of the  SECOND order derivatives finite difference operators 
%***               along x & y & lambda [1 -2  1] 
%***
%***    & Edge Detection using Gr_x & Gr_y Masks****

%**************
% D_xx  & Edge Detection using Gr_x Masks****

J_in_new_1=[1:R*(C-2)*W]+ 2*W*ceil([1:R*(C-2)*W]/(W*(C-2)))-2*W;
J_in_new_2=J_in_new_1+1*W;
J_in_new_3=J_in_new_1+2*W;

I_in_new_1=[[1:length(J_in_new_1)]';[1:length(J_in_new_1)]';[1:length(J_in_new_1)]'];


vect_Gr_x=Gr_x';
vect_Gr_x=vect_Gr_x(:);
dupliacte_vect_Gr_x=vect_Gr_x(:,ones(W,1))';
vect_duplicate_vect_Gr_x=dupliacte_vect_Gr_x(:); 

val_vect=[vect_duplicate_vect_Gr_x ; -2*vect_duplicate_vect_Gr_x; vect_duplicate_vect_Gr_x];

D_x=sparse(I_in_new_1,[J_in_new_1';J_in_new_2';J_in_new_3'],val_vect,R*(C-2)*W,R*C*W); % D_xx


%**************
% D_yy & Edge Detection using Gr_y Masks****

J_in_new_1=[1:(R-2)*C*W];
J_in_new_2=J_in_new_1+1*C*W;
J_in_new_3=J_in_new_1+2*C*W;

I_in_new_1=[[1:length(J_in_new_1)]';[1:length(J_in_new_1)]';[1:length(J_in_new_1)]'];


vect_Gr_y=Gr_y';
vect_Gr_y=vect_Gr_y(:);
dupliacte_vect_Gr_y=vect_Gr_y(:,ones(W,1))';
vect_duplicate_vect_Gr_y=dupliacte_vect_Gr_y(:); 

val_vect=[vect_duplicate_vect_Gr_y;-2*vect_duplicate_vect_Gr_y;vect_duplicate_vect_Gr_y];
D_y=sparse(I_in_new_1,[J_in_new_1';J_in_new_2';J_in_new_3'],val_vect,(R-2)*C*W,R*C*W);



elseif(derivative_order_xy==1)
    
    %%
%*** Construction of the FIRST order finite difference operators along x & y ****        
%***                          [-1 1]
%***     & Edge Detection using Gr_x & Gr_y Masks****


size_Gr_x=size(Gr_x);
size_Gr_y=size(Gr_y);

if(size_Gr_x(2)~=C || size_Gr_y(1)~=R )   
    error('Dimension mismatches: Gr_x & Gr_y should be of size [ R X (C-1) ] & [ (R-1) X C ] respectively');
end


%**************
% D_x & Edge Detection using Gr_x Masks****

I_in=[1:C*W*R]; 
I_in_1=[1:(R*C*W-W)];
J_in_1=[(W+1):C*R*W];
I_in_2=[I_in, I_in_1];
J_in_2=[I_in, J_in_1];

vect_Gr_x=Gr_x';
vect_Gr_x=vect_Gr_x(:);
dupliacte_vect_Gr_x=vect_Gr_x(:,ones(W,1))';
vect_duplicate_vect_Gr_x=dupliacte_vect_Gr_x(:); 

val_vect=[-vect_duplicate_vect_Gr_x ; vect_duplicate_vect_Gr_x(1:R*C*W-W)];
D_x=sparse(I_in_2,J_in_2,val_vect,W*R*C,W*R*C);

%**************
% D_y
I_in=[1:C*W*R]; 
% J_in=ceil(I_in/W);
J_in_1=[(W*C+1):C*R*W];
I_in_1=[1:(R*C*W-W*C)];
I_in_2=[I_in, I_in_1];
J_in_2=[I_in, J_in_1];

vect_Gr_y=Gr_y';
vect_Gr_y=vect_Gr_y(:);
dupliacte_vect_Gr_y=vect_Gr_y(:,ones(W,1))';
vect_duplicate_vect_Gr_y=dupliacte_vect_Gr_y(:); 


val_vect=[-vect_duplicate_vect_Gr_y ; vect_duplicate_vect_Gr_y(1:R*C*W-C*W)];
D_y  =sparse(I_in_2,J_in_2,val_vect,W*R*C,W*R*C);


%  Edge Detection using Gr_y Masks****


end

%% ************************************************************************


if(derivative_order_lambda==1)
    
    % %*** Construction of the FIRST order finite difference operator along
    %                         lambda  [1 -1]
    
    %**************
    % D_lbda Wihtout considering Edges
    temp_w=[-ones(1,W-1) 0; [0,ones(1,W-1)]]'; % [-1 1] !!!!
    temp_w=repmat(temp_w,R*C,1);
    
    D_lambda=spdiags(temp_w,[0 1],length(temp_w),length(temp_w)); % D_lbda;
    
    %*************************************************************************
    
elseif(derivative_order_lambda==2)
    
    
    % %*** Construction of the SECOND order finite difference operator along
    %                         lambda  [1 -2  1]
    
    % D_lbdalbda Wihtout considering Edges
    J_in_new_1=[1:R*C*(W-2)]+ 2*ceil([1:R*C*(W-2)]/(W-2))-2;
    J_in_new_2=J_in_new_1+1;
    J_in_new_3=J_in_new_1+2;
    
    I_in_new_1=[[1:length(J_in_new_1)]';[1:length(J_in_new_1)]';[1:length(J_in_new_1)]'];
    val_vect=ones(length(J_in_new_1),1);
    val_vect=[val_vect;-2*val_vect;val_vect]; % [1 -2  1]
    
    D_lambda=sparse(I_in_new_1,[J_in_new_1';J_in_new_2';J_in_new_3'],val_vect,R*C*(W-2),R*C*W); % D_lbdalbda
    
   
end

end