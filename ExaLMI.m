clc
clear all;
A=[0 1 0;9.8 0 1;0 -10 -10];
B=[0;1;10];
L=[1 2 3];
H=[1 6 7];
%%
setlmis([])             %reset all parameters  start with setlmis and end with getlmis
Y=lmivar(1,[3,1]);      %definning symmetric 3*3 matrix
gamma=lmivar(1,[1,1]);
%the first LMI Y>0
lmiterm([-1 1 1 Y],1,1) %Y>0
%lmiterm([1 1 1 0],0)
% the second LMI
lmiterm([2 1 1 Y],A,1,'s')             %Defining A for(i=1 j=1)---AX + X'A' 
lmiterm([2 1 1 0],B*L+L'*B');
lmiterm([2 1 2 0],eye(3));
lmiterm([2 1 3 Y],1,H');%YH'
lmiterm([2 2 1 0],eye(3));
lmiterm([2 2 2 0],-1*eye(3));
lmiterm([2 2 3 0],zeros(3,1));
lmiterm([2 3 1 Y],H,1);
lmiterm([2 3 2 0],zeros(1,3));
lmiterm([2 3 3 gamma],-1,eye(1));
%%
LMISYS = getlmis; 
[tmin,xfeas] = feasp(LMISYS); 
tmin
Y = dec2mat(LMISYS,xfeas,Y);  %getting the value of lmi variable
gamma=dec2mat(LMISYS,xfeas,gamma)
%%
%%
% % LMISYS = getlmis;                       %start with setlmis and end with getlmis
% % [lopt,xopt] = gevp(LMISYS,1);           %using generalized eign value problem for solver
% % Y = dec2mat(LMISYS,xopt,Y)            %getting the value of lmi variable
% % gamma = dec2mat(LMISYS,xopt,gamma)