clear all
close all
clc

n_from=2:50;
for ind=1:length(n_from)
%%% given data
m=2; %number of points in each space
n=n_from(ind);
d=ones(m,m)-eye(m); %distance matrices
h=ones(n,n)-eye(n);
mux=1/m*ones(m,1); %measures
muy=1/n*ones(n,1);

% %%% construct the tensor
A=[];
for ii=1:m
    for jj=1:n
        A=[A;compute_row(ii,jj,d,h)];
    end
end

evalues=eig(A);
ff=find(evalues<0);

neg_evalues(ind)=length(ff);
end


plot(n_from,neg_evalues,'k o-','LineWidth',1.5)
xlabel('size of second space')
ylabel('number of negative eigenvalues')
title('Delta_n spaces from panel (A)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions
function row=compute_row(ii,jj,d,h)
m=size(d,1);
n=size(h,1);
row=[];
for k=1:m
    for l=1:n
        row=[row abs(d(ii,k)-h(jj,l))];
    end
end
end
