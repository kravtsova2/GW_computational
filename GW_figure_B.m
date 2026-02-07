clear all
close all
clc

%%% load the data
load Lotka1.dat
load Lotka2.dat

%%% for the plot
npts=500;

ind1=randsample(size(Lotka1,1),npts);
ind2=randsample(size(Lotka2,1),npts);

Traj1=Lotka1(sort(ind1),:);
Traj2=Lotka2(sort(ind2),:);

figure(1)
plot3(Traj1(:,1),Traj1(:,2),Traj1(:,3),'bo--','Linewidth',0.5,'MarkerFaceColor','b','MarkerSize',5)
figure(2)
plot3(Traj2(:,1),Traj2(:,2),Traj2(:,3),'mo--','Linewidth',0.5,'MarkerFaceColor','m','MarkerSize',5)


m=50;
n_from=10:25;
for ind=1:length(n_from)

n=n_from(ind);

ind1=randsample(size(Lotka1,1),m);
ind2=randsample(size(Lotka2,1),n);

Traj1=Lotka1(sort(ind1),:);
Traj2=Lotka2(sort(ind2),:);

% figure(1)
% plot3(Traj1(:,1),Traj1(:,2),Traj1(:,3),'m *--','Linewidth',0.5)
% figure(2)
% plot3(Traj2(:,1),Traj2(:,2),Traj2(:,3),'b *--','Linewidth',0.5)

%%% compute matrices of intrinsic distances
d=my_geo(Traj1);
h=my_geo(Traj2);

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

figure(3)
plot(n_from,neg_evalues,'k o-','LineWidth',1.5)
xlabel('size of second space')
ylabel('number of negative eigenvalues')
title('Curves from panel (B) (up to 25 points in second curve)')


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


function D=my_geo(X)
%%% This function takes a trajectory X where each column of X gives values for a
%%% coordinate of time series. For example, time series X on the plane
%%% (with 2 coordinates) will have two columns. Number of rows is the
%%% length of time series

%%% The function outputs a matrix of intrinsic distances between ALL points
%%% along the time series. This matrix can be used for Gromov-Wasserstein
%%% distance computations

D=zeros(size(X,1),size(X,1));

for i=1:size(X,1)-1
D(i,i+1)=norm(X(i,:)-X(i+1,:));
end

for i=1:size(X,1)
    for j=i+2:size(X,1)
D(i,j)=D(i,j-1)+D(j-1,j);
    end
end

D=D+D';
end



