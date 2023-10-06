%**************************************************************************
%*     Exercise 17.3                                                      *
%*     filename: ch17pr03.m                                               *
%*     program listing number: 17.3                                       *
%*                                                                        *
%*     This program simulates two-dimensional percolation.                *
%*     Hoshen and Kopelman algorithm is used for cluster labeling.        *
%*                                                                        *
%*     Programed by Ryoichi Kawai for Computational Physics Course.       *
%*     Last modification:  03/17/2017.                                    *
%**************************************************************************
close all

N=32;
% p=0.59 is close to the transition point
p=0.4;

% Graphics setting
plot([0,N+1,N+1,0,0],[0,0,N+1,N+1,0]);
axis equal
axis([0 N+1 0 N+1])
hold on

% place atoms at random
r=rand(N,N);
lattice = r < p;

% draw clusters
for i=1:N
    for j=1:N
        if lattice(i,j)
            rectangle('Position',[i,j,1,1],'Curvature',[1 1],'FaceColor','b');
            drawnow
        end
    end
end
hold on


% Labeling: 1st pass (Making initial labels and map)

label=zeros(N,N);  % allocate array
remap=zeros(N*N);

new=0;
for j=1:N
    for i=1:N
        
        i1 = i-1;
        j1 = j-1;
        
        if lattice(i,j)>0
            if i1 > 0
                left=label(i1,j); % neighbor (left).
            else
                left=0;  % outside the box.
            end
            if j1 > 0
                down=label(i,j1); % neighbor (down).
            else
                down=0;  % outside the box.
            end
            
            if down==0 && left==0  % if both are unocupied
                new=new+1;         % create a new cluster.
                label(i,j) = new;
                remap(new)=new;
                
            elseif down*left>0  % both are occupied.
                
                if down==left        % if they are the same cluster.
                    label(i,j)=left; % join to the cluster.
                    
                else   % connecting two different clusters.
                    found = false;
                    while not(found)
                        if remap(left)==left
                            found = true;
                        else
                            left=remap(left);
                        end
                    end
                    found = false;
                    while not(found)
                        if remap(down)==down
                            found = true;
                        else
                            down=remap(down);
                        end
                    end
                    
                    if left==down  % they are again the same
                        label(i,j)=left;
                    else
                        nmax=max(left,down);
                        nmin=min(left,down);
                        label(i,j)=nmin;   % coalesce two clusters
                        remap(nmax)=nmin;  % add to the chain
                    end
                end
                
            elseif down>0 % only down neighbor is occupied
                label(i,j) = down;  % join to the neighbor
                
            else  % only the left neighbor is occupied
                label(i,j) = left;  % join to the neighbor
            end
        end
        
    end
end
% Labeling: 2nd pass (Collapse the label in the same cluster)

nmax = max(label(:));
for i=nmax:-1:1
    label(label==i)=remap(i);
end

% Labeling: 3rd pass (Make the label continuous)
%           This procedure is not essential.

j=0;
for i=1:nmax
    if remap(i)==i
        j=j+1;
        label(label==i) = j;
    end
end

% Identify the percolation
nmax = max(label(:));
size = zeros(nmax,1);
maxsize=0;
largest=1;
fprintf('Cluster   Size  Percolation\n')
for i=1:nmax
    percx = any(label(1,:)==i)*any(label(N,:)==i) >0;
    percy = any(label(:,1)==i)*any(label(:,N)==i) >0;
    size(i)=  sum(label(:)==i);
    if size(i) > maxsize
        largest = i;
        maxsize = size(i);
    end
    if percx || percy
        perc='YES';
    else
        perc=' NO';
    end
    fprintf('%5d:  %6d,    %s\n', i, size(i), perc)
end

fprintf('%d    %d\n', largest, size(largest))

for i=1:N
    for j=1:N
        if label(i,j) == largest
            rectangle('Position',[i,j,1,1],'Curvature',[1 1],'FaceColor','r');
            drawnow
        end
    end
end
figure

h=histogram(size,20);