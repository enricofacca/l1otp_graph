close all

folder = 'data_test_cases/'
folder = '/home/fh/runs_matlab_graph/'
folder = 'make_graph/'

%file=strcat('/home/fh/runs_matlab_graph/erdos_graph_p',level,'.dat')
file=strcat(folder,label_graph,'.dat')
if (exist(file, 'file') == 2 )
  [coord,topol]=read_graph(file);
  ntdens=size(topol,2);
  npot=size(coord,2);
else
  disp('graph not found')
end


sizecell=ones(ntdens,1);
file_sizecell=strcat(folder,label_problem,'_sizecell.dat')
if (exist(file_sizecell, 'file') == 2 )
  sizecell=read_td(file_sizecell);
else
  for i=1:ntdens
    sizecell(i)=norm(coord(:,topol(1,i))-coord(:,topol(2,i)));
  end
end

max(topol(:));
npot;


			  % set file related to spatial discretization
file_rhs=strcat(folder,label_problem,'_rhs.dat')
%file_rhs=strcat('/home/fh/runs_matlab_graph/erdos_rhs_p',level','nnz010.dat')

rhs=read_td(file_rhs);
[m,isource]=max(rhs(:));
npot=size(rhs,1);

optpot=zeros(size(rhs,1),1);


file_optpot=strcat(folder,label_problem,'_optpot.dat');
if (exist(file_optpot, 'file') == 2 )
    disp(strcat('Optpot form file',file_optpot) ) 
    optpot=read_td(file_optpot);
else
  optpot=zeros(npot,1);
  disp(strcat('Optpot not found') ) 
end



file=strcat(folder,label_problem,'_opttdens.dat');
if (exist(file, 'file') == 2 )
    disp(strcat('OptTdens form file',file) ) 
    opttdens=read_td(file);
else
  disp(strcat('OptTdens not found') ) 
  opttdens=zeros(ntdens,1);
end



tdens=ones(ntdens,1);
pot=zeros(npot,1);

return

nconnection=zeros(npot,1);
for i=1:ntdens
    n1=topol(1,i);
    n2=topol(2,i);
    nconnection(n1)=nconnection(n1)+1;
    nconnection(n2)=nconnection(n2)+1;
end

%sink=rhs(rhs(:)<0.0);

max(nconnection(:))
%lower_bound_tdens=min(abs(sink(:)))/(1.01*max(nconnection(:)));

%file=strcat('~/bbmuffe/runs/064_select/tdens00000005.dat');
%tdens=read_td(file);

%file=strcat('~/bbmuffe/runs/noselect_064_10/pot00000026.dat');
%pot=read_td(file);

bar_coord=zeros(2,ntdens);
for i= 1:ntdens
    n1=topol(1,i);
    n2=topol(2,i);
    bar_coord(:,i)=0.5*[coord(:,n1)+coord(:,n2)];    
end

optdens=zeros(ntdens,1);


factor=max(rhs);


for i= 1:ntdens
    n1=topol(1,i);
    n2=topol(2,i);
    x1=coord(1,n1);
    x2=coord(1,n2);
    y1=coord(2,n1);
    y2=coord(2,n2);
    xbar=bar_coord(1,i);
    ybar=bar_coord(2,i);
    
    if ( x1 == x2 ) 
        optdens(i)=0.0;
    elseif( x1~=x2 && y1~=y2 ) 
        optdens(i)=0.0;
    else
        if (ybar<0.25 || ybar>0.75 )
             optdens(i)=0.0;
        else
            if( xbar >= 0.125 && xbar <=0.375 )
                optdens(i) = factor*fix((xbar-0.125)/abs(x1-x2)+1);
            elseif ( xbar > 0.375 && xbar <0.625 )
                optdens(i) = factor*fix((0.375-0.125)/abs(x1-x2)+1);
            elseif ( xbar >= 0.625 && xbar <=0.875 )
                optdens(i) = factor*fix((0.875-xbar)/abs(x1-x2)+1);
            else
                optdens(i)=0;                
            end
        end
    end
    
    
    
    
end



%write2td('eik_064_opttdens.dat',optdens)


[m,isource]=max(rhs(:))


max(nconnection(:))
%lower_bound_tdens=min(abs(sink(:)))/(max(nconnection(:)))
coord(:,isource)

return
addpath('./dijkstra/')

nodes=[(1:npot);coord]';
segments=[(1:ntdens);topol]';
tic
[dist,path] = dijkstra(nodes,segments,isource);
toc
optpot=dist';

write2td('optpot.dat',dist')



bar_coord=zeros(2,ntdens);
for i= 1:ntdens
    n1=topol(1,i);
    n2=topol(2,i);
    bar_coord(:,i)=0.5*[coord(:,n1)+coord(:,n2)];    
end

optdens=zeros(ntdens,1);

factor=max(rhs);


for i= 1:ntdens
    n1=topol(1,i);
    n2=topol(2,i);
    x1=coord(1,n1);
    x2=coord(1,n2);
    y1=coord(2,n1);
    y2=coord(2,n2);
    xbar=bar_coord(1,i);
    ybar=bar_coord(2,i);
    
    if ( x1 == x2 ) 
        optdens(i)=0.0;
    elseif( x1~=x2 && y1~=y2 ) 
        optdens(i)=0.0;
    else
        if (ybar<0.25 || ybar>0.75 )
             optdens(i)=0.0;
        else
            if( xbar >= 0.125 && xbar <=0.375 )
                optdens(i) = factor*fix((xbar-0.125)/abs(x1-x2)+1);
            elseif ( xbar > 0.375 && xbar <0.625 )
                optdens(i) = factor*fix((0.375-0.125)/abs(x1-x2)+1);
            elseif ( xbar >= 0.625 && xbar <=0.875 )
                optdens(i) = factor*fix((0.875-xbar)/abs(x1-x2)+1);
            else
                optdens(i)=0;                
            end
        end
    end
    
    
    
    
end
