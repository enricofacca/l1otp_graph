function [myincidence]=myincidence(topol)
  ncell=size(topol,2);
  nnode=max(topol(:));
  ia=zeros(2*ncell,1);
  ja=zeros(2*ncell,1);  
  a=zeros(2*ncell,1);
  ia(1:2:2*ncell)=1:ncell;
  ia(2:2:2*ncell)=1:ncell;

  m=0;
  for icell=1:ncell;
    m=m+1;
    ja(m)=topol(1,icell);
    a(m)=1;
    m=m+1;
    ja(m)=topol(2,icell);
    a(m)=-1;
  end

  myincidence=sparse(ia,ja,a);
  myincidence=myincidence';


  
end


  
  
