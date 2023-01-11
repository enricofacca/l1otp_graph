function [ coord,topol ] = read_graph( file)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(file);

tline = fgetl(fid);
info=sscanf(tline,'%d  *');
nnode=info(1);

tline = fgetl(fid);
info=sscanf(tline,'%d %d *');
ncell=info(1);
nnodeincell=info(2);

coord=zeros(2,nnode);
topol=zeros(nnodeincell,ncell,'uint64');

for i=1:nnode
    tline = fgetl(fid);
    info=sscanf(tline,'%f %f *');
    coord(1:2,i)=info(1:2);
end

for i=1:ncell
    tline = fgetl(fid);
    info=sscanf(tline,'%d %d *');
    topol(1:2,i)=info(1:2);
end



end

