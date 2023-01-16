function [] = writegrid( file,coord,topol)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


fid = fopen(file,'w');

fprintf(fid, '%d \n', size(coord,2));
fprintf(fid, '%d 2 \n', size(topol,2));


for i=1:size(coord,2)
  fprintf(fid, '%18.12e %18.12e \n', coord(1,i),coord(2,i));
end

for i=1:size(topol,2)
  fprintf(fid,'%d %d %d \n',topol(1,i),topol(2,i),i);
end



end

