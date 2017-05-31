function success = det_save_cemv(output_file, edgelist)
% function success = det_save_cem_v1(edgelist)
% This function save an list of contour fragments (edgelist) to a file in
% the .cem format (version 1, simple contours, no appearance info) 
%  Returns:  edgelist - a cell array of edge lists in row, column, orientation coords in
%                       the form
%                      { [r1 c1 theta1   
%                         r2 c2 theta2   
%                         ...
%                         rN cN thetaN],
%                        [r1 c1 theta1   
%                         r2 c2 theta2   
%                         ...
%                         rN cN thetaN],
%                         ...
%                      }

% This is the format returned by Kovesi's edge linker
% (c) Nhon Trinh
% Date: Feb 6, 2009

%% Open file

fid = fopen(output_file, 'w');

if (fid == -1)
  fprintf(2, 'Could not open file %s for writing.\n', output_file);
  success = 0;
  return;
end;

%% Output header information
fprintf(fid, '# CONTOUR_EDGE_MAP : Logical-Linear + Shock_Grouping\n');
fprintf(fid, '# .cem files\n');
fprintf(fid, '#\n');
fprintf(fid, '# Format :\n');
fprintf(fid, '# Each contour block will consist of the following\n');
fprintf(fid, '# [BEGIN CONTOUR]\n');
fprintf(fid, '# EDGE_COUNT=num_of_edges\n');
fprintf(fid, '# [Pixel_Pos]  Pixel_Dir Pixel_Conf  [Sub_Pixel_Pos] Sub_Pixel_Dir Sub_Pixel_Conf\n');
fprintf(fid, '# ...\n');
fprintf(fid, '# ...\n');
fprintf(fid, '# [END CONTOUR]\n');
fprintf(fid, '\n');


%% Body

% contour count and edge count
num_contours = length(edgelist);
num_edges = 0;
for i = 1 : length(edgelist)
  num_edges = num_edges + size(edgelist{i}, 1);
end;

fprintf(fid, 'CONTOUR_COUNT=\n');
fprintf(fid, 'TOTAL_EDGE_COUNT=\n'); 


% position, orientation, and confidence level for each point
for ic = 1 : length(edgelist)
  sub_contour = edgelist{ic};
  
  edge_count = size(sub_contour, 1);
  
  % matlab starts indexing at 1, while coordinate axis start at 0
  % need to subtract 1 from [row, column] index to get the coordinates
%   sub_x = sub_contour(:, 2)-1; 
%   sub_y = sub_contour(:, 1)-1;
  sub_x = sub_contour(:, 1); 
  sub_y = sub_contour(:, 2);
  dir = sub_contour(:,3);
  str = sub_contour(:,4);
  %sub_dir = sub_contour(:, 3);
  %sub_conf = zeros(edge_count, 1);
  
  %pix_x = round(sub_x);
  %pix_y = round(sub_y);
  %pix_dir = sub_dir;
  %pix_conf = zeros(edge_count, 1);

  fprintf(fid, '[BEGIN CONTOUR]\n');
  fprintf(fid, 'EDGE_COUNT=%d\n', edge_count);
  for ip = 1 : edge_count
    
    % ' [Pixel_Pos]  Pixel_Dir Pixel_Conf  [Sub_Pixel_Pos] Sub_Pixel_Dir Sub_Pixel_Conf'
    % [%d, %d]\t%lf\t%lf\t[%lf, %lf]\t%lf\t%lf
    fprintf(fid, ' [%d, %d] %f %f [%f, %f] %f %f\n', ...
      0, 0, 0.0, 0.0, ...
      sub_x(ip), sub_y(ip), dir(ip), str(ip));
  end
  fprintf(fid, '[END CONTOUR]\n\n');  
end;
  
fclose(fid);
