function [level_handles,level_labels_out] = getlegendinfo;

% GETLEGENDINFO returns line plot info for creating legends.
%
%     [H,LABELS] = GETLEGENDINFO returns a vector of handles to the line
%     plots created by PLOTCLAW1, PLOTCLAW2 and PLOTCLAW3.  The handles are
%     in the vector H and optional labels are in the cell matrix
%     LABELS. This information can be passed directly to the Matlab LEGEND
%     command to add legends to your plots.   Additionally, this information
%     can be combined with other line plots you may add to the plot (1d
%     reference solutions, for example) to further describe your plot.
%
%     When you pass this data to the LEGEND command, the symbols you
%     set in ScatterStyle, PlotStyle, or LineStyle will be correctly
%     displayed on the legend.
%
%     The strings return in LABELS are default strings 'Level 1', 'Level 2'
%     and so on.
%
%     H = GETLEGENDINFO returns the vector of handles only.  The user can
%     then supply their own labels for the line plots created by Clawpack graphics.
%
%     Example :
%
%         % AFTERFRAME file for a scatter plot of AMR data :
%         [h_amr, labels_amr] = getlegendinfo;
%         % Use default labels in 'Level 1','Level 2','Level 3', etc.
%         legend(h_amr,labels_amr);
%
%     Example :
%
%         % AFTERFRAME file for a scatter plot of AMR data :
%         h_amr = getlegendinfo;
%         % Create your own labels
%         legend(h_amr,{'My Level 1', 'My Level 2', 'My Level 3'});
%
%     Example :
%
%         % AFTERFRAME file for a scatter plot of AMR data :
%         % Add 1d reference solution  line plot to legend
%         [data1d,tref] = readamrdata(1,Frame,'./qref/');
%         [qref,xref,p] = plotframe1ez(data1d,mq);
%         [h_amr,labels_amr] = getlegendinfo;
%         legend([h_amr,p],{labels_amr{:},'1d reference solution'});
%
%     See Also LEGEND, CELL, PAREN.

amrlines = get_lines;

level_handles = [];
level_labels = {};
for level = 1:length(amrlines),
  svec = amrlines{level};
  if (~isempty(svec))
    level_labels{level} = sprintf('Level %d',level);
    level_handles(level) = svec(1);
  end;
end;

if (nargout == 2)
  level_labels_out = level_labels;
end;
