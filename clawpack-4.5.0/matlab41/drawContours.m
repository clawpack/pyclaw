function h = drawContours(c,sval,sdir)

  st_idx = 1;
  line_cnt = 0;
  if (isempty(c))
    % No contour lines on this patch
    h = [];
    return;
  end;

  while (1)
    cval = c(1,st_idx);
    next_length = c(2,st_idx);
    line_cnt = line_cnt + 1;
    xdata = c(1,st_idx+1:st_idx+next_length);
    ydata = c(2,st_idx+1:st_idx+next_length);
    zdata = zeros(1,next_length) + sval;
    if (strcmp(sdir,'x'))
      % (y,z) passed in as (x,y)
      h(line_cnt) = line('XData',zdata,'YData',xdata,'ZData',ydata,'Color','k');
    elseif (strcmp(sdir,'y'))
      % (x,z) passed in as (x,y)
      h(line_cnt) = line('XData',xdata,'YData',zdata,'ZData',ydata,'Color','k');
    elseif (strcmp(sdir,'z'))
      % (x,y) passed in as (x,y)
      h(line_cnt) = line('XData',xdata,'YData',ydata,'ZData',zdata,'Color','k');
    end;
    set(h(line_cnt),'Tag','contourline');
    set(h(line_cnt),'UserData',struct('sval',sval,'sdir',sdir,'cval',cval));
    st_idx = st_idx + next_length + 1;
    if (st_idx > length(c))
      return;
    end;
  end;
