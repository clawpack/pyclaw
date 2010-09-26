
slicest = num2str(next_slice);
while(size(slicest,2))<5,
  slicest = ['0',slicest];
  end;
fname = ['slice' slicest];
printgif(fname);
