
% compare coarse grid solution with finer grid solution in subdirectory qref

setplot1

if length(mq)>1
   disp('can only compare one component at a time, set mq to a single value')
   mq = input('  mq = ');
   clf
   end

Frame = 0;
while Frame <= MaxFrames

   cd qref
   plotstyle = 'b-';
   plotframe1
   qref = q;
   mxref = mx;
   xref = x;
   cd ..

   hold on
   plotstyle = 'ro';
   plotframe1
   hold off

   ratio = mxref/mx;
   qref1 = coarsen(qref,ratio);
   err1 = sum(abs(q-qref1)) * dx

   queryframe
   end
