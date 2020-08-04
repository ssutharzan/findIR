% The following piece of code calls the 'findIR' function. The runtime of
% the function will be stored in the variable 'tm'.
tic
[iRStartPos,iREndPos,iRSeq]=findIR('atChr1.fa',4,1000);
tm=toc;
