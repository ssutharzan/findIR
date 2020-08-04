% The following piece of code will generate a data set of random nucleotide
% sequences and estimate the average runtimes of the function 'findIRSeq' on the
% data set.


% Generating the data set
seqSet = cell(100,100);
for i=1:100
    for j=1:100
        seqSet{i,j} = randseq(1000*i);
    end
end

tm = cell(100,100);


% Running the function on the data set
for k=1:100
    for l=1:100
        tic;
        [iRStartPos,iREndPos,iRSeq] = findIRSeq(seqSet{k,l},4,1000);
        tm{k,l}=toc;
        clear iRStartPos iREndPos iRSeq;
    end
end


% Estimating the average runtimes and the standarad devaiations
mean_time = zeros(1,100);
std_time = zeros(1,100);

for m=1:100
    mean_time(m) = mean(cell2mat(tm(m,:)));
    std_time(m) = std(cell2mat(tm(m,:)));
end


% Plotting the results
errorbar(mean_time,std_time)
xlabel('Sequence Length (x1000 nucleotides)');
ylabel('Average Runtitme (s)');