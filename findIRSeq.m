function [startPosVec,endPosVec,seqVec] = findIRSeq(inputSeq,minLength,maxLength)

% This function finds all the nested and overlapping inverted repeats in the given input file.
%
%
% Inputs
%
% seq           -   The char variable containing the input sequence.
%                   The sequence should contain only one sequence.
%                   Sequence can contain valid upper and lower case nucleotide 
%                   characters. Valid charcters are A,T,G,C,R,Y,S,W,K,M,B,D,H,V and N. 
%
% minLength 	- 	The minimum legnth of the inverted repeat to be detected,
%					in double. The lowest value should be 4 and the maximum value should be 1000.
%
% maxLength 	- 	The maximum legnth of the inverted repeat to be detected,
%					in double. The lowest value should be 4 and the maximum value should be 1000.
%
%
% Return variables 
%
% startPosVec   -   The vector containg the starting positions of the dectected inverted repeats,
%                   as double.
%
% endPosVec     -   The vector containg the ending positions of the dectected inverted repeats,
%                   as double.
%
% seqVec        -   The cell containg the detected inverted repeats sequences, as char.
%
%
% Example:
%
%    [iRstartPos,iRendPos,iRSeq] = findIRSeq('TAIR10_chr1.fas',4,1000)


% The following piece of code calculates the cumulative scores. Each nucleotype type is given a
% large prime number score to increase the uniqueness of the type.

if(minLength < 2)
    error('The minimum Inverted Repeat legnth to be detected should not be less than 2 nucleotides.');
elseif(maxLength > 1000)
    error('The maximum Inverted Repeat legnth to be detected should not be more than 1000 nucleotides.');
end

seq=lower(inputSeq);
seqSize = size(seq,2);
seedVec = zeros(1,seqSize);
scorePosVec = zeros(1,seqSize);

posScoreVec = [10007,10009,10037]; % posScoreVec - The vector containing the prime number positional scores.
                                   % The runtitme can be optimised by changing these scores.

for base=1:seqSize
    switch(seq(base))
        case 'a'
            scorePosVec(base+1) = posScoreVec(1);
        case 't'
            scorePosVec(base+1) = -posScoreVec(1);
        case 'g'
            scorePosVec(base+1) = posScoreVec(2);
        case 'c'
            scorePosVec(base+1) = -posScoreVec(2);
        case {'r','y','s','w','k','m','b','d','h','v','n'}
            scorePosVec(base+1) = posScoreVec(3);
        otherwise
            error('Input sequence contains invalid characters.');
    end
end

scoreCumVec = cumsum(scorePosVec);
clear scorePosVec;
clear seq;


% Sorts the scores and creates some needed matrices

[scoreCumVecSorted,scoreCumVecSortedIdx] = sort(scoreCumVec);
scoreCumVecSortedDiff = diff(scoreCumVecSorted);
sizeScoreCumVec = size(scoreCumVec,2)-1;


% The following piece of code finds unique scores which occur more than once.

pos = 0;

if(scoreCumVecSortedDiff(1) == 0 && scoreCumVecSortedDiff(2) ~= 0)
    pos = pos+1;
    uniqueCumScoreVec(pos) = scoreCumVec(scoreCumVecSortedIdx(1));
end

for score=2:sizeScoreCumVec
    if(score < sizeScoreCumVec && scoreCumVecSortedDiff(score) == 0 && scoreCumVecSortedDiff(score-1) ~= 0 && scoreCumVecSortedDiff(score + 1) ~= 0)
        if((scoreCumVecSortedIdx(score + 1) - scoreCumVecSortedIdx(score)) > 1 && (scoreCumVecSortedIdx(score + 1) - scoreCumVecSortedIdx(score)) < maxLength + 1)
            pos = pos +1;
            uniqueCumScoreVec(pos) = scoreCumVec(scoreCumVecSortedIdx(score));
        end
    elseif (scoreCumVecSortedDiff(score) == 0 && scoreCumVecSortedDiff(score -1) ~= 0)
        pos = pos + 1;
        uniqueCumScoreVec(pos) = scoreCumVec(scoreCumVecSortedIdx(score));
    elseif (score < sizeScoreCumVec && scoreCumVecSortedDiff(score) == 0 && scoreCumVecSortedDiff(score + 1) ~= 0)
        pos = pos + 1;
        uniqueCumScoreVec(pos) = scoreCumVec(scoreCumVecSortedIdx(score));
    end
end


% The following piece of code searches for the occurances of above found unique scores

if (pos > 0)
    iRIniCount = 0;
    sortedCumScorePos = 1;
    clear scoreCumVec;
    sizeUniqueCumScoreVec = size(uniqueCumScoreVec,2);
    initialIRSet = cell(1,sizeUniqueCumScoreVec);
    
    for i=1:sizeUniqueCumScoreVec
        for j=sortedCumScorePos:size(scoreCumVecSorted,2)
            if(scoreCumVecSorted(j) == uniqueCumScoreVec(i))
                sortedCumScorePos = j;
                break;
            end
        end
        sizeScoreCumVecSorted = size(scoreCumVecSorted,2);
        for k=sortedCumScorePos:sizeScoreCumVecSorted
            if(scoreCumVecSorted(k)>uniqueCumScoreVec(i))
                lastSortedCumScorePos = k-1;
                break;
            end
            if(k == size(scoreCumVecSorted,2) && scoreCumVecSorted(k) == uniqueCumScoreVec(i))
                lastSortedCumScorePos=k;
            end
        end
        sizeInitialIRset = lastSortedCumScorePos-sortedCumScorePos+1;
        
        for l=1:sizeInitialIRset
            initialIRSet{i}(l) = scoreCumVecSortedIdx(sortedCumScorePos + l - 1);
        end
        sortedCumScorePos = lastSortedCumScorePos;
    end
    
    
% The following piece of code constructs the possible inverted repeats and their seed points.
    
clear scoreCumVecSorted;
startPosVecIni = zeros(1,1000000);
endPosVecIni = zeros(1,1000000);
lenVecIni = zeros(1,1000000);

for m=1: sizeUniqueCumScoreVec
    iRPosVec = initialIRSet{m};
    if(size(iRPosVec,2) > 1)
        sizeIRPosVecLoopSize = size(iRPosVec,2)-1;
        for n=1:sizeIRPosVecLoopSize
            for j=n+1:size(iRPosVec,2)
                startPos = iRPosVec(n);
                endPos = iRPosVec(j) - 1;
                if  (endPos-startPos) < maxLength && mod((endPos-startPos+1),2) ==0
                    iRIniCount = iRIniCount+1;
                    startPosVecIni(iRIniCount) = startPos;
                    endPosVecIni(iRIniCount) = endPos;
                    lenVecIni(iRIniCount) = ((endPos-startPos + 1) / 2);
                    seedPosIni=startPos+((endPos-startPos + 1) / 2);
                    seedVec(seedPosIni) = seedVec(seedPosIni) + 1;
                else
                    break;
                end
            end
        end
    end
end

startPosVecIni = startPosVecIni(1:iRIniCount);
endPosVecIni = endPosVecIni(1:iRIniCount);
lenVecIni = lenVecIni(1:iRIniCount);

% The following piece of code validations the possible inverted repeats and finds the valid ones.

clear scoreCumVecSortedDiff;
clear pos_f_mat;
clear z_f_mat;

iRCount = 0;
lenSeqMat = endPosVecIni-startPosVecIni+1;

[~,lenSeqMatIdx] = sort(lenSeqMat,'descend');

startPosVec = zeros(1,size(lenSeqMatIdx,2));
endPosVec = zeros(1,size(lenSeqMatIdx,2));
sizeLenSeqMatIdx = size(lenSeqMatIdx,2);


for iniIR=1:sizeLenSeqMatIdx
    lenSeqMatCount = lenSeqMatIdx(iniIR);
    if(lenSeqMat(lenSeqMatCount) >= minLength)
        seedPosIni=startPosVecIni(lenSeqMatCount) + ((endPosVecIni(lenSeqMatCount) - startPosVecIni(lenSeqMatCount) + 1) / 2);
        if (lenVecIni(lenSeqMatCount) == seedVec(seedPosIni))
            iRCount = iRCount+1;
            startPosVec(iRCount) = startPosVecIni(lenSeqMatCount);
            endPosVec(iRCount) = endPosVecIni(lenSeqMatCount);
        end
        seedVec(seedPosIni) = seedVec(seedPosIni) - 1;
    end
end
startPosVec = startPosVec(1:iRCount);
endPosVec = endPosVec(1:iRCount);
else
    iRCount = 0;
end


% The following piece of code stores the valids inverted repeat Sequences.

if(iRCount > 0)
    seqVec = cell(1,iRCount);
    
    for iRPos=1:iRCount
        seqVec{iRPos} = inputSeq(startPosVec(iRPos):endPosVec(iRPos));
    end
else
    startPosVec = [];
    endPosVec = [];
    seqVec = {};
end
