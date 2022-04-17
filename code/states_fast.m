statMapP = zeros(8);
statMapN = zeros(8);

for iRow = 1:size(ethogram.behmat,1);
    for iCol = 1:size(ethogram.behmat,2)-1
        if all(~isnan(ethogram.behmat(iRow,iCol:iCol+1)))
            s1 = ethogram.behmat(iRow,iCol);
            s2 = ethogram.behmat(iRow,iCol+1);
            if ethogram.pulmat(iRow,iCol);
                statMapP(s1,s2) = statMapP(s1,s2)+1;
            else
                statMapN(s1,s2) = statMapN(s1,s2)+1;
            end
        end
    end
end
%%
figure(10);clf
imagesc(statMapP)
normSMP = bsxfun(@rdivide,statMapP,sum(statMapP));
normSMN = bsxfun(@rdivide,statMapN,sum(statMapN));

figure(11);clf
%normSMP = statMapP./sum(statMapP(:));
%normSMN = statMapN./sum(statMapN(:));
imm = [0 max(max(cat(1,normSMP,normSMN)))];
imagesc(normSMP,imm),colorbar
figure(12);clf
imagesc(normSMN,imm),colorbar


figure(13);clf
imagesc((normSMP-normSMN)*100),colorbar

