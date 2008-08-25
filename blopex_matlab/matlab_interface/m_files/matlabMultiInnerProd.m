function matrix = matlabMultiInnerProd (blockVectorX,maskX,...
    blockVectorY,maskY)
    

% if size(blockVector1In,1)~=size(blockVector1In,1)
%     error('multivector heights are different')
% end

matrix = blockVectorX(:,maskX)'*blockVectorY(:,maskY);
