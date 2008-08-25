function values = matlabMultiInnerProdDiag(blockVectorX,maskX,...
    blockVectorY,maskY)

values = dot(blockVectorX(:,maskX),blockVectorY(:,maskY))';
