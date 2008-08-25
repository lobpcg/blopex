function blockVectorZ = matlabMultiVectorByDiag(blockVectorX,maskX,...
    blockVectorY,maskY, values)

matrix = diag(values);
blockVectorY(:,maskY) = blockVectorX(:,maskX)*matrix;
blockVectorZ = blockVectorY;
