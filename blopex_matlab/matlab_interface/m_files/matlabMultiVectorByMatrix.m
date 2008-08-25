function blockVectorZ = matlabMultiVectorByMatrix(blockVectorX,maskX,...
    blockVectorY,maskY, matrix)

blockVectorY(:,maskY) = blockVectorX(:,maskX)*matrix;
blockVectorZ = blockVectorY;
