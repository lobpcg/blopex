function blockVectorZ = matlabCopyMultiVector(blockVectorX,maskX,...
    blockVectorY,maskY)

blockVectorY(:,maskY) = blockVectorX(:,maskX);
blockVectorZ = blockVectorY;
