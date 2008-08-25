function blockVectorZ = matlabMultiAxpy(blockVectorX,maskX,...
    blockVectorY,maskY,alpha)

if any(size(alpha)~=1)
    error('alpha must be scalar')
end
blockVectorY(:,maskY) = blockVectorY(:,maskY) + ...
  blockVectorX(:,maskX)*alpha;
blockVectorZ = blockVectorY;
