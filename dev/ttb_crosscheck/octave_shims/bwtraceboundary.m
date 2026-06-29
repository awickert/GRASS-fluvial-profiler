function C = bwtraceboundary(BW, P, dir, conn, maxlen, orientation)
% Octave shim for MATLAB's bwtraceboundary, backed by bwboundaries.
% Returns the boundary loop (closed) of the region containing P, rotated to
% start at P. Traversal start/direction differ from MATLAB but the undirected
% boundary EDGE SET is identical, which is all the downstream divide use needs.
  if nargin < 4 || isempty(conn), conn = 8; end
  b = bwboundaries(BW, conn, 'noholes');
  best = 1; bestd = inf;
  for k = 1:numel(b)
    L = b{k};
    d = min((L(:,1)-P(1)).^2 + (L(:,2)-P(2)).^2);
    if d < bestd, bestd = d; best = k; end
  end
  L = b{best};
  if size(L,1) > 1 && all(L(1,:) == L(end,:)), L(end,:) = []; end
  [~, i0] = min((L(:,1)-P(1)).^2 + (L(:,2)-P(2)).^2);
  L = [L(i0:end,:); L(1:i0-1,:)];
  C = [L; L(1,:)];
end
