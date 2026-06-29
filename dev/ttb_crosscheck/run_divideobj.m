% Run TopoToolbox DIVIDEobj (Scherler & Schwanghart, 2020) on flow routing
% supplied by the r.fluvial / dreich pipeline, and export the divide network for
% comparison against rivernetworkx.divides. See README.md.
%
% TopoToolbox (GPL) by Wolfgang Schwanghart and Dirk Scherler:
%   https://github.com/wschwanghart/topotoolbox
% This script does NOT vendor TopoToolbox; it loads a separately-obtained clone
% that has been made Octave-compatible by apply_octave_patches.sh.
%
% Edit TTB and SHIMS below to your paths, then:
%   octave --no-gui --quiet run_divideobj.m
more off;
TTB   = getenv('TTB');   if isempty(TTB),   TTB   = '/tmp/topotoolbox'; end
SHIMS = getenv('SHIMS'); if isempty(SHIMS), SHIMS = fullfile(pwd, 'octave_shims'); end
addpath(SHIMS); addpath(genpath(TTB));   % shims first so they take precedence
pkg load image statistics;

S = load('/tmp/ttb_M.mat');              % giver->receiver matrix from build_M.py
M = S.M; nr = double(S.nr); nc = double(S.nc); cs = double(S.cs); thr = double(S.thresh);
printf('M %dx%d, grid %dx%d, cs=%g, thresh=%g\n', size(M,1), size(M,2), nr, nc, cs, thr);

refmat = [0 -cs; cs 0; -cs cs];          % any consistent north-up refmat (grid-space compare)
FD = FLOWobj(M, 'cellsize', cs, 'size', [nr nc], 'refmat', refmat);   % routing == dreich's
A  = flowacc(FD);
ST = STREAMobj(FD, A >= thr);
D  = DIVIDEobj(FD, ST, 'outlets', false);   % outlets=false: edge-outlet getdivide is fragile on a clip
D  = divorder(D, 'topo');
printf('DIVIDEobj: %d nodes, ep=%d, jct=%d, maxorder=%g\n', ...
       numel(D.IX), numel(D.ep), numel(D.jct), max(D.order));

dlmwrite('/tmp/ttb_out_IX.txt',    D.IX(:),      'precision', '%d');
dlmwrite('/tmp/ttb_out_order.txt', D.order(:),   'precision', '%.10g');
dlmwrite('/tmp/ttb_out_meta.txt',  [nr nc cs],   'precision', '%g');
printf('EXPORTED %d divide nodes\n', numel(D.IX));
