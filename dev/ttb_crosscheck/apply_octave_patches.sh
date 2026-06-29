#!/usr/bin/env bash
# Make a TopoToolbox (v2) clone runnable under GNU Octave for the DIVIDEobj
# cross-check. Applies the minimal source edits that work around Octave/MATLAB
# incompatibilities. Each edit is a compatibility fix, NOT an algorithm change;
# see the rationale inline and in README.md.
#
# TopoToolbox (GPL) by Wolfgang Schwanghart and Dirk Scherler:
#   https://github.com/wschwanghart/topotoolbox
# We do not redistribute TopoToolbox; clone it yourself and point TTB at it.
#
#   git clone --depth 1 https://github.com/wschwanghart/topotoolbox.git /tmp/topotoolbox
#   TTB=/tmp/topotoolbox bash apply_octave_patches.sh
set -euo pipefail
TTB="${TTB:-/tmp/topotoolbox}"
echo "Patching TopoToolbox clone at: $TTB"

python3 - "$TTB" <<'PY'
import sys, re
ttb = sys.argv[1]

# 1) polygon2GRIDobj: force the non-polyshape (poly2mask) path. Octave lacks the
#    polyshape class; the else-branch rasterises with poly2mask, which Octave has.
p = ttb + '/@GRIDobj/polygon2GRIDobj.m'; s = open(p).read()
s = s.replace("if ~verLessThan('matlab','9.3')",
              "if false  % Octave: force non-polyshape poly2mask path")
open(p, 'w').write(s)

# 2) unique(): Octave does not implement the 3rd output (inverse index map ic)
#    for 'stable'/'rows'. Compute it via ismember (exact MATLAB values).
p = ttb + '/@DIVIDEobj/DIVIDEobj.m'; s = open(p).read()
s = s.replace("[e,~,f] = unique([T0(ix,1);T0(ix,2)],'rows','stable');",
              "tmpV=[T0(ix,1);T0(ix,2)]; e=unique(tmpV,'stable'); [~,f]=ismember(tmpV,e);")
s = s.replace("[e,~,f] = unique(st(:),'stable');",
              "e=unique(st(:),'stable'); [~,f]=ismember(st(:),e);")
open(p, 'w').write(s)

p = ttb + '/@GRIDobj/GRIDobj2polygon.m'; s = open(p).read()
s = s.replace("[uniquevals,~,DB2.Z(I)] = unique(DB.Z(I));",
              "uniquevals=unique(DB.Z(I)); [~,tmpic]=ismember(DB.Z(I),uniquevals); DB2.Z(I)=tmpic;")
open(p, 'w').write(s)

# 3) divorder.m uses 'do' as a variable/field name; 'do' is a reserved keyword in
#    Octave (do-until). Rename the variable/field to 'dvo' (no algorithm change).
p = ttb + '/@DIVIDEobj/divorder.m'; s = open(p).read()
s = re.sub(r'\bdo\b', 'dvo', s)
open(p, 'w').write(s)
print("patched: polygon2GRIDobj, DIVIDEobj (unique x2), GRIDobj2polygon (unique), divorder (do->dvo)")
PY
echo "Done. Also addpath the octave_shims/ dir (bwtraceboundary, verLessThan) BEFORE TopoToolbox."
