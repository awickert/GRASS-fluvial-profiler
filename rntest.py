import RiverNetwork as rn
reload(rn)
import RiverNetwork as rn

seg1 = rn.Segment(1, to_id=3)
seg2 = rn.Segment(2, to_id=3)
seg3 = rn.Segment(3)

segs = [seg1, seg2, seg3]

net = rn.Network(segs)

net.compute_fromstreams()
