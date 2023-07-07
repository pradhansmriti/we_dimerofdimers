mol load pdb cg_ABCD.pdb
set AB [atomselect top "segid P1"]
#set CD [atomselect top "segid P2"]
$AB moveby {-50 -50 10}
set all [atomselect top all]
$all writepdb basis.pdb
