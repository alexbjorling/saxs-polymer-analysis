source ../../../code/align_traj.tcl
align_traj 0 "all" 1 "all"

proc moveby {sel offset} {
  foreach coord [$sel get {x y z}] {
    lappend newcoords [vecadd $coord $offset]
  }
  $sel set {x y z} $newcoords
}

set mol1 [atomselect 1 "all"]
moveby $mol1 {20 0 0}

render snapshot render.png

