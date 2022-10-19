proc moveby {sel offset} {
  foreach coord [$sel get {x y z}] {
    lappend newcoords [vecadd $coord $offset]
  }
  $sel set {x y z} $newcoords
}

set mol1 [atomselect 1 "all"]
moveby $mol1 {0 0 23}
