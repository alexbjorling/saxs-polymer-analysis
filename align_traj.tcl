proc align_traj_alts {mol1 sel1 mol2 sel2 mol3 sel3} {
# align_traj_gro_alts mol1 sel1 mol2 sel2 mol3 sel3 moves molecule 1 onto one of the single-frame molecules 2 or 3, whichever fits best at each frame, throughout molecule 1's trajectory, based on a comparison of the atoms specified as selection expressions.  
    set reference [atomselect $mol1 $sel1]
    set compare1 [atomselect $mol2 $sel2]
    set compare2 [atomselect $mol3 $sel3]
    set num_steps [molinfo $mol1 get numframes]
    for {set iframe 0} {$iframe < $num_steps} {incr iframe} {
            $reference frame $iframe
# decide which fits best
# ... first move onto compare1 and calculate rmsd.
            set trans_mat_A [measure fit $reference $compare1]
            set movemol [atomselect $mol1 all frame $iframe]
            $movemol move $trans_mat_A
            set rmsd_A [measure rmsd $reference $compare1]
# ... now move onto compare2 and calculate rmsd.
            set trans_mat_B [measure fit $reference $compare2]
            set movemol [atomselect $mol1 all frame $iframe]
            $movemol move $trans_mat_B
            set rmsd_B [measure rmsd $reference $compare2]
# ... now either leave or move back
            if {$rmsd_A<$rmsd_B} {
                set trans_mat [measure fit $reference $compare1]
                set movemol [atomselect $mol1 all frame $iframe]
                $movemol move $trans_mat 
            }
    }
}

proc align_traj {mol1 sel1 mol2 sel2} {
# align_traj mol1 sel1 mol2 sel2 moves molecule 1 onto molecule 2 throughout their trajectories. If one trajectory runs out of frames it uses the last one.
# for example,
# > align_traj 0 "name BB" 1 "name CA"
# aligns a fg and a cg structure, based on backbone/C-alpha fitting.
    if {$sel2 == ""} {set sel2 $sel1}
    set group1 [atomselect $mol1 $sel1]
    set group2 [atomselect $mol2 $sel2]
    set num_steps_1 [molinfo $mol1 get numframes]
    set num_steps_2 [molinfo $mol2 get numframes]
    set num_steps_max $num_steps_1
#    echo $num_steps_1 $num_steps_2
    if {$num_steps_2 > $num_steps_1} {set num_steps_max $num_steps_2}
    for {set iframe 0} {$iframe < $num_steps_max} {incr iframe} {
        if {$iframe <= $num_steps_1-1} {$group1 frame $iframe}
        if {$iframe <= $num_steps_2-1} {$group2 frame $iframe}
#	echo [$group1 frame] [$group2 frame]
        # compute the transformation
        set trans_mat [measure fit $group1 $group2]
        # do the alignment
        set movemol [atomselect $mol1 all frame $iframe]
        $movemol move $trans_mat
    }
}

proc align_traj_on_itself {mol sel ref_frame} {
# align_traj_on_itself mol sel ref_frame aligns a molecule on the specified frame of its trajectory.
    set group [atomselect $mol $sel]
    set ref [atomselect $mol $sel]
    $ref frame $ref_frame
    set num_steps [molinfo $mol get numframes]
    for {set iframe 0} {$iframe < $num_steps} {incr iframe} {
	$group frame $iframe
        # compute the transformation
        set trans_mat [measure fit $group $ref]
        # do the alignment
        set movemol [atomselect $mol all frame $iframe]
        $movemol move $trans_mat
    }
}

