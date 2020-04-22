namespace eval ::hprefstruct:: {
        variable rosetta_bin /misc/chaoyi/src_builds/rosetta/rosetta_src_2018.42.60459_bundle/main/source/bin/rosetta_scripts.linuxgccrelease
        variable rosetta_db /misc/chaoyi/src_builds/rosetta/rosetta_src_2018.42.60459_bundle/main/database/
        variable molprobity_path /home/chaoyi/builds/Molprobity/build/bin
        variable AA_path /misc/chaoyi/src_builds/hpRefStruct/test/building_missing_loops/AAs_pdb
        proc check_install { } {
                if { ![file isfile $::hprefstruct::rosetta_bin] } { error "Cannot find $::hprefstruct::rosetta_bin" }
                if { ![file isdirectory $::hprefstruct::rosetta_db] } { error "Cannot find $::hprefstruct::rosetta_db" }
                if { ![file isdirectory $::hprefstruct::molprobity_path] } { error "Cannot find $::hprefstruct::molprobity_path" }
                if { ![file isdirectory $::hprefstruct::AA_path] } { error "Cannot find $::hprefstruct::AA_path" }
        }
        check_install
}

proc ::hprefstruct::identify_missing { pdb } {
        set ip [open $pdb r]
        set tmp [list]
        while { [gets $ip oline] >= 0 } {
                if { [regexp -- "REMARK 465" $oline] } {
                        if [string is integer -strict [string range $oline 25 25]] {
				if { ![info exists chn] } {
					set chn [string range $oline 19 19]
					set tmp2 [list $chn]
					eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
				} elseif { $chn == [string range $oline 19 19] } {
					eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
				} else {
					lappend tmp $tmp2
					set chn [string range $oline 19 19]
					set tmp2 [list $chn]
					eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
				}
                        }
                }
        }
	lappend tmp $tmp2
        close $ip
        return $tmp
}

proc ::hprefstruct::separate_chain { pdb } {
	set id [mol new $pdb]
	set ary [::hprefstruct::identify_missing $pdb]
	if { [llength $ary] == 0 } {
		error "Cannot find missing residue record in $pdb ..."
	} else {
		foreach i $ary {
			set chn [lindex $i 0]
			set sel [atomselect $id "protein and chain $chn"]
			set list1 [lrange $i 1 end]
			set l2 [::hprefstruct::find_cons [lsort -u -integer [$sel get resid]]]
			set tmp 0
			foreach {i j} $l2 {
				set flna [lindex [split $pdb .] 0].chain-$chn.$tmp.pdb
#				[atomselect $id "protein and chain $chn and resid $i to $j"] writepdb $flna
				# to do: use tcl to realize sed
#				exec sed -i "1d;\$d" $flna
				incr tmp
			}
			::hprefstruct::build_segments $chn $list1
#			::hprefstruct::sort_list [lindex [split $pdb .] 0] $chn $list1 $l2		
		}
	}
	mol delete $id
}

proc ::hprefstruct::build_segments { chn list1 } {
        set list1-1 [list]
        set list1-2 [list]
        foreach {i j} $list1 {
                lappend list1-1 $i
                lappend list1-2 $j
        }
        set l1 [::hprefstruct::find_cons ${list1-1}]
	set count 0
	foreach {i j} $l1 {
                set start [lsearch ${list1-1} $i]
                set end [lsearch ${list1-1} $j]
		for {set k $start} {$k <= $end} {incr k} {
			set pdb "AAs_pdb/[lindex ${list1-2} $k].pdb"
			file copy -force $pdb temp
			set resid [lindex ${list1-1} $k]
	                exec sed -i "s/ A / $chn /g" temp
                        exec sed -i "s/   1   /$resid   /g" temp
                        exec cat temp >> $chn.build.$count.pdb
		}
		incr count
	}
}

proc ::hprefstruct::sort_list { name chn list1 l2 } {
	set list1-1 [list]
        set list1-2 [list]
        foreach {i j} $list1 {
	        lappend list1-1 $i
                lappend list1-2 $j
	}
	set l1 [::hprefstruct::find_cons ${list1-1}]
	set loop [open "loop_chain-$chn.loop" w]
	file delete $chn.pdb
	if { [llength $l1] > [llength $l2] } {
		set start [lsearch ${list1-1} [lindex $l1 0]]
		set end [lsearch ${list1-1} [lindex $l1 1]]
		puts $loop "LOOP [lindex $l1 0] [expr [lindex $l1 1] + 1] [expr [lindex $l1 1] + 1] 0 1"
		foreach {i j} [lrange $l1 2 end-2] {
			puts $loop "LOOP [expr [lindex $l1 $i] - 1] [expr [lindex $l1 $j] + 1] [expr [lindex $l1 $j] + 1] 0 1"
		}
		puts $loop "LOOP [expr [lindex $l1 end-1] - 1] [lindex $l1 end] [expr [lindex $l1 end-1] - 1] 0 1"
		close $loop
		for {set i $start} {$i <= $end} {incr i} {
			set pdb "AAs_pdb/[lindex ${list1-2} $i].pdb"
			file copy -force $pdb temp
			set resid [lindex ${list1-1} $i]
			exec sed -i "s/ A / $chn /g" temp
			exec sed -i "s/   1   /$resid  /g" temp
			exec cat temp >> $chn.pdb
		}
		set l1_ [lrange $l1 2 end]
		set k 0
		foreach {i j} ${l1_} {
			exec cat $name.chain-$chn.$k.pdb >> $chn.pdb
			set start [lsearch ${list1-1} $i]
			set end [lsearch ${list1-1} $j]
			for {set l $start} {$l <= $end} {incr l} {
				set pdb "AAs_pdb/[lindex ${list1-2} $l].pdb"
				file copy -force $pdb temp
				set resid [lindex ${list1-1} $l]
                        	exec sed -i "s/ A / $chn /g" temp
                        	exec sed -i "s/   1   /  $resid   /g" temp
                        	exec cat temp >> $chn.pdb
			}
			incr k
		}
		
	}
	set id [mol new $chn.pdb]
	[atomselect $id all] writepdb $chn.pdb
	mol delete $id
}


proc ::hprefstruct::find_cons { list } {
	set tmp [lindex $list 0]
	set x [lindex $list 0]
	for { set i 1 } { $i < [llength $list] } { incr i } {
		if { [lindex $list $i] == [expr $x + 1] } {
			incr x
			continue
		} else {
			lappend tmp [lindex $list $i-1] [lindex $list $i]
			set x [lindex $list $i]
		}
	}
	lappend tmp [lindex $list end]
	return $tmp
}

proc ::hprefstruct::get_seq { pdb chain } {
        set ip [open $pdb r]
        set tmp [list]
        while { [gets $ip oline] >= 0 } {
                if { [regexp -- "REMARK 465" $oline] } {
                        if [string is integer -strict [string range $oline 25 25]] {
                                if { ![info exists chn] } {
                                        set chn [string range $oline 19 19]
                                        set tmp2 [list $chn]
                                        eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
                                } elseif { $chn == [string range $oline 19 19] } {
                                        eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
                                } else {
                                        lappend tmp $tmp2
                                        set chn [string range $oline 19 19]
                                        set tmp2 [list $chn]
                                        eval lappend tmp2 [string range $oline 22 25] [string range $oline 15 17]
                                }
                        }
                }
        }
}

proc ::hprefstruct::get_chain { pdb } {
	set id [mol new $pdb]
        set chn [lsort -u [[atomselect $id protein] get chain]]
	mol delete $id
	return $chn
}

proc ::hprefstruct::edit_chain { chain miss } {
	variable AA_path
	
}

set x [::hprefstruct::identify_missing 6vsb.pdb]
puts [lindex $x 0]
#puts [lindex $x 1]
#puts [lindex $x 2]
::hprefstruct::separate_chain 6vsb.pdb


