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

	} else {
		foreach i $ary {
			set chn [lindex $i 0]
			set sel [atomselect $id "protein and chain $chn"]
			set list1 [lrange $i 1 end]
			set list1-1 [list]
			set list1-2 [list]
			foreach i j $list1 {
				lappend list1-1 $i
				lappend list1-2 $j
			}
			set list2 [lsort -u -integer [$sel get resid]]
			set tmp 0
			foreach i j [::hprefstruct::find_cons $list2] {
				[atomselect $id "protein and chain $chn and resid $i to $j"] writepdb [lindex [split $pdb .] 0].chain-$chn.$tmp.pdb
				exec sed -i '1d;$d' [lindex [split $pdb .] 0].chain-$chn.$tmp.pdb
				incr tmp
			}
			foreach a ${list1-1} b $list2 {
				
		}
	}
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

#set x [::hprefstruct::identify_missing 6vsb.pdb]
#puts [lindex $x 0]
#puts [lindex $x 1]
#puts [lindex $x 2]
#separate_chain 6vsb.pdb


