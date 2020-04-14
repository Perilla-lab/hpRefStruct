proc rama { pdb } {
	puts "Ramalyze $pdb ..."
	set path /home/chaoyi/builds/Molprobity/build/bin
	set op [open "summary" a]
	puts $op "Ramalyze results:"
	close $op
	exec -ignorestderr -- $path/phenix.ramalyze ${pdb} outliers_only=True >> summary
	puts "Finished Ramalyze $pdb"
	set ip [open "summary" r]
	set temp1 [list]
    	while {[gets $ip oline] >= 0} {
        	if { [regexp -- "Ramalyze results" $oline] } {
			while { [gets $ip iline] >= 0 } {
                		if [string is integer -strict [string range $iline 2 5]] {
                    			lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                		}
			}
		}
	}
	unset -nocomplain path oline iline
	close $ip
	return [join [lsort -u -dictionary $temp1]]
}

proc rota { pdb } {
        puts "Rotalyze $pdb ..."
        set path /home/chaoyi/builds/Molprobity/build/bin
        set op [open "summary" a]
        puts $op "Rotalyze results:"
        close $op
        exec -ignorestderr -- $path/phenix.rotalyze ${pdb} outliers_only=True >> summary
        puts "Finished Rotalyze $pdb"
        set ip [open "summary" r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Rotalyze results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if [string is integer -strict [string range $iline 2 5]] {
                                        lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                                }
                        }
                }
        }
        unset -nocomplain path oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc cbeta { pdb } {
        puts "Cbetadev $pdb ..."
        set path /home/chaoyi/builds/Molprobity/build/bin
        set op [open "summary" a]
        puts $op "Cbetadev results:"
        close $op
        exec -ignorestderr -- $path/phenix.cbetadev ${pdb} outliers_only=True >> summary
        puts "Finished Cbetadev $pdb"
        set ip [open "summary" r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Cbetadev results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if [string is integer -strict [string range $iline 2 5]] {
                                        lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                                }
                        }
                }
        }
        unset -nocomplain path oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc clash { pdb } {
        puts "Clashscore $pdb ..."
        set path /home/chaoyi/builds/Molprobity/build/bin
        set op [open "summary" a]
        puts $op "Clashscore results:"
        close $op
        exec -ignorestderr -- $path/phenix.clashscore ${pdb} keep_hydrogens=True >> summary
        puts "Finished Clashscore $pdb"
        set ip [open "summary" r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Clashscore results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if [string is integer -strict [string range $iline 2 5]] {
                                        lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
					lappend temp1 "[string range $iline 18 18] [string range $iline 19 22]"

                                }
                        }
                }
        }
        unset -nocomplain path oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc omega { pdb } {
        puts "Omegalyze $pdb ..."
        set path /home/chaoyi/builds/Molprobity/build/bin
        set op [open "summary" a]
        puts $op "Omegalyze results:"
        close $op
        exec -ignorestderr -- $path/phenix.omegalyze ${pdb}  >> summary
        puts "Finished Omegalyze $pdb"
        set ip [open "summary" r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Omegalyze results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if [string is integer -strict [string range $iline 2 5]] {
                                        lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                                        lappend temp1 "[string range $iline 18 18] [string range $iline 19 22]"

                                }
                        }
                }
        }
        unset -nocomplain path oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

