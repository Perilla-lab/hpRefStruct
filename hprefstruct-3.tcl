#############################################################################
# Protein structure refinement based on Cryo-EM maps and Molprobity scores
#############################################################################
#     Functions
#############################################################################

namespace eval ::hprefstruct:: {
	variable rosetta_bin /home_orig/chaoyi/builds/rosetta_src_2020.08.61146_bundle/main/source/bin/rosetta_scripts.linuxgccrelease
	variable rosetta_db /home_orig/chaoyi/builds/rosetta_src_2020.08.61146_bundle/main/database/
	variable molprobity_path /home/chaoyi/builds/Molprobity/build/bin
	proc check_install { } {
		if { ![file isfile $::hprefstruct::rosetta_bin] } { error "Cannot find $::hprefstruct::rosetta_bin" }
		if { ![file isdirectory $::hprefstruct::rosetta_db] } { error "Cannot find $::hprefstruct::rosetta_db" }
		if { ![file isdirectory $::hprefstruct::molprobity_path] } { error "Cannot find $::hprefstruct::molprobity_path" }
	}
	check_install
}

proc hprs { args } {
	eval ::hprefstruct::hprs $args 
}

proc ::hprefstruct::rama { pdb } {
	variable molprobity_path
	variable report
	puts "Ramalyze $pdb ..."
	set op [open $report a]
	puts $op "\n##############################################"
	puts $op "Ramalyze results:"
	close $op
	exec -ignorestderr -- ${molprobity_path}/phenix.ramalyze ${pdb} outliers_only=True >> $report
	puts "Finished Ramalyze $pdb"
	set ip [open $report r]
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
	unset -nocomplain oline iline
	close $ip
	return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::rama2 { pdb } {
        variable molprobity_path
        variable report
        puts "Ramalyze $pdb ..."
        set op [open $report a]
        puts $op "\n##############################################"
        puts $op "Ramalyze results:"
        close $op
        exec -ignorestderr -- ${molprobity_path}/phenix.ramalyze ${pdb} >> $report
        puts "Finished Ramalyze $pdb"
        set ip [open $report r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Ramalyze results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if { [string is integer -strict [string range $iline 2 5]] && ([lindex [split [lindex $iline 2] :] 4] == "Allowed" || [lindex [split [lindex $iline 2] :] 4] == "OUTLIER" ) } {
                                       lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                                }
                        }
                }
        }
        unset -nocomplain oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}


proc ::hprefstruct::rota { pdb } {
	variable molprobity_path
	variable report
        puts "Rotalyze $pdb ..."
        set op [open $report a]
	puts $op "\n##############################################"
        puts $op "Rotalyze results:"
        close $op
        exec -ignorestderr -- ${molprobity_path}/phenix.rotalyze ${pdb} outliers_only=True >> $report
        puts "Finished Rotalyze $pdb"
        set ip [open $report r]
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
        unset -nocomplain oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::cbeta { pdb } {
	variable molprobity_path
	variable report
        puts "Cbetadev $pdb ..."
        set op [open $report a]
	puts $op "\n##############################################"
        puts $op "Cbetadev results:"
        close $op
        exec -ignorestderr -- ${molprobity_path}/phenix.cbetadev ${pdb} outliers_only=True >> $report
        puts "Finished Cbetadev $pdb"
        set ip [open $report r]
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
        unset -nocomplain oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::clash { pdb } {
	variable molprobity_path
	variable report
        puts "Clashscore $pdb ..."
        set op [open $report a]
	puts $op "\n##############################################"
        puts $op "Clashscore results:"
        close $op
        exec -ignorestderr -- ${molprobity_path}/phenix.clashscore ${pdb} keep_hydrogens=True >> $report
        puts "Finished Clashscore $pdb"
        set ip [open $report r]
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
        unset -nocomplain oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::omega { pdb } {
	variable molprobity_path
	variable report
        puts "Omegalyze $pdb ..."
        set op [open $report a]
	puts $op "\n##############################################"
        puts $op "Omegalyze results:"
        close $op
        exec -ignorestderr -- ${molprobity_path}/phenix.omegalyze ${pdb}  >> $report
        puts "Finished Omegalyze $pdb"
        set ip [open $report r]
        set temp1 [list]
        while {[gets $ip oline] >= 0} {
                if { [regexp -- "Omegalyze results" $oline] } {
                        while { [gets $ip iline] >= 0 } {
                                if [string is integer -strict [string range $iline 2 5]] {
                                        lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                                        lappend temp1 "[string range $iline 16 16] [string range $iline 17 20]"
                                }
                        }
                }
        }
        unset -nocomplain oline iline
        close $ip
        return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::combine { list1 } {
	set temp1 [list]
 	foreach x $list1 {
		foreach {chn id} $x {
			lappend temp1 "${chn} ${id}"
		}
	}
    	unset -nocomplain  list1 list2
    	return [join [lsort -u -dictionary $temp1]]
}

proc ::hprefstruct::structural_analysis { pdb } {
	set name [lindex [split $pdb .] 0]
	variable iter 
	variable report summary.$name.iter$iter.rep
	set ip [open $report w]
	puts $ip "##############################################"
	puts $ip "#Summary for structure $pdb, iteration $iter"
	puts $ip "##############################################\n"
	close $ip

	set ll [list [::hprefstruct::rama2 $pdb] [::hprefstruct::clash $pdb] [::hprefstruct::rota $pdb] [::hprefstruct::cbeta $pdb] [::hprefstruct::omega $pdb]]
	set x [::hprefstruct::combine $ll]
	set ip [open $report a]
	puts $ip "\n##############################################"
	puts $ip "The list of problematic residues identified in iteration $iter are"
	puts $ip "Chain	Residue"
	foreach {chn id} [join $x] {
		puts $ip "$chn	$id"
	}
	puts $ip "##############################################"
	close $ip
	unset -nocomplain ll
	return $x
}

proc ::hprefstruct::greppdb { jobname } {
	variable iter
    set input [open ${jobname}_run[expr ${iter} - 1]_0001.pdb r]
    set output [open ${jobname}.pdb w]
    while {[gets $input line] >= 0} {
        if {[lindex $line 0] == "ATOM"} {
            puts $output $line
        } elseif {[lindex $line 0] == "ENDMDL"} {
            break
        }
    }
    puts "copy ${jobname}_run[expr ${iter} - 1]_0001.pdb to ${jobname}.pdb"
    close $input
    close $output
}

proc ::hprefstruct::count { temppdb } {
    set input [open $temppdb r]
    set temp1 [list]
    while { [gets $input line] >= 0 } {
        lappend temp1 [string range $line 21 21]
    }
    return [lsort -unique $temp1]
    close $input
}

proc ::hprefstruct::renumber { chain resid jobname } {
    set temp1 [open ${jobname}.pdb r]
    set id 0
    while { [gets $temp1 line] >= 0 } {
        if { $id == 0 } {
            incr id
            set temp2 [string range $line 21 21]
            set temp3 [string range $line 22 25]
        } elseif { [string range $line 21 21] != $temp2 || [string range $line 22 25] != $temp3 } {
            incr id
            set temp2 [string range $line 21 21]
            set temp3 [string range $line 22 25]
        } elseif { [string range $line 21 21] == $chain && [string range $line 22 25] == $resid } {
            return $id
            break
        }
    }
    close $temp1
    unset temp2 temp3
}

proc ::hprefstruct::find_resolution { pdb } {
    set input [open $pdb r]
    while {[gets $input line] >= 0} {
        if [regexp -- "REMARK   2 RESOLUTION\." $line] {
            close $input
            return [lindex $line 3]
            break
        }
    }
    close $input
}

proc ::hprefstruct::find_nbs2 { jobname cutoff reslist } {
    set id [mol new ${jobname}.pdb]
    set temp1   [list]
    foreach { chn res } $reslist {
        set tempsel1 [atomselect $id "chain $chn and resid $res"]
        lappend temp1 [$tempsel1 get index]
        $tempsel1 delete
    }
    set temp2   [join $temp1]
    set nbsel   [atomselect $id "same residue as within $cutoff of (index $temp2)"]
    set nbs     [$nbsel get index]
    $nbsel delete

    set temp3   [list]
    foreach index $nbs {
        set tmp [atomselect $id "index $index"]
        set ch [$tmp get chain]
        set re [$tmp get resid]
        if { [info exists $ch.$re] } {
        } else {
            set $ch.$re 0
            lappend temp3 $ch $re
        }
        $tmp delete
    }
    mol delete all
    return $temp3
}

proc ::hprefstruct::printH { chains } {
    variable iter
    set op [open run${iter}.xml w]
    puts $op "<ROSETTASCRIPTS>\n\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"dens\" weights=\"beta_cart\">\n\t\t\t<Reweight scoretype=\"elec_dens_fast\" weight=\"%%denswt%%\"/>\n\t\t\t<Set scale_sc_dens_byres=\"R:0.76,K:0.76,E:0.76,D:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88\"/>\n\t\t</ScoreFunction>\n\t</SCOREFXNS>\n\t<MOVERS>\n\t\t<SetupForDensityScoring name=\"setupdens\"/>\n\t\t<LoadDensityMap name=\"loaddens\" mapfile=\"%%map%%\"/>\n\t\t<FastRelax name=\"relaxcart\" scorefxn=\"dens\" repeats=\"3\" cartesian=\"1\" >\n\t\t\t<MoveMap>"
    set temp 1
    foreach chn $chains {
        puts $op "\t\t\t<Chain number=\"$temp\" chi=\"0\" bb=\"0\"/> Chain=$chn"
        incr temp
    }
    close $op
}

proc ::hprefstruct::printM { list bb jobname } {
    variable iter
    set op [open run${iter}.xml a]
	set chi 1
    foreach { chn res } $list {
        set resid [renumber $chn $res ${jobname}]
        puts $op "\t\t\t\t<Span begin=\"$resid\" end=\"$resid\" chi=\"${chi}\" bb=\"${bb}\"/> Chain=$chn Resid=$res"
    }
    puts $op "\t\t\t</MoveMap>\n\t\t</FastRelax>\n\t\tReportFSC name=\"report\" testmap=\"%%testmap%%\" res_low=\"10.0\" res_high=\"%%reso%%\"/>\n\t</MOVERS>\n\t<PROTOCOLS>\n\t\t<Add mover=\"setupdens\"/>\n\t\t<Add mover=\"loaddens\"/>\n\t\t<Add mover=\"relaxcart\"/>\n\t\tAdd mover=\"report\"\n\t</PROTOCOLS>\n\t<OUTPUT scorefxn=\"dens\"/>\n</ROSETTASCRIPTS>"
    close $op
    puts "Outputed run${iter} Rosetta xml file ..."
}

proc ::hprefstruct::printBash { pdb density resolution weight } {
    variable rosetta_bin
    variable rosetta_db
    variable iter
    set op [open rosetta.run${iter}.sh w]
    puts $op "#!/bin/bash\n$rosetta_bin \\\n\
        \t-database $rosetta_db \\\n\
        \t-in::file::s $pdb \\\n\
        \t-parser::protocol run${iter}.xml \\\n\
        \t-parser::script_vars denswt=${weight} reso=${resolution} map=${density} \\\n\
        \t-ignore_unrecognized_res \\\n\
        \t-edensity::mapreso $resolution \\\n\
        \t-edensity::cryoem_scatterers \\\n\
        \t\t-beta \\\n\
        \t-out::suffix _run${iter} \\\n\
        \t-crystal_refine
    "
    close $op
}

proc ::hprefstruct::op_results { jobname } {
    variable iter
    set output [open run${iter}.molprobity a]
    puts $output "\n\t===== Refine run${iter} parameters =====\t"
    if { $iter == 1 } {
        puts $output "  Initial pdb: ${jobname}_run0_0001.pdb"
        puts $output "  Refined pdb: ${jobname}_run1_0001.pdb"
    } else {
        puts $output "  Initial pdb: ${jobname}_run[expr ${iter} - 1]_0001.pdb"
        puts $output "  Refined pdb: ${jobname}_run${iter}_0001.pdb"
    }
    puts $output "  xml file: run${iter}.xml"
    puts $output "  bash file: rosetta.run${iter}.sh"
    close $output
}

proc ::hprefstruct::close_all_files { } {
    foreach channel [file channels "file*"] {
        close $channel
    }
}

#############################################################################
##     Execution
##############################################################################

proc ::hprefstruct::hprs { args } {

    if { [llength $args] % 2 } {
        puts "usage: hprs ?-arg var? ..."
        puts "\t-p pdb file"
        puts "\t-m map file (in .map or .mrc format)"
        puts "\t-r resolution of the map"
        puts "\t-j job name"
        puts "\t-n number of refinement iterations"
        puts "\t-w weight of the density "
        error "Error: odd number of arguments: $args"
    } elseif { [llength $args] == 0 } {
        puts "usage: hprs ?-arg var? ..."
        puts "\t-p pdb file"
        puts "\t-m map file (in .map or .mrc format)"
        puts "\t-r resolution of the map"
        puts "\t-j job name"
        puts "\t-n number of refinement iterations"
        puts "\t-w weight of the density "
        error "Error: empty argument ..."
    }

    foreach { name value } $args {
        switch -- $name {
            -p { set inputpdb $value }
            -m { set density $value }
            -r { set resolution $value }
            -j { set jobname $value }
            -n { set n $value }
            -w { set weight $value }
	    -s { set scheme $value }
            default { error "unknown argument: $name $value" }
        }
    }

    if { [file exists $inputpdb] == 0 } {
        error "Error: missing input pdb file ..."
    } elseif { [file exists $density] == 0 } {
        error "Error: missing input map file ..."
    } elseif { [info exists resolution] == 0 } {
        if { [::hprefstruct::find_resolution $inputpdb] == 0 } {
            error "Error: missing map resolution argument ...\nError: rerun refinement with -r arg? option to input the resolution"
        } else {
            set resolution [::hprefstruct::find_resolution $inputpdb]
        }
    } elseif { [info exists weight] == 0 } {
        set weight 5
    } elseif { [info exists jobname] == 0 } {
        set jobname temp
    }

    if { [info exists n] == 0 } {
	set start 1
    } else {
	set start $n
    }

	if { [info exists scheme] == 0 } {
		for { set i $start } { $i < 10 } { incr i } {
			variable iter $i
			display update ui
		        if { $i == 1 } {
		            file copy -force ${inputpdb} ${jobname}_run0_0001.pdb
		        }

		 	::hprefstruct::greppdb ${jobname}
			set problist [::hprefstruct::structural_analysis ${jobname}.pdb]
	        	::hprefstruct::printH [count ${jobname}.pdb] 
		        if { $i <= 3 } {
       			     ::hprefstruct::printM $problist 0 ${jobname} 
	      		} elseif { $i <= 6 } {
		            ::hprefstruct::printM $problist 1 ${jobname} 
	        	} elseif { $i <= 9 } {
		            ::hprefstruct::printM [find_nbs2 ${jobname} 2.0 $problist] 1 ${jobname} 
		        } else {
	        	    ::hprefstruct::printM [find_nbs2 ${jobname} 3.0 $problist] 1 ${jobname} 
		        }

		        ::hprefstruct::printBash ${jobname}.pdb $density $resolution $weight

	        	puts "Start run${i} refinement ..."
		        exec sh rosetta.run${i}.sh > run${i}.log
		        puts "Finished run${i} refinement ..."

	        	::hprefstruct::op_results ${jobname} 
		        puts "Output run${i} refinement parameters and results ..."
			::hprefstruct::close_all_files
    		}
	} else {
		variable iter $start
		display update ui
		if { $start == 1 } {
                	file copy -force ${inputpdb} ${jobname}_run0_0001.pdb
                }
		::hprefstruct::greppdb ${jobname}
                set problist [::hprefstruct::structural_analysis ${jobname}.pdb]
                ::hprefstruct::printH [count ${jobname}.pdb]
		
		set bb [lindex $scheme 0]
		set cutoff [lindex $scheme 1]
		
		::hprefstruct::printM [find_nbs2 ${jobname} $cutoff $problist] $bb ${jobname}
		
                puts "Start run${start} refinement ..."
                exec sh rosetta.run${start}.sh > run${start}.log
                puts "Finished run${start} refinement ..."

                ::hprefstruct::op_results ${jobname}
                puts "Output run${start} refinement parameters and results ..."
                ::hprefstruct::close_all_files
	}
}

