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

proc ::hprefstruct::count { temppdb } {
    set input [open $temppdb r]
    set temp1 [list]
    while { [gets $input line] >= 0 } {
        lappend temp1 [string range $line 21 21]
    }
    return [lsort -unique $temp1]
    close $input
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
        lappend temp1 [[atomselect $id "chain $chn and resid $res"] get index]
    }
    set temp2   [join $temp1]
    set nbs     [[atomselect $id "same residue as within $cutoff of (index $temp2)"] get index]

    set temp3   [list]
    foreach index $nbs {
        set ch [[atomselect $id "index $index"] get chain]
        set re [[atomselect $id "index $index"] get resid]
        if { [info exists $ch.$re] } {
        } else {
            set $ch.$re 0
            lappend temp3 $ch $re
        }
    }
    mol delete all
    return $temp3
}

proc ::hprefstruct::readline2 { line } {
        set l1 [lindex [split [lindex [split $line (] 1] )] 0]
        set l2 [lindex [split [lindex [split $line (] 2] )] 0]
        set l3 [lindex [split [lindex [split [lindex [split $line (] 2] )] 1] !] 0]
        set id1 [::hprefstruct::renumber A [lindex $l1 1] refine]
        set id2 [::hprefstruct::renumber A [lindex $l2 1] refine]
        # AtomPair Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum FLAT_HARMONIC x0 sd tol
        return "AtomPair [lindex $l1 4] $id1 [lindex $l2 4] $id2 FLAT_HARMONIC [format %.2f [expr [lindex $l3 0] + ([lindex $l3 2] - [lindex $l3 1])/2]] 0.02 [format %.2f [expr abs(([lindex $l3 2] + [lindex $l3 1])/2)]]"
        unset l1 l2 l3
}

proc ::hprefstruct::get_tors { pdb } {
        set id [mol new $pdb]
        set chn [lsort -u [[atomselect $id protein] get chain]]
        set op [open "torsion.cst" w]
        puts $chn
        foreach i $chn {
                set resid [lsort -u -int [[atomselect $id "protein and chain $i"] get resid]]
                foreach j $resid {
                        puts "resid $j"
                        # phi
                        if {$j != [lindex $resid 0]} {
                                set phi [[atomselect $id "protein and chain $i and resid $j and alpha"] get phi]
                                # Dihedral Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum Atom3_Name Atom3_ResNum Atom4_Name Atom4_ResNum Func_Type Func_Def
                                puts $op "Dihedral C [::hprefstruct::renumber A [expr $j - 1] refine] N [::hprefstruct::renumber A $j refine] CA [::hprefstruct::renumber A $j refine] C [::hprefstruct::renumber A $j refine] CIRCULARHARMONIC [expr $phi*3.1415926/180] 0.02"
                        }
                        # psi
                        if {$j != [lindex $resid end]} {
                                set psi [[atomselect $id "protein and chain $i and resid $j and alpha"] get psi]
                                # Dihedral Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum Atom3_Name Atom3_ResNum Atom4_Name Atom4_ResNum Func_Type Func_Def
                                puts $op "Dihedral N [::hprefstruct::renumber A $j refine] CA [::hprefstruct::renumber A $j refine] C [::hprefstruct::renumber A $j refine] N [::hprefstruct::renumber A [expr $j + 1] refine] CIRCULARHARMONIC [expr $psi*3.1415926/180] 0.02"
                        }
                }
        }
        close $op
        mol delete $id
}

proc ::hprefstruct::readtbl { tbl } {
        set ip [open $tbl]
        set data [read $ip]
        close $ip
        set op [open [lrange [split $tbl .] 0 end-1].cst w]
        foreach line [split $data \n] {
                if { [lindex $line 0] == "assign" } {
                        puts $op "[::hprefstruct::readline2 $line]"
                }
        }
        unset data
        close $op
}
                 

proc ::hprefstruct::printH { chains } {
    variable iter
    set op [open run${iter}.xml w]
    puts $op "<ROSETTASCRIPTS>\n\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"r15_cart\" weights=\"ref2015\">\n\t\t\t<Reweight scoretype=\"pro_close\" weight=\"0.0\"/>\n\t\t\t<Reweight scoretype=\"cart_bonded\" weight=\"0.625\"/>\n\t\t\t<Reweight scoretype=\"atom_pair_constraint\" weight=\"%%weight%%\"/>\n\t\t\t<Reweight scoretype=\"dihedral_constraint\" weight=\"%%weight%%\"/>\n\t\t</ScoreFunction>\n\t</SCOREFXNS>\n\t<MOVERS>\n\t\t<FastRelax cartesian=\"1\" name=\"relaxcart\" repeats=\"5\" scorefxn=\"r15_cart\">\n\t\t\t<MoveMap>"
    set temp 1
    foreach chn $chains {
        puts $op "\t\t\t<Chain number=\"$temp\" chi=\"0\" bb=\"0\"/> Chain=$chn"
        incr temp
    }
    close $op
}

proc ::hprefstruct::printxml2 { } {
	variable iter
	set op [open run${iter}.xml w]
	puts $op "<ROSETTASCRIPTS>\n\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"r15_cart\" weights=\"ref2015\">\n\t\t\t<Reweight scoretype=\"pro_close\" weight=\"0.0\"/>\n\t\t\t<Reweight scoretype=\"cart_bonded\" weight=\"0.625\"/>\n\t\t\t<Reweight scoretype=\"atom_pair_constraint\" weight=\"%%weight%%\"/>\n\t\t\t<Reweight scoretype=\"dihedral_constraint\" weight=\"%%weight%%\"/>\n\t\t</ScoreFunction>\n\t</SCOREFXNS>\n\t<MOVERS>\n\t\t<FastRelax cartesian=\"1\" name=\"relaxcart\" repeats=\"5\" scorefxn=\"r15_cart\">\n\t\t\t<MoveMap>\n\t\t\t\t<Chain bb=\"1\" chi=\"1\" number=\"1\"/>\n\t\t\t</MoveMap>\n\t\t</FastRelax>\n\t\t<AddConstraints name=\"add_NMR\">\n\t\t\t<FileConstraintGenerator filename=\"%%nmr%%\" name=\"NMR_cst\"/>\n\t\t</AddConstraints>\n\t\t<AddConstraints name=\"add_tor\">\n\t\t\t<FileConstraintGenerator filename=\"%%talos%%\" name=\"TOR_cst\"/>\n\t\t</AddConstraints>\n\t\t<RemoveConstraints constraint_generators=\"NMR_cst\" name=\"rm_NMR\"/>\n\t</MOVERS>\n\t<PROTOCOLS>\n\t\t<Add mover=\"add_NMR\"/>\n\t\t<Add mover=\"add_tor\"/>\n\t\t<Add mover=\"relaxcart\"/>\n\t</PROTOCOLS>\n\t<OUTPUT scorefxn=\"r15_cart\"/>\n</ROSETTASCRIPTS>"
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
    puts $op "\t\t\t</MoveMap>\n\t\t</FastRelax>\n\t\t<AddConstraints name=\"add_NMR\">\n\t\t\t<FileConstraintGenerator filename=\"%%nmr%%\" name=\"NMR_cst\"/>\n\t\t</AddConstraints>\n\t\t<AddConstraints name=\"add_tor\">\n\t\t\t<FileConstraintGenerator filename=\"%%talos%%\" name=\"TOR_cst\"/>\n\t\t</AddConstraints>\n\t\t<RemoveConstraints constraint_generators=\"NMR_cst\" name=\"rm_NMR\"/>\n\t</MOVERS>\n\t<PROTOCOLS>\n\t\t<Add mover=\"add_NMR\"/>\n\t\t<Add mover=\"add_tor\"/>\n\t\t<Add mover=\"relaxcart\"/>\n\t</PROTOCOLS>\n\t<OUTPUT scorefxn=\"r15_cart\"/>\n</ROSETTASCRIPTS>"
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

proc ::hprefstruct::printBash2 { pdb k NMR TALOS } {
    variable rosetta_bin
    variable rosetta_db
    variable iter
    set op [open rosetta.run${iter}.sh w]
    puts $op "#!/bin/bash\n$rosetta_bin \\\n\
        \t-database $rosetta_db \\\n\
        \t-in::file::s $pdb \\\n\
        \t-parser::protocol run${iter}.xml \\\n\
        \t-parser::script_vars weight=$k nmr=${NMR} talos=${TALOS} \\\n\
        \t-ignore_unrecognized_res \\\n\
        \t-out::suffix _run${iter} \\\n\
        \t-crystal_refine
    "
    close $op
}

###############################
#
###############################

proc ::hprefstruct::hprs { args } {
    for {set i 0} {$i < [llength $args]} {incr i} {
            switch -- [lindex $args $i] {
                -em { set mode 1 }
                -EM { set mode 1 }
                -nmr { set mode 2 }
                -NMR { set mode 2 }
                -p {
                        set inputpdb [lindex $args $i+1]
                        incr i
                }
                -m {
                        set density [lindex $args $i+1]
                        incr i
                }
                    -r {
                            set resolution [lindex $args $i+1]
                        incr i
                }
                -j { set jobname [lindex $args $i+1]
                        incr i
                }
                    -n { set n [lindex $args $i+1]
                        incr i
                }
                    -w { set weight [lindex $args $i+1]
                        incr i
                }
                    -s { set scheme [lindex $args $i+1]
                        incr i
                }
			-cst-dist { set dist [lindex $args $i+1] 
				incr i
			} 
			-cst-talos { set talos [lindex $args $i+1]                                             
				incr i
			}

                default { error "Unknown input argument: [lindex $args $i]" }
            }
    }

    if { [file exists $inputpdb] == 0 } {
        error "Error: missing input pdb file ..."
    } elseif { [info exists jobname] == 0 } {
        set jobname temp
    }


    if { $mode == 1 } {
            if { [file exists $density] == 0 } {
                error "Error: missing input map file ..."
            } elseif { [info exists resolution] == 0 } {
                if { [::hprefstruct::find_resolution $inputpdb] == 0 } {
                    error "Error: missing map resolution argument ...\nError: rerun refinement with -r arg? option to input the resolution"
                } else {
                    set resolution [::hprefstruct::find_resolution $inputpdb]
                }
            } elseif { [info exists weight] == 0 } {
                set weight 5
            }
    }

    if { [info exists n] == 0 } {
        set start 1
    } else {
        set start $n
    }

    if { $mode == 2} {
	    
	for {set i $start} {$i <= 5} {incr i} {
		variable iter $i

	    if { $i == 1 } {
		    file copy -force ${inputpdb} ${jobname}_run0_0001.pdb        
	    
	    ::hprefstruct::greppdb ${jobname} 

	    ::hprefstruct::printxml2
	    ::hprefstruct::printBash2 ${jobname}.pdb $weight $dist $talos
	    
	    puts "Start run${i} refinement ..."
	    exec sh rosetta.run${i}.sh > run${i}.log
	    puts "Finished run${i} refinement ..."
	} else {
		
	
	    ::hprefstruct::greppdb ${jobname} 
			set problist [::hprefstruct::structural_analysis ${jobname}.pdb]
		::hprefstruct::printH [count ${jobname}.pdb] 
		if {$i == 2} {
           	 ::hprefstruct::printM $problist 0 ${jobname} 
		} elseif {$i == 3} {
           	 ::hprefstruct::printM $problist 1 ${jobname} 
		} elseif {$i == 4} {
		            ::hprefstruct::printM [find_nbs2 ${jobname} 2.0 $problist] 1 ${jobname} 
		} elseif {$i == 5} {
		            ::hprefstruct::printM [find_nbs2 ${jobname} 3.0 $problist] 1 ${jobname} 
		}
		        ::hprefstruct::printBash2 ${jobname}.pdb $weight $dist $talos
	    puts "Start run${i} refinement ..."
	    exec sh rosetta.run${i}.sh > run${i}.log
	    puts "Finished run${i} refinement ..."
    }
	}
    }
}
