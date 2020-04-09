#############################################################################
# Protein structure refinement based on Cryo-EM maps and Molprobity scores
#############################################################################
#     Functions
#############################################################################

# molprobity score file parser, return geometery restraint outliers residue list
proc read_mpb { jobname iter } {
    puts "phenix.molprobity ${jobname}.pdb ..."
    set input [open "|phenix.molprobity ${jobname}.pdb" r]
    set output [open run${iter}.molprobity w]
    set temp1 [list]
    while {[gets $input oline] >= 0} {
        if { [regexp -- "-Bond lengths-" $oline] || [regexp -- "-Bond angles-" $oline] || [regexp -- "-Dihedral angles-" $oline] || [regexp -- "-Chiral volumes-" $oline] || [regexp -- "-Planar groups-" $oline] } {
            gets $input
            gets $input
            while { [gets $input iline] >= 0 } {
                if { [llength $iline] == 0 } {
                    break
                } elseif [string is integer -strict [string range $iline 4 7]] {
                    lappend temp1 "[string range $iline 3 3] [string range $iline 4 7]"
                } else {
                    error "Error when reading molprobity score output ..."
                }
            }
        } elseif [regexp -- "Summary" $oline] {
            gets $input
            puts $output "\t===== Molprobity score summary: before refine run${iter} =====\t"
            while { [gets $input iline] >= 0 } {
                if { [llength $iline] == 0 } {
                    break
                } else {
                    puts $output $iline
                }
            }
        }
    }
    flush $output
    close $input
    close $output
    unset -nocomplain bl ba da cv pg oline iline
    file delete molprobity_coot.py molprobity_probe.txt molprobity.out
    puts "Finished run${iter} molprobity analysis ..."
    return [join [lsort -u -dictionary $temp1]]
}

# molprobity score analysis
proc mpb { jobname iter } {
    #  variable iter
    #  variable jobname
    puts "phenix.molprobity ${jobname}_run${iter}_0001.pdb ..."
    set input [open "|phenix.molprobity ${jobname}_run${iter}_0001.pdb" r]
    set output [open run${iter}.molprobity a]
    while {[gets $input oline] >= 0} {
        if [regexp -- "Summary" $oline] {
            gets $input
            puts $output "\n\t===== Molprobity score summary: after refine run${iter} =====\t"
            while { [gets $input iline] >= 0 } {
                if { [llength $iline] == 0 } {
                    break
                } else {
                    puts $output $iline
                }
            }
        }
    }
    flush $output
    close $output
    close $input
    file delete molprobity_coot.py molprobity_probe.txt molprobity.out
    #  puts "Finished run${iter} molprobity analysis ..."
}

# return clash residue list
proc read_clash { jobname iter } {
    #  variable iter
    #  variable jobname
    if { $iter == 1 } {
        puts "phenix.clashscore model=${jobname}.pdb ..."
        set input [open "|phenix.clashscore model=${jobname}.pdb" r]
    } else {
        puts "phenix.clashscore model=${jobname}.pdb keep_hydrogens=True ..."
        set input [open "|phenix.clashscore model=${jobname}.pdb keep_hydrogens=True" r]
    }
    set temp1 [list]
    set pattern "Bad Clashes"
    while {[gets $input oline] >= 0} {
        if [regexp -- "Bad Clashes" $oline] {
            while { [gets $input iline] >= 0 } {
                if [regexp -- "clashscore" $iline] {
                    break
                } elseif { [string is integer -strict [string range $iline 2 5]] || [string is integer -strict [string range $iline 19 22]] } {
                    lappend temp1 "[string range $iline 1 1] [string range $iline 2 5]"
                    lappend temp1 "[string range $iline 18 18] [string range $iline 19 22]"
                } else {
                    error "Error when reading clash score output ..."
                }
            }
        }
    }
    close $input
    unset -nocomplain pattern oline iline
    puts "Finished run${iter} clashscore calculation ..."
    return [join [lsort -u -dictionary $temp1]]
}

# combine two residue lists
proc combine { list1 list2 } {
    set temp1 [list]
    foreach {chn id} $list1 {
        lappend temp1 "${chn} ${id}"
    }
    if { [llength $list2] != 0 } {
        foreach {chn id} $list2 {
            lappend temp1 "${chn} ${id}"
        }
    }
    unset -nocomplain  list1 list2
    #  return $temp1
    return [join [lsort -u -dictionary $temp1]]
}

# create a temperary pdb file
proc greppdb { jobname iter } {
    #  variable jobname
    set input [open ${jobname}_run[expr ${iter} - 1]_0001.pdb r]
    set output [open ${jobname}.pdb w]
    while {[gets $input line] >= 0} {
        if {[lindex $line 0] == "ATOM"} {
            puts $output $line
        } elseif {[lindex $line 0] == "ENDMDL"} {
            break
        }
    }
    flush $output
    close $input
    close $output
}

# count the number of chains
proc count { temppdb } {
    set input [open $temppdb r]
    set temp1 [list]
    while { [gets $input line] >= 0 } {
        lappend temp1 [string range $line 21 21]
    }
    return [lsort -unique $temp1]
    close $input
}

# renumbering the resid using the global resid with a temperary pdb file
proc renumber { chain resid jobname} {
    #  variable jobname
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
    flush $temp1
    close $temp1
    unset temp2 temp3
}

# find neighor residues using temp pdb file
proc find_nbs { chain resid } {
    variable jobname
    set cf 3.0
    set temp1 [mol new ${jobname}.pdb]
    
}

# find map resolution in pdb file
proc find_resolution { pdb } {
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

# output refinement results
proc op_results { jobname iter } {
    #  variable iter
    #  variable jobname
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
    flush $output
    close $output
}

# find the nerighor residues using VMD built-in atomselection
proc find_nbs2 { jobname cutoff reslist } {
#    variable jobname
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

## print Rosetta xml script
# print the head part
proc printH { chains iter } {
    #  variable iter
    set op [open run${iter}.xml w]
    puts $op "<ROSETTASCRIPTS>\n\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"dens\" weights=\"beta_cart\">\n\t\t\t<Reweight scoretype=\"elec_dens_fast\" weight=\"%%denswt%%\"/>\n\t\t\t<Set scale_sc_dens_byres=\"R:0.76,K:0.76,E:0.76,D:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88\"/>\n\t\t</ScoreFunction>\n\t</SCOREFXNS>\n\t<MOVERS>\n\t\t<SetupForDensityScoring name=\"setupdens\"/>\n\t\t<LoadDensityMap name=\"loaddens\" mapfile=\"%%map%%\"/>\n\t\t<FastRelax name=\"relaxcart\" scorefxn=\"dens\" repeats=\"3\" cartesian=\"1\" >\n\t\t\t<MoveMap>"
    set temp 1
    foreach chn $chains {
        puts $op "\t\t\t<Chain number=\"$temp\" chi=\"0\" bb=\"0\"/> Chain=$chn"
        incr temp
    }
    flush $op
    close $op
}

# print the main part
proc printM { list chi bb jobname iter } {
    #  variable iter
    set op [open run${iter}.xml a]
    foreach { chn res } $list {
        set resid [renumber $chn $res ${jobname}]
        puts $op "\t\t\t\t<Span begin=\"$resid\" end=\"$resid\" chi=\"${chi}\" bb=\"${bb}\"/> Chain=$chn Resid=$res"
    }
    puts $op "\t\t\t</MoveMap>\n\t\t</FastRelax>\n\t\tReportFSC name=\"report\" testmap=\"%%testmap%%\" res_low=\"10.0\" res_high=\"%%reso%%\"/>\n\t</MOVERS>\n\t<PROTOCOLS>\n\t\t<Add mover=\"setupdens\"/>\n\t\t<Add mover=\"loaddens\"/>\n\t\t<Add mover=\"relaxcart\"/>\n\t\tAdd mover=\"report\"\n\t</PROTOCOLS>\n\t<OUTPUT scorefxn=\"dens\"/>\n</ROSETTASCRIPTS>"
    flush $op
    close $op
    puts "Outputed run${iter} Rosetta xml file ..."
}

## printout Bash script
proc printBash { pdb density resolution weight iter } {
    #   variable iter
    #   variable density
    #   variable pdb
    variable binpath
    variable dbpath
    #   variable resolution
    #   variable weight
    set op [open rosetta.run${iter}.sh w]
    puts $op "#!/bin/bash\n$binpath \\\n\
        \t-database $dbpath \\\n\
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
    flush $op
    close $op
}

#############################################################################
#     Execution
#############################################################################

# input

set binpath  /misc/chaoyi/src_builds/rosetta/rosetta_src_2018.42.60459_bundle/main/source/bin/rosetta_scripts.linuxgccrelease
set dbpath   /misc/chaoyi/src_builds/rosetta/rosetta_src_2018.42.60459_bundle/main/database/

# main function to execute
proc proref { args } {
    variable rosetta_bin_path
    variable rosetta_db_path

    # print usage if wrong input parameters
    if { [llength $args] % 2 } {
        puts "usage: proref ?-arg var? ..."
        puts "\t-p pdb file"
        puts "\t-m map file (in .map or .mrc format)"
        puts "\t-r resolution of the map"
        puts "\t-j job name"
        puts "\t-n number of refinement iterations"
        puts "\t-w weight of the density "
        error "Error: odd number of arguments: $args"
    } elseif { [llength $args] == 0 } {
        puts "usage:tclsh proref ?-arg var? ..."
        puts "\t-p pdb file"
        puts "\t-m map file (in .map or .mrc format)"
        puts "\t-r resolution of the map"
        puts "\t-j job name"
        puts "\t-n number of refinement iterations"
        puts "\t-w weight of the density "
        error "Error: empty argument ..."
    }

    # input control
    foreach { name value } $args {
        switch -- $name {
            -p { set inputpdb $value }
            -m { set density $value }
            -r { set resolution $value }
            -j { set jobname $value }
            -n { set n $value }
            -w { set weight $value }
            default { error "unknown argument: $name $value" }
        }
    }

    # check if input files exist or arguments are defined
    if { [file exists $inputpdb] == 0 } {
        error "Error: missing input pdb file ..."
    } elseif { [file exists $density] == 0 } {
        error "Error: missing input map file ..."
    } elseif { [info exists resolution] == 0 } {
        if { [find_resolution $inputpdb] == 0 } {
            error "Error: missing map resolution argument ...\nError: rerun refinement with -r arg? option to input the resolution"
        } else {
            set resolution [find_resolution $inputpdb]
        }
    } elseif { [info exists weight] == 0 } {
        set weight 5
    } elseif { [info exists jobname] == 0 } {
        set jobname temp
    }

    # restart funtion is good
    if { [info exists n] == 0 } {
        set start 1
    } else {
        set start $n
    } 

    for { set iter $start } { $iter < 13 } { incr iter } {
        # initiate iter, select and modify the input pdb file name
        if { $iter == 1 } {
            file copy -force ${inputpdb} ${jobname}_run0_0001.pdb
        }

        # grep temp pdb file and set pdb for the initial refinement
        greppdb ${jobname} ${iter}

        # collect outlier lists
        set list1 [read_mpb ${jobname} ${iter}]
        set list2 [read_clash ${jobname} ${iter}]
#	set list2 [list]
        # print rosetta script and rosetta bash files
        printH [count ${jobname}.pdb] ${iter}
        if { $iter <= 3 } {
            printM [combine $list1 $list2] 1 0 ${jobname} ${iter}
        } elseif { $iter <= 6 } {
            printM [combine $list1 $list2] 1 1 ${jobname} ${iter}
        } elseif { $iter <= 9 } {
            printM [find_nbs2 ${jobname} 2.0 [combine $list1 $list2]] 1 1 ${jobname} ${iter}
	} else {
            printM [find_nbs2 ${jobname} 3.0 [combine $list1 $list2]] 1 1 ${jobname} ${iter}
	}

        printBash ${jobname}.pdb $density $resolution $weight $iter

        # execute Rosetta refinement
        puts "Start run${iter} refinement ..."
        exec sh rosetta.run${iter}.sh > run${iter}.log
        puts "Finished run${iter} refinement ..."

        # print refinement parameters and results
        op_results ${jobname} ${iter}
        mpb ${jobname} ${iter}
        puts "Output run${iter} refinement parameters and results ..."
    }
}


