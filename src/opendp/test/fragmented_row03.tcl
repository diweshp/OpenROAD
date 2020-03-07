source helpers.tcl
read_lef Nangate45.lef
read_def fragmented_row03.def
detailed_placement
catch {check_placement} error
puts $error

set def_file [make_result_file fragmented_row03.def]
write_def $def_file
diff_file $def_file fragmented_row03.defok
