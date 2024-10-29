# exclude horizontal edges
from openroad import Design, Tech
import helpers
import ppl_aux

tech = Tech()
tech.readLef("Nangate45/Nangate45.lef")
design = helpers.make_design(tech)
design.readDef("gcd.def")

ppl_aux.place_pins(
    design,
    hor_layers="metal3",
    ver_layers="metal2",
    corner_avoidance=0,
    min_distance=0.12,
    group_pins=[
        "clk req_msg[31] req_msg[30] req_msg[29] req_msg[28] req_msg[27] req_msg[26] req_msg[25] req_msg[24] req_msg[23] req_msg[22] req_msg[21] req_msg[20]"
    ],
)

def_file = helpers.make_result_file("group_pins2.def")
design.writeDef(def_file)
helpers.diff_files("group_pins2.defok", def_file)


# # exclude horizontal edges
# source "helpers.tcl"
# read_lef Nangate45/Nangate45.lef
# read_def gcd.def

# place_pins -hor_layers metal3 -ver_layers metal2 -corner_avoidance 0 -min_distance 0.12 -group_pins {clk req_msg[31] req_msg[30] req_msg[29] req_msg[28] req_msg[27] req_msg[26] req_msg[25] req_msg[24] req_msg[23] req_msg[22] req_msg[21] req_msg[20]}

# set def_file [make_result_file group_pins2.def]

# write_def $def_file

# diff_file group_pins2.defok $def_file
