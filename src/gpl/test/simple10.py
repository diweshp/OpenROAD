from openroad import Design, Tech
import helpers
import gpl_aux

tech = Tech()
tech.readLef("./sky130hd.lef")
design = Design(tech)
design.readDef("./simple10.def")

gpl_aux.global_placement(design, density=0.75)

def_file = helpers.make_result_file("simple10.def")
design.writeDef(def_file)
helpers.diff_files(def_file, "simple10.defok")
