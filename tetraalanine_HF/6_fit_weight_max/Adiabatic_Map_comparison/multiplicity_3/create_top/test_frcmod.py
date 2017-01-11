from parmed.amber import * 
import parmed

base = AmberParm("fit.prmtop","fit.rst7")
parmed.tools.writeFrcmod(base,"test.frcmod").execute()
