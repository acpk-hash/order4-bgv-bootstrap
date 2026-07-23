#!/usr/bin/env python3
"""Convert butterfly_stages.json to a flat text file for the C++ prototype."""
import json
from pathlib import Path
OUT = Path("/home/user/experiments/order4bgv/tower_linear_transform/encrypted_prototype")
exp = json.load(open(OUT/"butterfly_stages.json"))
D = exp["D"]; d = exp["d"]
with open(OUT/"stages.txt","w") as f:
    f.write(f"{exp['p']} {D} {d}\n")
    f.write(" ".join(map(str,exp["G"]))+"\n")
    allst = exp["ntt_stages"]+exp["intt_stages"]
    f.write(f"{len(allst)}\n")
    for st in allst:
        offs = st["offsets"]
        f.write(f"{len(offs)}\n")
        for o in offs:
            f.write(f"{o}\n")
            dv = st["diags"][str(o)]
            for i in range(D):
                f.write(" ".join(map(str,dv[i]))+"\n")
print("wrote stages.txt")
