# Run: uvicorn main:app

from fastapi import FastAPI, HTTPException
import jsonpickle
from pydantic import BaseModel
from vina import Vina

app = FastAPI()

class Request(BaseModel):
    prot_pdbqt_file: str
    lig_pdbqt_str: str
    pocket_center_json_pickle: str
    box_size_json_pickle: str

RECEPTORS_DIR = "data/test_set"

@app.post("/")
async def dock(data: Request):
    
    v = Vina(sf_name="vina", seed=0, verbosity=0)
    v.set_receptor(f"{RECEPTORS_DIR}/{data.prot_pdbqt_file}")
    v.set_ligand_from_string(data.lig_pdbqt_str)

    pocket_center = jsonpickle.decode(data.pocket_center_json_pickle)
    box_size = jsonpickle.decode(data.box_size_json_pickle)
    v.compute_vina_maps(center=pocket_center, box_size=box_size)
    score = v.score()[0]

    return { "score": score }