import os
import json
from vina import Vina
import jsonpickle

def lambda_handler(event, context):

    try:

        if "body" in event:
            # API Gateway
            data = json.loads(event["body"])
        else:
            # Lambda CLI invoke
            data = event

        prot_pdbqt_file = data["prot_pdbqt_file"]
        lig_pdbqt_str = data["lig_pdbqt_str"]
        pocket_center = jsonpickle.decode(data["pocket_center_json_pickle"])
        box_size = jsonpickle.decode(data["box_size_json_pickle"])

        v = Vina(sf_name="vina", seed=0, verbosity=0)
        v.set_receptor(f"./data/{prot_pdbqt_file}")
        v.set_ligand_from_string(lig_pdbqt_str)

        v.compute_vina_maps(center=pocket_center, box_size=box_size)
        
        score = v.score()[0]
        
        return {
            "statusCode": 200,
            "body": json.dumps({"score": score})
        }

    except Exception as e:
        return {
            "statusCode": 400,
            "body": f"{str(type(e))}: {str(e)}"
        }