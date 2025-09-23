#main file. 1

from flask import Flask
from flasgger import Swagger
from ro5 import ro5_compute
from flask import jsonify
from flask import request

app = Flask(__name__)
Swagger(app)


#test
@app.get("/htest")
def health():
    """Test check."""
    return {"status": "good"}

#ethanol: CCO
#curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CCO","vmax":1}' | jq

#aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
#curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","vmax":1}' | jq
@app.post("/ro5")
def ro5_res():
    """
    Lipinski rule of Five
    ---
    parameters:
      - in: body
        name: body
        schema:
          type: object
          properties:
            smiles:
              type: array
              items: {type: string}
        vmax: {type: integer, default: 1}
    responses:
      200:
        description: Ro5 results
    """
    data = request.get_json(force=True, silent=True) or {}


    smiles_list = data.get("smiles", [])
    vmax = int(data.get("vmax", 1))
    items = []

    #if single
    if isinstance(smiles_list, str):
        smiles_list = [smiles_list]


    for s in smiles_list:
        r = ro5_compute(str(s).strip())
        if "error" not in r:
            r["vmax"] = vmax
            r["passes_ro5"] = r["violations"] <= vmax
        items.append(r)
    return jsonify({"items": items})


    



if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000)

