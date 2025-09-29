#main file. 1

from flask import Flask
from flasgger import Swagger
from ro5 import ro5_compute
from flask import jsonify
from flask import request

from flask_cors import CORS


import statistics

app = Flask(__name__)
Swagger(app)
CORS(app)

#test
@app.get("/htest")
def health():
    """Test check."""
    return {"status": "good"}


def compute_summary(items):
  valid = [i for i in items if "error" not in i]
  n = len(valid)

  def vals(key):
    return [i[key] for i in valid]

  def stats(v):
    #if 0
    if len(v) == 0:
       return {"n": 0, "mean": None, "stdev": None}
    #if 1
    if len(v) == 1:
       return {"n": 1, "mean": round(v[0], 3), "stdev": 0.0}
    #ekse
    return {
        "n": len(v),
         "mean": round(statistics.mean(v), 3),
        "stdev": round(statistics.stdev(v), 3),
    }


  #pass and fail stuff
  pass_count = sum(1 for i in valid if i.get("passes_ro5"))
  fail_count = n - pass_count
  def pct(k): 
    if n:
      return round(100.0 * k / n, 2)
    else:
      return None


  #final return
  return {
      "mwt":  stats(vals("mwt")),
      "logp": stats(vals("logp")),
      "hbd":  stats(vals("hbd")),
      "hba":  stats(vals("hba")),
      "pass_fail": {
            "n": n,
            "pass": {"count": pass_count, "pct": pct(pass_count)},
            "fail": {"count": fail_count, "pct": pct(fail_count)},
        },
  }

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
    
    summary = compute_summary(items)
    return jsonify({"items": items, "summary": summary})


    



if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000)

