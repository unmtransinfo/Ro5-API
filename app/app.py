#main file. 1

from flask import Flask
from flasgger import Swagger
from ro5 import ro5_compute
from flask import jsonify
from flask import request

from flask_cors import CORS


import statistics
import math

import io, csv, base64
import os

RENDER_LIMIT = 1 #5000
#works but just shows in console.log. not in page.
REJECT_LIMIT = 100 


swagger_config = {
  "headers": [],
  "specs": [
    {
      "endpoint": "apispec_1",
      "route": "/apispec_1.json",
      "rule_filter": lambda rule: True, 
      "model_filter": lambda tag: True}
  ],
  "static_url_path": "/flasgger_static",
  "swagger_ui": True,
  "specs_route": "/apidocs/"
}



app = Flask(__name__)
Swagger(app, config=swagger_config)
CORS(app)

#test
@app.get("/htest")
def health():
  """Test check."""
  return {"status": "good"}


#helper func FOR box/histogram
def hist(values, bins=20):
  n = len(values)
  if n == 0:
    return {"bins": [], "counts": []}


  lo, hi = min(values), max(values)

  #if value same then just return
  if lo == hi:
    return {"bins": [round(lo,3), round(hi,3)], "counts": [n]}
  
  #step edges counts
  step = (hi - lo) / bins
  edges = [lo + i*step for i in range(bins+1)]
  counts = [0]*bins

  #which bin it belongs to
  for v in values:
    idx = min(int((v - lo) / step), bins-1)
    counts[idx] += 1


  return {
    "bins": [round(x,3) for x in edges],
    "counts": counts
  }
#helper func #2
def box(values):
  n = len(values)

  if n == 0:
    return {"min": None, "q1": None, "median": None, "q3": None, "max": None}

  if n == 1:
    x = round(values[0],3)
    return {"min": x, "q1": x, "median": x, "q3": x, "max": x}
    
  
  qs = statistics.quantiles(values, n=4, method="inclusive")

  return {
    "min": round(min(values),3),
    "q1": round(qs[0],3),
    "median": round(statistics.median(values),3),
    "q3": round(qs[2],3),
    "max": round(max(values),3),
  }


#for csv builder(double check later)
def build_csv(items):
  headers = [
        "smiles","mwt","logp","hbd","hba","violations","passes_ro5","vmax",
        "mwt_violation","hbd_violation","hba_violation","logp_violation","error"
    ]

  out = io.StringIO()
  w = csv.DictWriter(out, fieldnames=headers, extrasaction="ignore")
  w.writeheader()

  for it in items:
    row = {h: it.get(h, "") for h in headers}
    w.writerow(row)
  
  #get everything
  data = out.getvalue().encode("utf-8")


  return {
    "filename": "ro5_results.csv",
    "mime": "text/csv",
    "content": base64.b64encode(data).decode("ascii"),
}


#MAIN FUNCTION FOR COMPUTING SUMMARY
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
      "distributions": {
        "mwt":  {"hist": hist(vals("mwt")),"box": box(vals("mwt"))},
        "logp": {"hist": hist(vals("logp")),"box": box(vals("logp"))},
        "hbd":  {"hist":hist(vals("hbd")),"box": box(vals("hbd"))},
        "hba":  {"hist": hist(vals("hba")),"box": box(vals("hba"))},
      }
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
  n_input = len(smiles_list)
  if n_input >= REJECT_LIMIT:
    return jsonify({"error": f"Too many stuff: {n_input}. Right now max is {REJECT_LIMIT}"}), 413

  for s in smiles_list:
      r = ro5_compute(str(s).strip())
      if "error" not in r:
          r["vmax"] = vmax
          r["passes_ro5"] = r["violations"] <= vmax
      items.append(r)
    
  summary = compute_summary(items)

  res123 = {"items": items, "summary": summary}

  if n_input >= RENDER_LIMIT:
    res123["note"] = "Click the button below to download the summary file."
    #res123["items"] = []
    res123["download"] = build_csv(items)

  return jsonify(res123)


    



if __name__ == "__main__":
  port = int(os.getenv("APP_PORT", 8000))
  app.run(host="0.0.0.0", port=port)

