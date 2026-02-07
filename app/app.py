#main file. 1

from flask import Flask
from flasgger import Swagger
from ro5 import ro5_compute
from flask import jsonify
from flask import request

from flask_cors import CORS

from services.summary_service import compute_summary
from utils.csv_utils import build_csv

import os

RENDER_LIMIT = 1 #5000
REJECT_LIMIT = 1500 

app = Flask(__name__)

# Load
URL_PREFIX = os.getenv("URL_PREFIX", "").strip("/")
API_BASE = f"/{URL_PREFIX}/api" if URL_PREFIX else "/api"

IN_PROD = app.config.get("FLASK_ENV", "") == "production"
print("IN_PROD value:", IN_PROD)

prefix = os.getenv("SWAGGER_PREFIX", "").rstrip("/")
# Configure swagger config 
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
  "specs_route": "/apidocs/",
  "ui_params": {"url_prefix": prefix},
}




Swagger(app, config=swagger_config)
CORS(app)

#test
@app.get(f"{API_BASE}/htest")
def health():
  """Test check."""
  return {"status": "good"}





#ethanol: CCO
#curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CCO","vmax":1}' | jq

#aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
#curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","vmax":1}' | jq
@app.post(f"{API_BASE}/ro5")
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
  names_list = data.get("names", [])

  vmax = int(data.get("vmax", 1))
  items = []

  #if single
  if isinstance(smiles_list, str):
    smiles_list = [smiles_list]
  if isinstance(names_list, str):
    names_list = [names_list]
  
  n_input = len(smiles_list)
  if n_input >= REJECT_LIMIT:
    return jsonify({"error": f"Too many compounds provided: {n_input}. Current max is {REJECT_LIMIT}"}), 413

  for idx, s in enumerate(smiles_list):
      r = ro5_compute(str(s).strip())
      if "error" not in r:
          r["vmax"] = vmax
          r["passes_ro5"] = r["violations"] <= vmax
      ##attach name if provided
      if idx < len(names_list) and names_list[idx]:
        r["name"] = names_list[idx].strip()

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

