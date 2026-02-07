import io, csv, base64
#for csv builder
def build_csv(items):
  headers = [
        "smiles",
        "name",
        "mwt","logp","hbd","hba","violations","passes_ro5","vmax",
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