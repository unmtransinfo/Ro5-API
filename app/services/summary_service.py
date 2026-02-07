import statistics
from utils.stats import hist, box

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