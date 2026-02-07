import statistics

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