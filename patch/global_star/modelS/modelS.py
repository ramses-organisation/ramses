
def ReadModelS():
  f = open('data/cptrho.l5bi.d.15c', 'r')
  lines = f.readlines()
  f.close()

  lines = lines[6:]
  r = []
  c = []
  rho = []
  p = []
  gamma_1 = []
  T = []

  for i in lines:
      columns = i.split()
      r.append(float(columns[0]))
      c.append(float(columns[1]))
      rho.append(float(columns[2]))
      p.append(float(columns[3]))
      gamma_1.append(float(columns[4]))
      T.append(float(columns[5]))

  return [r,c,rho,p,gamma_1,T]
