import numpy as numpy


def make_opc_to_twox():

  year =31557600
  dt=0.25e-2*year
  tmax=1.1e4*year
  t1=0.0
  t2=t1+140.*year
  n1=int(t1/dt)
  n2=int(t2/dt)
  index = np.arange(n1,n2+1)
  t = dt * index
  pco2 = 28.35 * np.exp(np.log(1.01)*t/year)

  pco2[pco2>2.0*28.35] = 2.0*28.35

  fout = open('pco2_opc_2x.txt','w')
  for i in range(0,len(pco2)):
    fout.write(str(pco2[i])+'\n')
  fout.close()

  return

def make_opc_to_fourx():

  year =31557600
  dt=0.25e-2*year
  tmax=1.1e4*year
  t1=0.0
  t2=t1+140.*year
  n1=int(t1/dt)
  n2=int(t2/dt)
  index = np.arange(n1,n2+1)
  t = dt * index
  pco2 = 28.35 * np.exp(np.log(1.01)*t/year)

  pco2[pco2>4.0*28.35] = 4.0*28.35

  fout = open('pco2_opc_4x.txt','w')
  for i in range(0,len(pco2)):
    fout.write(str(pco2[i])+'\n')
  fout.close()

  return
