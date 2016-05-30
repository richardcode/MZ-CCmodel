import numpy as np

def load_carbondata(d_file):

  f_in = open(d_file,'r')
  f_line = f_in.read()
  data = f_line.split('\n')[:-1]
  data = [f.split('  ') for f in data]
  for i in range(0,len(data)):
    for j in range(0,len(data[0])):
      data[i][j] = float(data[i][j].replace('D','E'))
  data = np.array(data)
  return data

def load_thermaldata(d_file,nlevs=101,ntime=1001,layered=True):

  f_in = open(d_file,'r')
  f_line = f_in.read()
  data = f_line.split('\n')[:-1]
  data = [float(f) for f in data]
  data = np.array(data)
  if layered== True:
    data = data.reshape(nlevs,ntime)
  return data