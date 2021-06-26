fp = open('run', 'w')
for i in range(801):
  fp.write(f'python animation.py {i}\n')
fp.close()
