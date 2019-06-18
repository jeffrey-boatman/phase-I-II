import paramiko
import getpass
import sys

un = input('username: ')
pw = getpass.getpass()
nice = input('nice: ')

nice_options = list(range(1, 20))
nice_options = [str(num) for num in nice_options]
if str(nice) not in nice_options:
  nice = '19'

servers = ['carbon', 'cesium', 'chromium' , 'potassium', 'silicon']
# servers = ['carbon', 'cesium', 'potassium' , 'silicon']

command = "nohup nice +" + nice + " R-3.5 CMD BATCH --vanilla apply-mems.R &"

def runsim(server):
  ssh = paramiko.SSHClient()
  ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
  ssh.connect(server + '.biostat.umn.edu', username = un, password = pw)
  ssh.exec_command("cd sim/phase-I-II; " + command)
  ssh.close()

try:
  for server in servers:
    runsim(server)
    print("jobs submitted to %s." %server)
except:
  print('ssh failed. Check username/password, code directory and/or internet connection.')
