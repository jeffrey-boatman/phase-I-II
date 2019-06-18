import paramiko
import getpass
import sys


un = input('username: ')
pw = getpass.getpass()
servers = ['carbon', 'cesium', 'chromium', 'potassium', 'silicon']

# def get_cmd(server):
#     ssh = paramiko.SSHClient()
#     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     ssh.connect(server + '.biostat.umn.edu', username = un, password = pw)
#     stdin, stdout, stderr = ssh.exec_command("ps -u" + un)
#     jobs = stdout.readlines()
#     jobs = [line.split(' ') for line in jobs]
#     ssh.close()
#     pids2kill = [job[0] for job in jobs if job[-1] == 'R\n']
#     cmds = ['kill ' + pid + ';' for pid in pids2kill]
#     kill_command = ''
#     for cmd in cmds:
#       kill_command = kill_command + cmd
#     ssh.close
#     return(kill_command)


# def kill_pids(server):
#     ssh = paramiko.SSHClient()
#     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#     ssh.connect(server + '.biostat.umn.edu', username = un, password = pw)
#     stdin, stdout, stderr = ssh.exec_command("ps -u" + un)
#     jobs = stdout.readlines()
#     jobs = [line.split(' ') for line in jobs]
#     pids2kill = [job[1] for job in jobs if job[-1] == 'R\n']
#     cmds = ['kill ' + pid + ';' for pid in pids2kill]
#     kill_command = ''
#     for cmd in cmds:
#       kill_command = kill_command + cmd
#     stdin, stdout, stderr = ssh.exec_command(kill_command)
#     ssh.close()
def kill_pids(server):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server + '.biostat.umn.edu', username = un, password = pw)
    stdin, stdout, stderr = ssh.exec_command("pkill -u" + un)
#     jobs = stdout.readlines()
#     jobs = [line.split(' ') for line in jobs]
#     pids2kill = [job[1] for job in jobs if job[-1] == 'R\n']
#     cmds = ['kill ' + pid + ';' for pid in pids2kill]
#     kill_command = ''
#     for cmd in cmds:
#       kill_command = kill_command + cmd
#     stdin, stdout, stderr = ssh.exec_command(kill_command)
    ssh.close()


try:
  for server in servers:
    kill_pids(server)
    print("jobs killed on %s." %server)
except:
  print('Could not kill jobs. Check username/password, code directory and/or internet connection.')

# jobs = get_pids('carbon')
# jobs = [line.split(' ') for line in jobs]
# pids2kill = [job[0] for job in jobs if job[-1] == 'R']
# cmds = ['kill ' + pid + ';' for pid in pids2kill]
# big_command = ''
# for cmd in cmds:
#   big_command = big_command + cmd




# print('Connecting to servers...')
# try:
#     reports = [server_report(server) for server in servers]
#     for report in reports:
#         report_str = '{0:12s} {1:s}'.format(report[0], report[1])
#         print(report_str)
# except:
#     print('ssh failed. Check username/password and/or internet connection.')

