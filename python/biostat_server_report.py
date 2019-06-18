import paramiko
import getpass
import sys

un = input('username: ')
pw = getpass.getpass()
servers = ['carbon', 'cesium', 'chromium', 'potassium', 'silicon']

def server_report(server):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server + '.biostat.umn.edu', username = un, password = pw)
    stdin, stdout, stderr = ssh.exec_command("uptime")
    report = [server] + stdout.readlines()
    ssh.close()
    report = [item.rstrip() for item in report]	
    return report

print('Connecting to servers...')
try:
    reports = [server_report(server) for server in servers]
    for report in reports:
        report_str = '{0:12s} {1:s}'.format(report[0], report[1])
        print(report_str)
except:
    print('ssh failed. Check username/password and/or internet connection.')

