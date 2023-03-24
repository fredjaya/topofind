import subprocess

def run_command(rid, cmd):
    """
    General function to run command-line processes.
    """
    # Popen lets you access the I/O pipes.
    # stdout and stderr options specify which pipes you want to capture.
    process = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # communicate() waits for the process to complete and returns a tuple of
    # the captured pipes.
    stdout, stderr = process.communicate()
    exit_code = process.returncode
    return stdout, stderr, exit_code
