#!/usr/bin/env python
import subprocess
import sys
import logging

LOG = logging.getLogger("__name__")

COMPLETE_STATUS = "COMPLETED"
RUNNING_STATUSES = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]

COMPLETED_SIGNAL = "success"
RUNNING_SIGNAL = "running"
FAILED_SIGNAL = "failed"

jobid = sys.argv[1]

job_status = FAILED_SIGNAL

try:
    sacct_output = subprocess.check_output(["sacct", "-P", "-b", "-j", str(jobid), "-n"])
    for sacct_line in sacct_output.decode().strip().split("\n"):
        if sacct_line.split("|")[0] == str(jobid):
            if COMPLETE_STATUS in sacct_line:
                job_status = COMPLETED_SIGNAL
                break
            elif any(status in sacct_line for status in RUNNING_STATUSES):
                job_status = RUNNING_SIGNAL
                break
            else:
                job_status = FAILED_SIGNAL
                break

except Exception as e:
    LOG.error("failed to run sacct")
    LOG.error(e)
    job_status = FAILED_SIGNAL

print(job_status)
