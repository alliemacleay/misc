import sys
import os
import time
from twilio.rest import TwilioRestClient
import subprocess

"""
Watch a job running with platform LSF and send a text when the job finishes

To get notifications about an LSF job
- - - - - - - - - - - - - - - - - - - - -
BSUB="bsub -u am282 -e lsf/e.log -o lsf/o.log -q medium"
$BSUB 'ls -a' | python /data/ngscid-clinical/home/ngscid/sandbox/am282/apps/misc/track.py &
"""


def wait_for_finish(job, frequency, timeout):
    """
    Sleeps for a defined interval then checks whether a job id is still running
    :param job: job_id eg. 415883
    :param frequency: seconds to wait between queries
    :param timeout: maximum seconds to keep looping
    :return: True|False whether the job completed
    """
    job_done = False
    start = time.time()
    duration = 0
    while not job_done and duration < timeout:
        time.sleep(frequency)
        cmd = ["bjobs", job]
        response = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.readlines()
        duration = time.time() - start
        job_done, status = is_job_done(response)
    print "{} - {}".format(job_num, status)
    return job_done, status


def is_job_done(res):
    """
    Parses a bjobs response to whether or not the job is done
    :param res: platform LSF response array
    :return: True if the job is done, False if it is not
    """
    status = ''
    job_done = False
    if 'is not found\n' in res[0]:  # response: "Job <452341> is not found\n"
        job_done = True
    else:
        if len(res) == 2:
            # response: JOBID   USER    STAT  QUEUE      FROM_HOST   EXEC_HOST   JOB_NAME   SUBMIT_TIME
            #           536475  ngscid  DONE  medium     eris1n2     cn012       ls -a      Apr 18 14:55
            status = res[1].split()[2]
            if status in ['DONE', 'EXIT']:
                job_done = True
    return job_done, status


def twilio_alert(twilio_account, twilio_pw, msg):
    """ Sent Twilio Text message """
    client = TwilioRestClient(twilio_account, twilio_pw)
    client.messages.create(to=to_number, from_=from_number, body=msg)

if __name__ == "__main__":
    frequency = 60  # number of seconds to sleep before checking for the job again
    timeout = 60 * 60 * 10  # ten hours
    twilio_account = os.getenv("TWILIO_ACT")
    twilio_pw = os.getenv("TWILIO_PW")
    from_number = os.getenv("TWILIO_NUM")
    to_number = os.getenv("MY_PHONE")
    if '' in [twilio_account, twilio_pw, from_number, to_number]:
        raise ValueError("Make sure TWILIO_ACT, TWILIO_PW, TWILIO_NUM, and MY_PHONE are set in the environment")

    if not sys.stdin.isatty():
        job_num = sys.stdin.read().split('<')[1].split('>')[0]
    else:
        job_num = sys.argv[1]
    print "Waiting for {} to finish.  When done a text will be sent to {}".format(job_num, to_number)
    done, status = wait_for_finish(job_num, frequency, timeout)
    if not done:
        message = "Job {} has timed out after {} seconds".format(job_num, timeout)
    else:
        message = "Job {} has finished with status of {}".format(job_num, status)
    print message
    twilio_alert(twilio_account, twilio_pw, message)

