import sys
import os
import time
from twilio.rest import TwilioRestClient
import smtplib
from email.mime.text import MIMEText
import subprocess

"""
This works with a Twilio account or to an ATT phone number

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


def email(to_number, subject, message, cc=None):
    """
    http://lifehacker.com/5506326/how-can-i-send-an-email-via-text-message
    Alltel: phonenumber@message.alltel.com
    AT&T: phonenumber@txt.att.net
    T-Mobile: phonenumber@tmomail.net
    Virgin Mobile: phonenumber@vmobl.com
    Sprint: phonenumber@messaging.sprintpcs.com
    Verizon: phonenumber@vtext.com
    Nextel: phonenumber@messaging.nextel.com
    US Cellular: phonenumber@mms.uscc.net
    """
    server, port = 'phsmgout.partners.org', 25
    text_type = 'html'
    to_number = to_number.replace('+', '')
    if to_number[0] == '1':
        to_number = to_number[1:]
    send_to = '{}@txt.att.net'.format(to_number)
    sendfrom = 'CID@MGH.HARVARD.EDU'
    msg = MIMEText(message, text_type)
    msg['From'] = sendfrom
    msg['To'] = ','.join(send_to)
    msg['Subject'] = subject
    smtp = smtplib.SMTP(timeout=timeout)
    if cc is not None:
        send_to.extend(cc)
        msg['Cc'] = ','.join(cc)
    try:
        smtp.connect(server, port)
    except Exception as e:
        print('Email must be sent from behind the Partners firewall:\nError:{}'.format(e.message))
        raise e
    smtp.sendmail(sendfrom, send_to, msg.as_string())
    smtp.quit()


if __name__ == "__main__":
    frequency = 60  # number of seconds to sleep before checking for the job again
    timeout = 60 * 60 * 10  # ten hours
    twilio_account = os.getenv("TWILIO_ACT", "")
    twilio_pw = os.getenv("TWILIO_PW", "")
    from_number = os.getenv("TWILIO_NUM", "")
    to_number = os.getenv("MY_PHONE", "")
    if '' in [twilio_account, twilio_pw, from_number, to_number]:
        raise ValueError("Make sure MY_PHONE is set in the environment")

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
    if '' not in [twilio_account, twilio_pw, from_number]:
        twilio_alert(twilio_account, twilio_pw, message)
    else:
        email(to_number, 'Job status', message)

