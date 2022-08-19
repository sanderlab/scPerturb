import requests
import re
import time
import os
import subprocess
import pandas as pd
import yaml

with open("config.yaml", "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

DIR = config['DIR']
table = config['table']

# limit seems to be 3 (10 with a key) requests per second, see
# https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
def get_run_info(GSE):
    time.sleep(1)
    output = DIR + GSE+'/' + GSE + '_runinfo.csv'
    os.makedirs( DIR + GSE+'/', exist_ok=True)
    if not os.path.isfile(output):
        f = requests.get("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + GSE)
        m = re.search('PRJN?A[0-9]+', f.text)
        bioproject = m.group(0)
        print(f'Identified {GSE} as bioproject {bioproject}.')
        # "esearch -db sra -query '{bioproject} [PRJA]' | efetch -format runinfo > {output}"

        # This does not work, piping is empty?
        #ps = subprocess.Popen(("esearch", "-db", "sra", "-query", "'{"+bioproject+"} [PRJA]'"), stdout=subprocess.PIPE)
        #subprocess.Popen(("efetch", "-format", "runinfo"), stdin=ps.stdout, stdout=open(output, "wb"))

        # This is bad practice and poses a security risk (use of shell=True):
        subprocess.call("esearch -db sra -query '{"+bioproject+"} [PRJA]' | efetch -format runinfo", stdout=open(output, "wb"), shell=True)

if __name__ == "__main__":
    tab = pd.read_csv(table)
    GSEs = list(tab['GSE'][~pd.isna(tab['GSE'])])
    GSEs = [x for x in GSEs if x!='GSE116297']  # skip this is huge
    for GSE in GSEs:
        get_run_info(GSE)
