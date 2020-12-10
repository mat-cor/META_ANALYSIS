#!/usr/bin/env python3
import shlex
from subprocess import Popen, PIPE,call,check_output
import argparse,datetime,subprocess


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build Docker file and push to gcr")

    parser.add_argument("--image", type= str,
                        help="name of image",default = 'saige')
    parser.add_argument("--version", type= str,
                        help="version value, e.g.0.001",required = True)
    parser.add_argument("--push",action = 'store_true')
    parser.add_argument("--args",type = str,default = '')
    args = parser.parse_args()

    
    basic_cmd = 'docker build -t eu.gcr.io/covid19-hgi-ganna/' + args.image +':' +args.version
    cmd = basic_cmd + ' -f Dockerfile_SAIGE_GWAS ..' + ' ' + args.args
    print(cmd)
    call(shlex.split(cmd))

    if args.push:
        cmd = 'gcloud docker -- push eu.gcr.io/covid19-hgi-ganna/' + args.image +':' + args.version
        call(shlex.split(cmd))