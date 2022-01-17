#!/usr/bin/from

from __future__ import print_function
import requests
import sys
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Symbol to Ensembl gene')
parser.add_argument('-g', '--gene', default="ENSG00000139618", dest='gene', help='gene symbol')
args = parser.parse_args()


server = "https://grch37.rest.ensembl.org"
ext = "/lookup/id/"
ext += args.gene
ext += "?"

r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
descript = "None"
if 'description' in decoded.keys():
  descript = decoded['description']
print(decoded['id'], decoded['display_name'], descript, decoded['strand'], decoded['assembly_name'], sep="\t")
