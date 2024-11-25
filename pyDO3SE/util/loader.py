import json
import csv

def json_loader(fp):
    with open(fp) as f:
        return json.load(f)

def csv_loader(fp):
    with open(fp) as f:
        spamreader = csv.reader(f)
        header = next(spamreader)
        out = []
        for row in spamreader:
            out.append({k:v for k, v in zip(header, row)})
        return out
