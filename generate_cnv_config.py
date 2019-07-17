import os
import argparse
import json

"""
Generates a config file of the following format:
{
    "mappings": "/global/scratch2/sweram19/mapping/mapping",
    "input_samples": {
        "SRS396822_Donald_Pan-troglodytes-verus-x-troglodytes_male": [
            "SRR747947.bam",
            "SRR747948.bam",
            "SRR747949.bam",
            "SRR747950.bam"
            ]
        "SRS394730_Vincent_Pan-troglodytes-schweinfurthii_male": [
            "SRR726241.bam",
            "SRR726242.bam",
            "SRR726243.bam"
            ]
        }
}

To be used with the following Snakefile:
### insert Snakefile path here ###
"""

__author__ = "Swetha Ramesh"


def make_json(mappings_path):
    j_out = {"mappings": mappings_path}
    samples = {}
    cwd = os.getcwd()
    dirs = os.listdir(mappings_path)
    for d in dirs:
        files = []
        for f in os.listdir(os.path.join(mappings_path, d)):
            files.append(f)
        samples[d] = files
    j_out['input_samples'] = samples
    FOUT = open("{cwd}/config.json".format(cwd=cwd),'w')
    FOUT.write(json.dumps(j_out, indent=4, separators=(",", ": ")))
    FOUT.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mappings_path", required=True)
    o = parser.parse_args()

    make_json(o.mappings_path)
