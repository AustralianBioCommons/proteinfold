#!/usr/bin/env python

import os
from matplotlib import pyplot as plt
import argparse
from collections import OrderedDict

def generate_output_images(msa_path, plddt_paths, name, out_dir):
    msa = []
    with open(msa_path, 'r') as in_file:
        for line in in_file:
            msa.append([int(x) for x in line.strip().split()])

    seqid = []
    for sequence in msa:
        matches = [1.0 if first == other else 0.0 for first, other in zip(msa[0], sequence)]
        seqid.append(sum(matches) / len(matches))

    seqid_sort = sorted(range(len(seqid)), key=seqid.__getitem__)

    non_gaps = []
    for sequence in msa:
        non_gaps.append([float(num != 21) if num != 21 else float('nan') for num in sequence])

    sorted_non_gaps = [non_gaps[i] for i in seqid_sort]
    final = []
    for sorted_seq, identity in zip(sorted_non_gaps, [seqid[i] for i in seqid_sort]):
        final.append([value * identity if not isinstance(value, str) else value for value in sorted_seq])



    ##################################################################
    plt.figure(figsize=(14, 14), dpi=100)
    ##################################################################
    plt.title("Sequence coverage")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    
    column_counts = [0] * len(msa[0])
    for col in range(len(msa[0])):
        for row in msa:
            if row[col] != 21:
                column_counts[col] += 1
                
    plt.plot(column_counts, color='black')
    plt.xlim(-0.5, len(msa[0]) - 0.5)
    plt.ylim(-0.5, len(msa) - 0.5)
    
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}seq_coverage.png")
    
    ##################################################################
    plddt_per_model = OrderedDict()
    plddt_paths_srt = plddt_paths
    plddt_paths_srt.sort()
    for plddt_path in plddt_paths_srt:
        with open(plddt_path, 'r') as in_file:
            plddt_per_model[os.path.basename(plddt_path)[:-4]] = [float(x) for x in in_file.read().strip().split()]

    plt.figure(figsize=(14, 14), dpi=100)
    plt.title("Predicted LDDT per position")
    for model_name, value_plddt in plddt_per_model.items():
        plt.plot(value_plddt, label=model_name)
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png")
    
    # split into figures
    i = 0
    for model_name, value_plddt in plddt_per_model.items():
        plt.figure(figsize=(14, 14), dpi=100)
        plt.title("Predicted LDDT per position")
        plt.plot(value_plddt, label=model_name)
        plt.ylim(0, 100)
        plt.ylabel("Predicted LDDT")
        plt.xlabel("Positions")
        plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT_{i}.png")
        i += 1
    
    
    ##################################################################

    
    ##################################################################
    """
    num_models = 5 # columns
    num_runs_per_model = math.ceil(len(model_names)/num_models)
    fig = plt.figure(figsize=(3 * num_models, 2 * num_runs_per_model), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(num_runs_per_model, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    fig.tight_layout()
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}PAE.png")
    """
    ##################################################################



print("Starting..")
parser = argparse.ArgumentParser()
parser.add_argument('--msa',dest='msa',required=True)
parser.add_argument('--plddt',dest='plddt',required=True, nargs="+")
parser.add_argument('--pdb',dest='pdb',required=True, nargs="+")
parser.add_argument('--name',dest='name')
parser.add_argument('--output_dir',dest='output_dir')
parser.add_argument('--html_template',dest='html_template')
parser.set_defaults(output_dir='')
parser.set_defaults(name='')
args = parser.parse_args()


generate_output_images(args.msa, args.plddt, args.name, args.output_dir)

print("generating html report...")
structures = args.pdb
structures.sort()

alphfold_template = open(args.html_template, "r").read()
alphfold_template = alphfold_template.replace(f"*sample_name_here*", args.name)
i = 0
for structure in structures:
    alphfold_template = alphfold_template.replace(f"*_data_ranked_{i}.pdb*", open(structure, "r").read().replace("\n", "\\n"))
    i += 1

with open(f"{args.output_dir}/{args.name}_alphafold.html", "w") as out_file:
    out_file.write(alphfold_template)
