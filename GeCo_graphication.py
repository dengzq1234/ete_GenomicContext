import argparse
import random
import math

import pandas as pd
from json import loads


from ete3 import Tree, TreeStyle, add_face_to_node, Face,\
                    FaceContainer, CircleFace, TextFace
from faces import ArrowFace
from qt import (QBrush, QPen, QColor, QPixmap, QPainter, Qt, QRect,
                QPainterPath)


from get_context import launch_analysis


def get_operons_from_json(inputfile):
    with open(inputfile, "r") as handler:
        ops = loads(handler.read())
    return eval(ops)


def get_unique_notation(operons, notation, level):
    unique = {}
    for central_gene in operons.values():
        for pos, gene in central_gene['neighborhood'].items():
            if (abs(int(pos)) <= nside):
                unique = { **unique,
                           **get_notation(gene, notation, level)
                         }
    return unique


def get_notation(gene, notation, level):
    unique = {}
    if notation == "eggNOG" and gene[notation] != {}:
        items = gene[notation][level]
    else:
        items = gene[notation]
    for f, d in items.items():
        if f != "" and f != "scores":
            try:
                unique[f] = d['description']
            except:
                unique[f] = ""
    return unique


def get_palette(colors_file, unique_notation):
    with open(colors_file, "r") as handle:
        colors = eval(handle.read())

    while len(colors) < len(unique_notation):
        colors.extend(colors)
    palette = {
        unique_notation[i] : colors[i] for i in
             range(len(unique_notation))
    }
    return palette


def arrow_layout(node):
    if node.is_leaf():
        operon = operons[node.name]['neighborhood']
        for pos, gene in operon.items():
            pos = int(pos)
            if (abs(pos)) <= nside:
                unigene = str(gene['unigene'])
                if unigene != "nan":
                    colors = []
                    for n in get_notation(gene, notation, level).keys():
                        if n != "":
                            colors.append(palette[n])
                    if len(colors) == 0:
                        colors = ["#cfcfcf"]
                    try:
                        strand = gene['strand']
                    except:
                        strand = "+"
                    geneFace = ArrowFace(30, 20, strand, colors)
                    add_face_to_node(geneFace, node,
                                     column=pos+nside, position="aligned")


def style_tree(ts, unique_notation, palette=False):
    ts.show_branch_support = True
    ts.branch_vertical_margin = 10
    ts.show_scale = True

    ts.layout_fn = arrow_layout
    if palette:
        ts.legend.add_face(TextFace("\t\t\t\t", fsize=12), column=0)
        ts.legend.add_face(TextFace(notation + " legend", fsize=12), column=1)
        ts.legend.add_face(CircleFace(5, "#cfcfcf"), column=0)
        ts.legend.add_face(TextFace("No data"), column=1)
        for name, color in palette.items():
            ts.legend.add_face(CircleFace(5, color), column=0)
            description = unique_notation[name].strip()
            # description = description.split(" ")
            # for i in range(7, len(description), 7):
                # description[i] += "\n"
            # description = " ".join(description).strip()
            ts.legend.add_face(TextFace(name + ": " + description),
                               column=1)
        ts.legend_position = 4
    return ts


def arg_parser():

    parser = argparse.ArgumentParser(description='Genomic context \
                                    visualization using ETE.')
    parser.add_argument('--cluster', type=str, required=True,
                        help="Cluster to \
                         visualize: (XXX_XXX_XXX)")
    parser.add_argument("--operons", type=str, help="File containing \
                            operon data in JSON format")
    parser.add_argument("--output", type=str, help="Output file")
    parser.add_argument('--tree', type=str, help="Filepath to Newick file")
    parser.add_argument('--nside', type=int, help="Number of neighbor genes \
                            up/downstream of central gene. Default: 10. Max: 20")
    parser.add_argument('--notation', type=str, required = True,
                        help="Functional notation. E.g. KEGG, eggNOG...")
    parser.add_argument('--level', type=str, help="Level to specify in \
                                    certain notation. E.g. eggNOG, level 2")
    args = parser.parse_args()
    return args


def main():

    args = arg_parser()
    cluster = args.cluster
    #cluster = "001_754_949"
    if args.operons:
        operon_file = args.operons
    else:
        operon_file = "data/" + cluster + ".txt"

    if args.tree:
        tree_file = args.tree
    else:
        tree_file = "data/" + cluster + "_newick.txt"

    if args.output:
        output_file = args.output
    else:
        output_file = "results/" + cluster + "_GeCo"

    colors_file = "colors.txt"

    global nside, notation, level, operons, palette
    if args.nside:
        nside = args.nside
    else:
        nside = 10
    notation = args.notation
    level = args.level

    # chunk of neigh gene vis
    launch_analysis(cluster,
                          nside,
                          30,
                          True)
    operons = get_operons_from_json(operon_file)
    unique_notation = get_unique_notation(operons, notation, level)
    palette = get_palette(colors_file, list(unique_notation.keys()))

    t = Tree(tree_file)
    ts = TreeStyle()
    ts = style_tree(ts, unique_notation, palette)

    t.render(output_file + "_1.png", dpi=1200, tree_style=ts)


main()


# E.g.
# ipython GeCo_graphication.py -- --cluster 001_756_949 --notation KEGG
