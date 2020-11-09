import argparse
import random
import math

import pandas as pd
from json import loads


from ete3 import Tree, TreeStyle, add_face_to_node, Face,\
                    FaceContainer, CircleFace, TextFace, AttrFace, RectFace, SeqGroup, SeqMotifFace
from ete3.treeview import random_color

from faces import ArrowFace
import colorsys

from qt import (QBrush, QPen, QColor, QPixmap, QPainter, Qt, QRect,
                QPainterPath)


from get_context import (launch_analysis, 
                            Cluster2TreeAlg, 
                            Unigene2Cluster, 
                            Cluster2Members, 
                            Unigene2Taxa, 
                            Unigene2Biome, 
                            Unigene2Trembl, 
                            Unigene2Emapper)


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

# Pick two colors. Values from 0 to 1. See "hue" at
# http://en.wikipedia.org/wiki/HSL_and_HSV
def color_gradient(hue, intensity, granularity):
    min_lightness = 0.35 
    max_lightness = 0.9
    base_value = intensity

    # each gradient must contain 100 lightly descendant colors
    colors = []   
    rgb2hex = lambda rgb: '#%02x%02x%02x' % rgb
    l_factor = (max_lightness-min_lightness) / float(granularity)
    l = min_lightness
    while l <= max_lightness:
        l += l_factor
        rgb =  rgb2hex(tuple(map(lambda x: int(x*255), 
                                 colorsys.hls_to_rgb(hue, l, base_value))))
        colors.append(rgb)
        
    colors.append("#ffffff")
    return colors

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


# set layout for biome distribution and taxa
SPACER = 15 # for 15 biomes
def layout(node):
     # set space for biome
    distF = AttrFace("dist", fsize=8, formatter="%0.1g", fgcolor="#555555")
    rednameF = AttrFace("name", fsize=8, fgcolor="indianred")
    bluenameF = AttrFace("name", fsize=8, fgcolor="steelblue")
    redgradient = color_gradient(0.95, 0.6, 10)
    node.img_style["size"] = 0

    if not node.up:
        return

    if node.name.startswith("GMGC"):
        #print(node.name)
        name = node.name.split('.')[1]
    else:
        name = node.name

    if node.is_leaf():
        colidx = 0

        # for gradient
        if gene2trembl[name]:
            ident_trembl = gene2trembl[name][1]
            color_idx = int(ident_trembl*10)
            color = redgradient[color_idx]


            identF = TextFace("%.1f" % (ident_trembl*100), fgcolor="black")
            identF.background.color = color
            add_face_to_node(identF, node, column=colidx, position="aligned")
            identF.margin_left = SPACER
            colidx += 1

            hitF = TextFace(gene2trembl[name][0])        
            add_face_to_node(hitF, node, column=colidx, position="aligned")
            hitF.margin_left= SPACER
            colidx += 1
        else:
            ident_trembl = 0
            color_idx = int(ident_trembl)
            color = redgradient[color_idx]


            identF = TextFace(ident_trembl, fgcolor="black")
            identF.background.color = color
            add_face_to_node(identF, node, column=colidx, position="aligned")
            identF.margin_left = SPACER
            colidx += 1

            hitF = TextFace(" ") #gene2trembl[name]
            add_face_to_node(hitF, node, column=colidx, position="aligned")
            hitF.margin_left= SPACER
            colidx += 1

        # emapper annotation
        
        if gene2emapper[name]:
            emapperF = TextFace(gene2emapper[name][0])
            add_face_to_node(emapperF, node, column=colidx, position="aligned")
            emapperF.margin_left= SPACER
            colidx +=1

            if gene2emapper[name][1] != "NA|NA|NA":
                ogsF = TextFace(gene2emapper[name][1])
            else:
                ogsF = TextFace(gene2emapper[name][2])
            add_face_to_node(ogsF, node, column=colidx, position="aligned")
            ogsF.margin_left= SPACER
            colidx +=1

        else:
            emapperF = TextFace(" ")
            add_face_to_node(emapperF, node, column=colidx, position="aligned")
            emapperF.margin_left= SPACER
            colidx +=1

            ogsF = TextFace(" ")
            add_face_to_node(ogsF, node, column=colidx, position="aligned")
            ogsF.margin_left= SPACER
            colidx +=1
        

        # set space for taxa
        
        if gene2taxa[name]:
            taxaF = TextFace(gene2taxa[name][2])
            add_face_to_node(taxaF, node, column=colidx, position="aligned")
            taxaF.margin_left= SPACER
            colidx +=1
        else:
            taxaF = TextFace(" ")
            add_face_to_node(taxaF, node, column=colidx, position="aligned")
            taxaF.margin_left= SPACER
            colidx +=1
        

        # set color of biome distribution
        biome_keys= [
                #'am',
                #'iso',
                'bu',
                'was',
                'soil',
                'mar',
                'fw',
                'pig',
                'cat',
                'dog',
                'mous',
                'gut',
                'or',
                'nos',
                'skin',
                'vag'
            ]
        colors = random_color(num=len(biome_keys), l=0.5, s=0.35, h=60)
        biome2color = {biome:colors[i] for i, biome in enumerate(biome_keys)}
        biome2face = {}
        for b in biome_keys: 
            biome2face[b] = RectFace(5, 10, fgcolor=biome2color[b], bgcolor=biome2color[b])
        biome2face["white"] = RectFace(5, 10, fgcolor="white", bgcolor="white")

        for i, b in enumerate(biome_keys):
            if gene2biomes[name] and (gene2biomes[name][i] > 0):
                #print(gene2biomes[name])
                paF = TextFace(gene2biomes[name][i], fsize=9, ftype="Courier", fgcolor = "white" )
                paF.background.color = biome2color[b]
            else: 
                paF = TextFace("  ", fsize=9, ftype="Courier", fgcolor = "white" )
                paF.background.color = "#D3D3D3"

            paF.margin_left = 1
            paF.margin_right = 2
            #paF.border.color = "black"

            add_face_to_node(paF, node, colidx, position="aligned")        
            colidx +=1 
        
        # not integrate algnment yet
        # seqF = SeqMotifFace(alg.get_seq(node.name), seq_format="compactseq", scale_factor=0.2)      
        # seqF.margin_left = SPACER
        # add_face_to_node(seqF, node, column=colidx, position="aligned")

        # for neigh
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
                    # add_face_to_node(geneFace, node,
                    #                  column=pos+nside, position="aligned")
                    add_face_to_node(geneFace, node,
                                     column=colidx+pos+nside, position="aligned")
                
        if 0 < ident_trembl < 0.95:
            color_name = "SteelBlue"
        elif ident_trembl == 0:
            color_name = "#8B0000"
        else: 
            color_name = "#333333"
            name = gene2trembl[name][0]

        nameF =  TextFace(name, fgcolor=color_name)
        add_face_to_node(nameF, node, column=0, position="branch-right")

    add_face_to_node(distF, node, column=0, position="branch-top")
    return

def style_tree(ts, unique_notation, palette=False):
    ts.show_branch_support = True
    ts.show_leaf_name = False
    ts.branch_vertical_margin = 10
    ts.show_scale = True
    
    #ts.layout_fn = arrow_layout
    ts.layout_fn = layout
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
    parser.add_argument('--aln', type=str, help="Filepath to Alignment file")
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
        
    if args.aln:
        aln_file = args.aln

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

    # analysis of neigh genes and visualization
    launch_analysis(cluster,
                          nside,
                          30,
                          True)
    operons = get_operons_from_json(operon_file)
    unique_notation = get_unique_notation(operons, notation, level)
    palette = get_palette(colors_file, list(unique_notation.keys()))

    # get cluster members
    gmgcfam_clm = Cluster2Members(cluster)

    # get taxa, biomes, trembl, eggnog annotation
    global gene2taxa, gene2biomes, gene2trembl, gene2emapper
    gene2taxa = {}
    gene2biomes = {}
    gene2trembl = {}
    gene2emapper = {}
    for unigene in gmgcfam_clm:
        gene2taxa[unigene] = Unigene2Taxa(unigene)
        gene2biomes[unigene] = Unigene2Biome(unigene)
        gene2trembl[unigene] = Unigene2Trembl(unigene)
        gene2emapper[unigene] = Unigene2Emapper(unigene)

    # get cluster family tree and alignment
    gmgcfam_nw, gmgcfam_aln = Cluster2TreeAlg(cluster)

    # build tree with set style
    t = Tree(tree_file)
    t.standardize()
    ts = TreeStyle()

    ts = style_tree(ts, unique_notation, palette)

    try:
        t.set_outgroup(t.get_midpoint_outgroup())
    except:
        pass

    t.render(output_file + ".png", dpi=2000, tree_style=ts)


main()


# E.g.
# ipython GeCo_graphication.py -- --cluster 001_756_949 --notation KEGG
