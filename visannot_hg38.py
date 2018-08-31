import os
import logging
import multiprocessing

from tqdm import tqdm
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def make_box(ax, cx, cy, width, height, facecolor, edgecolor='white'):
    pat = patches.Rectangle(
        (cx, cy),               # low left point
        width, height,
        facecolor=facecolor,
        edgecolor=edgecolor,
    )
    ax.add_patch(pat)
    return pat


def calc_figsize(df_gene, trpt_height):
    """
    :param trpt_height: the box height for plotting a transcript
    """
    num_trpts = df_gene.transcript_id.unique().shape[0]
    fig_height = trpt_height * num_trpts
    return (30, fig_height)


def calc_coord_range(df_gene):
    beg = df_gene.start.min() - 1000
    end = df_gene.end.max() + 1000
    return beg, end


def calc_xlim(df_gene):
    beg, end = calc_coord_range(df_gene)
    return 0, end - beg


def calc_ylim(df_gene):
    return 0, df_gene.transcript_id.unique().shape[0]


def is_utr(row):
    return row.feature == "five_prime_utr" or row.feature == "three_prime_utr"


def calc_component_color(row):
    """
    :param row: a pandas dataframe row encoding the information of a particular
    component of the transcript
    """
    if is_utr(row):
        return dict(facecolor='green', edgecolor='white')
    else:
        dd = dict(
            protein_coding='blue',
            processed_transcript='magenta',
            retained_intron='orange',
            nonsense_mediated_decay='cyan',
        )
        return dict(facecolor=dd.get(row.transcript_biotype, 'black'),
                    edgecolor='white')


def calc_left_bottom_y_coordinate(row, nth_trpt, trpt_height):
    """
    :param nth_trpt: basically counting the index of the current transcript
    being plotted
    """
    if is_utr(row):
        return trpt_height * nth_trpt + trpt_height / 4
    else:
        return trpt_height * nth_trpt


def calc_component_height(row, trpt_height):
    """
    :param nth_trpt: basically counting the index of the current transcript
    being plotted
    """
    if is_utr(row):
        return trpt_height / 2
    else:
        return trpt_height


def plot_connector(ax, beg_x, end_x, coord_y):
    """
    :param currnet_row: the current transcript component being plotted
    currently
    """
    xs = [beg_x, end_x]
    ys = [coord_y, coord_y]
    ax.plot(xs, ys, lw=0.5, color='black')


def plot_transcript(ax, df_trpt, nth_trpt, trpt_height):
    previous_end_offset = None
    # three other possible values are CDS, start_codon and stop_codon.
    # exon contains CDS and UTRs
    # start_codon and stop_cond could be too small to see
    for fea in ['exon', "five_prime_utr", "three_prime_utr"]:
        sub_grp = df_trpt.query(f'feature == "{fea}"').sort_values('beg_offset')
        for key, row in sub_grp.iterrows():
            color_dd = calc_component_color(row)
            cy = calc_left_bottom_y_coordinate(row, nth_trpt, trpt_height)
            hght = calc_component_height(row, trpt_height)
            make_box(
                ax,
                cx=row.beg_offset,
                cy=cy,
                width=row.end_offset - row.beg_offset + 1,
                height=hght,
                **color_dd
            )

            # plot a line connecting transcript components
            if fea == 'exon':
                if previous_end_offset is not None:
                    y = cy + trpt_height / 2
                    plot_connector(ax, previous_end_offset, row.beg_offset, y)
                previous_end_offset = row.end_offset

    tid = df_trpt.transcript_id.unique()[0]
    y = trpt_height * nth_trpt + trpt_height / 2
    ax.text(0, y, tid, va='center', ha='right', fontsize=15)


def plot_transcripts(ax, df_gene, trpt_height):
    gene_info_tuple = get_gene_info_tuple(df_gene)

    # order by start position
    trpt_ids = df_gene.groupby('transcript_id').apply(lambda g: g.start.min()).sort_values().index.values

    dfs_trpts = df_gene.groupby('transcript_id')
    for k, tid in enumerate(trpt_ids):
        df_trpt = dfs_trpts.get_group(tid)
        logging.info(f'working on {gene_info_tuple}: {tid} ...')
        plot_transcript(ax, df_trpt, k, trpt_height)


def get_gene_info_tuple(df_gene):
    for i in ['seqname', 'gene_name', 'gene_id', 'gene_biotype', 'strand']:
        assert df_gene[i].unique().shape[0] == 1

    seqname = df_gene.seqname.unique()[0]
    gene_name = df_gene.gene_name.unique()[0]
    gene_id = df_gene.gene_id.unique()[0]
    gene_biotype = df_gene.gene_biotype.unique()[0]
    strand = df_gene.strand.unique()[0]
    return seqname, gene_name, gene_id, gene_biotype, strand


def decorate_ax(ax, df_gene):
    xlim = calc_xlim(df_gene)
    ylim = calc_ylim(df_gene)

    gene_info_tuple = get_gene_info_tuple(df_gene)
    seqname, gene_name, gene_id, gene_biotype, strand = gene_info_tuple
    ax.set_title(f'{seqname}:{xlim[0]}-{xlim[1]}, {gene_name}, {gene_id}, {gene_biotype}, ({strand})')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.grid()

    # so it doesn't overlap with the transcript ids
    ax.yaxis.tick_right()


def gen_output_png(root_dir, df_gene):
    seqname, gene_name, gene_id, _, _ = get_gene_info_tuple(df_gene)
    return os.path.join(
        root_dir, seqname, gene_name[:1], gene_name[:2], gene_name, f'{gene_id}.png')


def plot_gene(df_gene, trpt_height, outpng):
    """
    :param df_gene: the subset of gtf dataframe related to a gene
    """
    coord_range = calc_coord_range(df_gene)
    df_gene['beg_offset'] = df_gene['start'] - coord_range[0]
    df_gene['end_offset'] = df_gene['end'] - coord_range[0]

    figsize = calc_figsize(df_gene, trpt_height)
    logging.info(f'figsize: {figsize}')
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    plot_transcripts(ax, df_gene, trpt_height)
    decorate_ax(ax, df_gene)

    save_plot(fig, outpng)


def save_plot(fig, outpng):
    outdir = os.path.dirname(outpng)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    plt.savefig(outpng, bbox_inches='tight', dpi=250)
    logging.info(f'plotted {outpng}')
    fig.clf()
    plt.close(fig)


def load_gtf(pkl):
    adf = pd.read_pickle(pkl)

    # some cleaning
    for c in ['cds_end_NF', 'cds_start_NF', 'mRNA_end_NF', 'mRNA_start_NF',
              'CCDS', 'basic',
              'exon_number',
              'exon_version', 'protein_version', 'transcript_version']:
        adf[c] = adf[c].fillna(0).astype(int)

    for c in ['seqname']:
        adf[c] = adf[c].astype(str)

    assert adf.gene_id.unique().shape[0] == adf[['gene_id', 'seqname']].drop_duplicates().shape[0]
    assert adf.gene_id.unique().shape[0] == adf[['gene_id', 'gene_biotype']].drop_duplicates().shape[0]

    bdf = adf.query('feature != "gene"')\
             .query('feature != "transcript"')\
             .drop(['source', 'gene_source'], axis=1)
    return bdf


if __name__ == "__main__":
    pkl = '../gtf2csv/data/Homo_sapiens.GRCh38.92.pkl'
    root_outdir = './figs/hg38'
    trpt_height = 1
    num_cpus = 24

    if not os.path.exists(root_outdir):
        raise OSError(f"{root_outdir} doesn't exist yet!")

    df_gtf = load_gtf(pkl)

    grp_by_cols = ['seqname', 'gene_name', 'gene_id', 'gene_biotype', 'strand']
    logging.info(f'start grouping transcripts by {grp_by_cols} ...')
    args = []
    for key, grp in tqdm(df_gtf.groupby(grp_by_cols)):
        outpng = gen_output_png(root_outdir, grp)
        args.append((grp, trpt_height, outpng))

    with multiprocessing.Pool(num_cpus) as p:
        p.starmap(plot_gene, args)
