"""
Protein class for handling protein data and performing various operations.
This is a work-in-progress.

Classes
-------
Protein: Class for handling protein data and performing various operations

Functions
---------
generate_pair_sequence: Generate pair sequence from two proteins, segmenting the sequence into regions with low or high pLDDT
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from Bio import SeqIO
from Bio.PDB import PDBParser
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter1d

from .._settings import get_setting
from ._math import get_3d_avg, normalize, smooth, square_grad
from .data import (
    get_lddt_from_uniprot_id,
    get_seq_from_uniprot_id,
    get_uniprot_from_gene_name,
)

pio.templates.default = "plotly_white"

bioparser = PDBParser()


def generate_pair_sequence(P1, P2, output_dir):
    """generate pair sequence from row, segmenting the sequence into regions with low or high pLDDT

    Parameters
    ----------
    P1: str
        First protein
    P2: str
        Second protein
    output_dir: str
        Output directory
    """
    protein_a = Protein(P1)
    protein_b = Protein(P2)
    low_or_high_plddt_region_sequence_a = protein_a.low_or_high_plddt_region_sequence
    low_or_high_plddt_region_sequence_b = protein_b.low_or_high_plddt_region_sequence
    for i, seq_a in enumerate(low_or_high_plddt_region_sequence_a):
        for j, seq_b in enumerate(low_or_high_plddt_region_sequence_b):
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            filename = (
                Path(output_dir)
                / f"{protein_a.gene_name}_{i}_{protein_b.gene_name}_{j}.fasta"
            )
            with Path(filename).open("w") as f:
                f.write(
                    f">{protein_a.gene_name}_{i}.{protein_b.gene_name}_{j}\n{str(seq_a.seq)}:{str(seq_b.seq)}\n"
                )


class Protein:
    """Protein class for handling protein data and performing various operations"""

    def __init__(
        self,
        gene_name,
        homodimer=False,
        use_es=False,
        esm_folder="/manitou/pmg/users/xf2217/demo_data/esm1b/esm1b_t33_650M_UR90S_1",
        af2_folder="/manitou/pmg/users/xf2217/demo_data/af2",
        xml_dir=Path(get_setting("cache_dir")) / "uniprot_xml",
        window_size=10,
    ):
        """
        Parameters
        ----------
        gene_name (str): gene name
        homodimer (bool): whether to use homodimer structure
        use_es (bool): whether to use ES score
        esm_folder (str): path to ESM1b model
        af2_folder (str): path to AF2 model
        xml_dir (str): path to uniprot XML directory
        window_size (int): size of mutation window
        """
        self.gene_name = gene_name
        self.window_size = window_size
        self.esm_folder = esm_folder
        self.af2_folder = af2_folder
        self.xml_dir = xml_dir
        # create if not exists
        self.xml_dir.mkdir(parents=True, exist_ok=True)
        self.uniprot_id = get_uniprot_from_gene_name(gene_name)
        self.plddt = get_lddt_from_uniprot_id(self.uniprot_id)
        self.length = len(self.plddt)
        self.sequence = get_seq_from_uniprot_id(self.uniprot_id)
        self.smoothed_plddt_gaussian = self.get_smoothed_plddt_gaussian()
        self.smoothed_plddt = self.get_smooth_plddt()
        self.domains = self.get_domain_from_uniprot()

        if use_es:
            self.grad = self.get_plddt_grad()
            self.pairwise_distance = self.get_pairwise_distance(
                dimer=False, af2_folder=self.af2_folder
            )
            self.esm = self.get_esm()
            self.es = smooth(
                self.get_final_score_gated_grad_3d(), window_size=self.window_size
            )

    def get_domain_from_uniprot(self):
        """Get domain information from uniprot"""
        from .uniprot import UniProtAPI

        uniprot = UniProtAPI()
        return uniprot.get_domains(self.uniprot_id, xml_dir=self.xml_dir)

    def get_smooth_plddt(self, window_size=10):
        result = np.convolve(
            self.plddt, np.ones(window_size) / window_size, mode="same"
        )
        return result / np.max(result)

    def get_smoothed_plddt_gaussian(self, sigma=2):
        result = gaussian_filter1d(self.plddt, sigma=sigma)
        return result / np.max(result)

    def get_plddt_grad(self):
        grad = square_grad(self.smoothed_plddt_gaussian)
        grad = np.clip(grad, np.quantile(grad, 0.2), np.quantile(grad, 0.8))
        return grad

    def get_esm(
        self,
        format="pos",
    ):
        """Use UNIPROT_ID to get ESM score from ESM1b model"""
        # read PAX5 (PAX5) | Q02548.csv
        # variant,score,pos
        # M1K,-8.983,1
        # M1R,-8.712,1
        df = pd.read_csv(
            f"{self.esm_folder}/{self.uniprot_id}_LLR.csv",
            index_col=0,
        )
        if format == "long":
            # melt to long format, column, row, value
            df = df.reset_index().melt(id_vars="index")
            df["variant"] = df["variable"].str.replace(" ", "") + df["index"].astype(
                str
            )
            df["pos"] = df["variable"].apply(lambda x: int(x.split(" ")[1]))
            df = df.rename({"value": "esm"}, axis=1)
            df = df[["variant", "pos", "esm"]]
            return df
        elif format == "wide":
            return df
        elif format == "pos":
            return normalize(-df.mean(0).values)
        # df['ALT'] = df.variant.apply(extract_alt_from_mut)
        # df['REF'] = df.variant.apply(extract_wt_from_mut)
        # df['esm'] = df['esm'].astype(float)
        # df['pos'] = df['pos'].astype(int)
        # df = df[['esm', 'ALT', 'pos']].pivot_table(index='ALT', columns='pos', values='esm').fillna(0)
        # df = normalize(-df.mean(0).values)
        # return df

    @property
    def low_plddt_region(self, threshold=0.6):
        idx = np.where(self.smoothed_plddt < threshold)[0]
        # get regions from idx, join two regions if they are close (<30aa apart)
        regions = []
        for i in idx:
            if len(regions) == 0:
                regions.append([i, i + 1])
            else:
                if i - regions[-1][1] < 30:
                    regions[-1][1] = i
                else:
                    regions.append([i, i + 1])

        for i in regions:
            if i[1] - i[0] < 30:
                regions.remove(i)
        return regions

    @property
    def low_plddt_region_sequence(self, threshold=0.6):
        regions = self.low_plddt_region
        sequences = []
        for i, region in enumerate(regions):
            s = SeqIO.SeqRecord(
                self.sequence[region[0] : region[1]],
                id=self.gene_name + "_" + str(i),
                name=self.gene_name,
                description=self.gene_name
                + "_"
                + str(i)
                + ": "
                + str(region[0])
                + "-"
                + str(region[1]),
            )
            sequences.append(s)
        return sequences

    @property
    def low_or_high_plddt_region(self, threshold=0.6):
        """
        Get regions with low plddt and high plddt, define breakpoint by merge regions if they are close (<30bp apart)

        Parameters
        ----------
        threshold: float
            Threshold for pLDDT, default is 0.6

        Note
        -----
        if the last region is close to the end, merge it with the second last region
        If the second last region has already been removed, ignore it
        """
        # get regions list from low_plddt_region, keep also high plddt regions
        region_breakpoint = list(np.concatenate(self.low_plddt_region))
        # add 0 if not in the list
        if region_breakpoint[0] != 0:
            region_breakpoint = [0] + region_breakpoint
        if region_breakpoint[-1] != self.length:
            region_breakpoint.append(self.length)
        # merge regions if they are close (<30bp apart)
        # special case: if the last region is close to the end, merge it with the second last region
        # (0, 15, 60, 200, 240, 264)
        len_breakpoint = len(region_breakpoint)
        region_breakpoint_output = region_breakpoint.copy()
        for i in range(len_breakpoint - 1):
            # len(region_breakpoint) = 6; range(5) = 0, 1, 2, 3, 4
            if region_breakpoint[i + 1] - region_breakpoint[i] < 30:
                if i == len_breakpoint - 2:  # i = 4, second last region
                    if region_breakpoint[i] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i])
                else:
                    if region_breakpoint[i + 1] in region_breakpoint_output:
                        region_breakpoint_output.remove(region_breakpoint[i + 1])
        regions = []
        for i in range(len(region_breakpoint_output) - 1):
            regions.append(
                [region_breakpoint_output[i], region_breakpoint_output[i + 1]]
            )

        return regions

    @property
    def low_or_high_plddt_region_sequence(self, threshold=0.6):
        """
        Get sequence of low or high pLDDT regions

        Parameters
        ----------
        threshold: float
            Threshold for pLDDT, default is 0.6
        """
        # get regions list from low_plddt_region, keep also high plddt regions
        sequences = []
        for i, region in enumerate(self.low_or_high_plddt_region):
            s = SeqIO.SeqRecord(
                self.sequence[region[0] : region[1]],
                id=self.gene_name + "_" + str(i),
                name=self.gene_name,
                description=self.gene_name
                + "_"
                + str(i)
                + ": "
                + str(region[0])
                + "-"
                + str(region[1]),
            )
            sequences.append(s)
        return sequences

    @property
    def negative_charge_residue(self):
        """Get negative charge residues"""
        return np.isin(np.array(self.sequence), ["D", "E"])

    @property
    def positive_charge_residue(self):
        """Get positive charge residues"""
        return np.isin(np.array(self.sequence), ["R", "K", "H"])

    def get_final_score_gated_grad_3d(self, interaction_threshold=20):
        """Get final score gated by gradient and ES score"""
        f = self.grad * self.esm
        pairwise_interaction = self.pairwise_distance < interaction_threshold
        f[(self.smoothed_plddt < 0.5)] = 0
        f = get_3d_avg(f, pairwise_interaction)
        f = normalize(f)

        return f

    def get_pairwise_distance(
        self, dimer=False, af2_folder="/manitou/pmg/users/xf2217/demo_data/af2"
    ):
        """Get pairwise distance between residues"""
        if dimer:
            structure = bioparser.get_structure(
                "homodimer", f"{af2_folder}/dimer_structures/" + self.gene_name + ".pdb"
            )
            model = structure[0]
            chain = model["A"]
            residues = [r for r in model.get_residues()]
            whole_len = len(residues)
            chain_len = len(chain)
            distance = np.zeros((whole_len, whole_len))
            for i, residue1 in enumerate(residues):
                for j, residue2 in enumerate(residues):
                    # compute distance between CA atoms
                    try:
                        d = residue1["CA"] - residue2["CA"]
                        distance[i][j] = d
                        distance[j][i] = d
                    except KeyError:
                        continue
            distance = np.fmin(
                distance[0:chain_len, 0:chain_len],
                distance[0:chain_len, chain_len:whole_len],
            )

        else:
            if Path(
                f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy"
            ).exists():
                distance = np.load(
                    f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy"
                )
            else:
                # make sure structures folder exists
                if not Path(f"{af2_folder}/structures/").exists():
                    Path(f"{af2_folder}/structures/").mkdir(parents=True, exist_ok=True)

                # download pdb file from AFDB to structures/
                import urllib.request

                url = (
                    "https://alphafold.ebi.ac.uk/files/AF-"
                    + self.uniprot_id
                    + "-F1-model_v6.pdb"
                )
                urllib.request.urlretrieve(
                    url,
                    f"{af2_folder}/structures/AF-"
                    + self.uniprot_id
                    + "-F1-model_v6.pdb",
                )

                # https://alphafold.ebi.ac.uk/files/AF-Q02548-F1-model_v6.pdb
                with Path(
                    f"{af2_folder}/structures/AF-"
                    + self.uniprot_id
                    + "-F1-model_v6.pdb"
                ).open() as f:
                    structure = bioparser.get_structure("monomer", f)
                model = structure[0]
                chain = model["A"]
                residues = [r for r in model.get_residues()]
                whole_len = len(residues)
                chain_len = len(chain)
                distance = np.zeros((whole_len, whole_len))
                for i, residue1 in enumerate(residues):
                    for j, residue2 in enumerate(residues):
                        # compute distance between CA atoms
                        try:
                            d = residue1["CA"] - residue2["CA"]
                            distance[i][j] = d
                            distance[j][i] = d
                        except KeyError:
                            continue
            # make sure pairwise_interaction folder exists
            if not Path(f"{af2_folder}/pairwise_interaction/").exists():
                Path(f"{af2_folder}/pairwise_interaction/").mkdir(
                    parents=True, exist_ok=True
                )

            np.save(
                f"{af2_folder}/pairwise_interaction/" + self.uniprot_id + ".npy",
                distance,
            )
        return distance

    def plotly_plddt(
        self,
        pos_to_highlight=None,
        to_compare=None,
        filename=None,
        show_low_plddt=True,
        show_domain=True,
        domains_to_show=[
            "region of interest",
            "DNA-binding region",
            "short sequence motif",
        ],
    ):
        """Plot pLDDT using Plotly, highlight specified positions and domains

        Parameters
        ----------
        pos_to_highlight: list of int
            Positions to highlight
        to_compare: list of float
            Data to compare with
        filename: str
            Filename to save the plot
        show_low_plddt: bool
            Whether to show low pLDDT regions
        show_domain: bool
            Whether to show domains from Uniprot
        domains_to_show: list of str
            Types of domains to show. Default is ["region of interest", "DNA-binding region", "short sequence motif"].
            They are Uniprot feature types.

        Returns
        -------
        fig: plotly.graph_objects.Figure
            Plotly figure
        """
        # Initialize Plotly Figure
        fig = go.Figure()

        # Highlight low pLDDT regions
        if show_low_plddt:
            for i, region in enumerate(self.low_or_high_plddt_region):
                segment_name = f"{self.gene_name}_{str(i)} pLDDT={np.mean(self.smoothed_plddt[region[0]:region[1]]):.2f}"
                fig.add_trace(
                    go.Scatter(
                        x=[region[0], region[1], region[1], region[0]],
                        y=[0.8, 0.8, 1, 1],
                        fill="toself",
                        fillcolor="grey",
                        opacity=0.2,
                        line=dict(color="black"),
                        hoverinfo="text",
                        hovertext=[segment_name],
                        mode="lines",
                        showlegend=True,
                        legendgroup="pLDDT Segments",
                        legendgrouptitle=dict(text="pLDDT Segments"),
                        name=segment_name,
                    )
                )
        # Highlight specified positions
        if pos_to_highlight is not None:
            pos_to_highlight = np.array(pos_to_highlight) - 1
            fig.add_trace(
                go.Scatter(
                    x=pos_to_highlight,
                    y=self.smoothed_plddt[pos_to_highlight],
                    mode="markers",
                    marker=dict(color="orange", size=8),
                )
            )

        # Show domains if applicable, color by feature_type
        if show_domain:
            # Create a color mapping for each unique feature_type
            feature_types = self.domains.query(
                "feature_type.isin(@domains_to_show)"
            ).feature_description.unique()
            # use tab20 color map
            colormap = plt.get_cmap("Set3").colors
            colors = colormap[1 : len(feature_types) + 1]
            # convert to rgb string
            colors = ["rgb" + str(i) for i in colors]
            color_mapping = dict(zip(feature_types, colors))

            for i, domain in self.domains.iterrows():
                if domain.feature_type in domains_to_show:
                    # Get the color for this feature_type
                    color = color_mapping[domain.feature_description]
                    # Add domain feature as filled rectangle
                    fig.add_trace(
                        go.Scatter(
                            x=[
                                domain.feature_begin,
                                domain.feature_end,
                                domain.feature_end,
                                domain.feature_begin,
                            ],
                            y=[0.6, 0.6, 0.8, 0.8],
                            fill="toself",
                            fillcolor=color,
                            opacity=0.2,
                            line=dict(color=color),
                            hoverinfo="text",
                            hovertext=[domain.feature_description],
                            mode="lines",
                            showlegend=True,
                            legendgroup=domain.feature_type,
                            legendgrouptitle=dict(text=domain.feature_type),
                            name=domain.feature_description,
                        )
                    )

        # Plot main pLDDT line
        fig.add_trace(
            go.Scatter(
                y=self.smoothed_plddt,
                mode="lines",
                name="pLDDT",
                line=dict(color="orange"),
            )
        )

        # Plot secondary comparison line if specified
        if to_compare is not None:
            if isinstance(to_compare, dict):
                for name, (color, data) in to_compare.items():
                    fig.add_trace(
                        go.Scatter(
                            y=data,
                            mode="lines",
                            name=name,
                            line=dict(color=color),
                        )
                    )
            else:
                fig.add_trace(go.Scatter(y=to_compare, mode="lines", name="To Compare"))

        # Plot ES if available and no comparison data is specified
        elif hasattr(self, "es"):
            fig.add_trace(
                go.Scatter(y=self.es, mode="lines", name="ES", line=dict(color="blue"))
            )

        # Additional plot settings
        fig.update_layout(
            title=f"{self.gene_name} pLDDT",
            xaxis_title="Residue",
            yaxis_title="pLDDT",
            xaxis=dict(
                tickvals=np.arange(0, self.length, 100),
                ticktext=np.arange(1, self.length + 1, 100),
            ),
            font=dict(
                family="Arial",
            ),
        )

        # Save figure if filename is provided
        if filename is not None:
            fig.write_image(filename)

        return fig

    def plot_plddt_manuscript(self, to_compare=None, filename=None):
        """Plot pLDDT for manuscript

        Parameters
        ----------
        to_compare: list of float
            Data to compare with
        filename: str
            Filename to save the plot
        """
        plt.figure(figsize=(12, 3))
        plt.plot(self.plddt)
        if to_compare is not None:
            plt.plot(to_compare)
        # highlight low plddt region
        for region in self.low_plddt_region:
            plt.axvspan(region[0], region[1], ymax=1, ymin=0, color="grey", alpha=0.1)

        # highlight domain, color by feature_type
        cmap = plt.get_cmap("tab20").colors
        # map feature_type to color
        feature_type_to_color = {}
        self.domains = self.domains.query(
            '(feature_type=="domain") or (feature_type=="region of interest")'
        )
        # for i, t in enumerate(obj.domains.feature_type.unique()):
        #     feature_type_to_color[t] = cmap[i+5]
        feature_type_to_color["domain"] = cmap[5]
        feature_type_to_color["region of interest"] = cmap[6]

        y_span = 0.1  # 0.8/len(obj.domains.feature_type.unique())
        for i, domain in self.domains.iterrows():
            idx = np.where(self.domains.feature_type.unique() == domain.feature_type)[
                0
            ][0]
            plt.axvspan(
                domain.feature_begin,
                domain.feature_end,
                ymax=idx * y_span + y_span,
                ymin=idx * y_span,
                color=feature_type_to_color[domain.feature_type],
                alpha=0.2,
            )
        # add legend of domain color
        legend_elements = []
        for i in self.domains.feature_type.unique():
            legend_elements.append(Patch(facecolor=feature_type_to_color[i], label=i))
            # reverse the order of legend
        legend_elements = legend_elements[::-1]
        # add "low plddt region" to legend
        legend_elements.append(Patch(facecolor="grey", label="low pLDDT region"))
        # add number index of low or high plddt region on top of the plot
        # for i, region in enumerate(obj.low_or_high_plddt_region):
        #     plt.text(region[0], 0.9, f"{i}", fontsize=12)
        # legend outside the plot
        plt.legend(
            handles=legend_elements,
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
        )
        plt.ylabel(f"{self.gene_name} pLDDT")
        plt.xlabel("Residue")
        # set xlim to the length of protein
        plt.xlim(0, len(self.plddt))
        # plt.ylabel("pLDDT")
        # plt.tight_layout()
        if filename is not None:
            plt.savefig(filename, dpi=300, bbox_inches="tight")
        plt.show()

    def plot_plddt(
        self,
        pos_to_highlight=None,
        to_compare=None,
        filename=None,
        show_low_plddt=True,
        show_domain=True,
        domains_to_show=["region of interest", "DNA-binding region", "splice variant"],
    ):
        """Plot pLDDT using matplotlib

        Parameters
        ----------
        pos_to_highlight: list of int
            Positions to highlight
        to_compare: list of float
            Data to compare with
        filename: str
            Filename to save the plot
        show_low_plddt: bool
            Whether to show low pLDDT regions
        show_domain: bool
            Whether to show domains from Uniprot
        domains_to_show: list of str
            Types of domains to show. Default is ["region of interest", "DNA-binding region", "short sequence motif"].
            They are Uniprot feature types.
        """
        fig, ax = plt.subplots(figsize=(20, 5))
        plt.plot(self.smoothed_plddt, label="pLDDT", color="orange")
        if to_compare is not None:
            plt.plot(to_compare)
        elif hasattr(self, "es"):
            plt.plot(self.es, label="ES", color="blue")
        if show_low_plddt:
            # highlight low plddt region
            for region in self.low_plddt_region:
                plt.axvspan(
                    region[0], region[1], ymax=1, ymin=0.8, color="red", alpha=0.2
                )

        if pos_to_highlight is not None:
            pos_to_highlight = np.array(pos_to_highlight) - 1
            if hasattr(self, "es"):
                plt.scatter(
                    pos_to_highlight, self.es[pos_to_highlight], color="blue", s=50
                )
            plt.scatter(
                pos_to_highlight,
                self.smoothed_plddt[pos_to_highlight],
                color="orange",
                s=50,
            )

        if show_domain:
            if domains_to_show is None:
                domains_to_show = self.domains.feature_type.unique()
            else:
                domains_to_show = np.array(domains_to_show)

            domains = self.domains.query("feature_type.isin(@domains_to_show)")
            # highlight domain, color by feature_type
            cmap = plt.get_cmap("Set3").colors
            # map feature_type to color
            feature_type_to_color = {}
            for i, t in enumerate(domains_to_show):
                feature_type_to_color[t] = cmap[i]

            y_span = 0.8 / len(domains_to_show)

            for i, domain in domains.query(
                "feature_type.isin(@domains_to_show)"
            ).iterrows():
                idx = np.where(domains_to_show == domain.feature_type)[0][0]
                plt.axvspan(
                    domain.feature_begin - 1,
                    domain.feature_end - 1,
                    ymax=idx * y_span + y_span,
                    ymin=idx * y_span,
                    color=feature_type_to_color[domain.feature_type],
                    alpha=0.2,
                )
            # add legend of domain color
            legend_elements = []
            for i in domains_to_show:
                legend_elements.append(
                    Patch(facecolor=feature_type_to_color[i], label=i)
                )
                # reverse the order of legend
            legend_elements = legend_elements[::-1]
            # legend outside the plot, append line plot legend in the end
            legend_elements.append(plt.Line2D([0], [0], color="blue", label="ES"))
            legend_elements.append(plt.Line2D([0], [0], color="orange", label="pLDDT"))
            plt.legend(
                handles=legend_elements,
                bbox_to_anchor=(1.05, 1, 0.1, 0.1),
                loc="upper left",
                borderaxespad=0.0,
                fontsize=16,
            )
        # add number index of low or high plddt region on top of the plot
        for i, region in enumerate(self.low_or_high_plddt_region):
            plt.text(region[0], 0.9, f"{i}", fontsize=12)

        plt.title(f"{self.gene_name} pLDDT")
        plt.xlabel("Residue")
        plt.ylabel("pLDDT")
        # set x ticks to start from 1
        plt.xticks(np.arange(0, self.length, 100), np.arange(1, self.length + 1, 100))
        plt.tight_layout()
        if filename is not None:
            plt.savefig(filename, dpi=300, bbox_inches="tight")
        return fig, ax
