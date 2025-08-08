from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.community as nxcom
import numpy as np
import plotly.graph_objects as go


def get_top_edge_weight(G, n=1000):
    """Get the top n edge weights from the graph.

    Parameters
    ----------
    G: nx.Graph
        The graph to get the top edge weights from.
    n: int, optional
        The number of edge weights to get, defaults to 1000.

    Returns
    -------
    np.array
        The top n edge weights.
    """
    edge_weight_list = []
    for u, v, d in G.edges(data=True):
        edge_weight_list.append(d["weight"])
    w = np.array(edge_weight_list)
    w = np.absolute(w)
    w = w[np.argsort(w)[-n:][0]]
    return w


def set_node_community(G, communities):
    """Add community to node attributes.

    Parameters
    ----------
    G: nx.Graph
        The graph to add the community to.
    communities: list
        The communities to add to the graph.
    """
    for c, v_c in enumerate(communities):
        for v in v_c:
            # Add 1 to save 0 for external edges
            G.nodes[v]["community"] = c + 1


def set_edge_community(G):
    """Find internal edges and add their community to their attributes.

    Parameters
    ----------
    G: nx.Graph
        The graph to add the community to.
    """
    for (
        v,
        w,
    ) in G.edges:
        if G.nodes[v]["community"] == G.nodes[w]["community"]:
            # Internal edge, mark with community
            G.edges[v, w]["community"] = G.nodes[v]["community"]
        else:
            # External edge, mark as 0
            G.edges[v, w]["community"] = 0


def get_color(i, r_off=1, g_off=1, b_off=1):
    """Assign a color to a vertex.

    Parameters
    ----------
    i: int
        The index of the vertex.
    r_off: int, optional
        The offset for the red channel, defaults to 1.
    g_off: int, optional
        The offset for the green channel, defaults to 1.
    b_off: int, optional
        The offset for the blue channel, defaults to 1.
    """
    n = 16
    low, high = 0.5, 0.9
    span = high - low
    r = low + span * (((i + r_off) * 3) % n) / (n - 1)
    g = low + span * (((i + g_off) * 5) % n) / (n - 1)
    b = low + span * (((i + b_off) * 7) % n) / (n - 1)
    return (r, g, b)


def plot_comm(G, figsize=(10, 10), title="Network structure", savefig=False):
    """Plot the community structure of the graph.

    Parameters
    ----------
    G: nx.Graph
        The graph to plot.
    figsize: tuple, optional
        The size of the figure, defaults to (10, 10).
    title: str, optional
        The title of the plot, defaults to "Network structure".
    savefig: str, optional
        The path to save the plot, defaults to False.
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    # Set community color for internal edges
    external = [(v, w) for v, w in G.edges if G.edges[v, w]["community"] == 0]
    internal = [(v, w) for v, w in G.edges if G.edges[v, w]["community"] > 0]
    internal_color = [
        "#fddbc7" if (G.edges[e]["weight"] > 0) else "#d1e5f0" for e in internal
    ]
    external_color = [
        "#ef8a62" if (G.edges[e]["weight"] > 0) else "#67a9cf" for e in external
    ]
    node_color = [get_color(G.nodes[v]["community"]) for v in G.nodes]

    weights = [np.absolute(G[u][v]["weight"]) for u, v in G.edges]
    # external edges
    pos = nx.spring_layout(G)
    nx.draw(
        G,
        pos=pos,
        node_size=0,
        edgelist=external,
        width=weights,
        ax=ax,
        edge_color=external_color,
        node_color=node_color,
        with_labels=True,
    )
    # internal edges
    nx.draw(
        G,
        pos=pos,
        node_size=100,
        edgelist=internal,
        ax=ax,
        width=weights,
        edge_color=internal_color,
        node_color=node_color,
        with_labels=True,
    )
    ax.set_title(title)
    plt.tight_layout()
    # add padding in the margins
    plt.margins(0.3, 0.3)
    if savefig:
        save_path = Path(savefig) / f"{title}.pdf"
        plt.savefig(save_path, dpi=300)
    plt.show()
    plt.close()


def plotly_networkx_digraph(
    G: nx.DiGraph, hoverinfo_dict=None, node_weights_dict=None
) -> go.Figure:
    """Plot the graph using plotly.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to plot.
    hoverinfo_dict: dict, optional
        The hover information for the nodes, defaults to None.
    node_weights_dict: dict, optional
        The node weights for the nodes, defaults to None.
    """
    pos = nx.spring_layout(G)

    edge_count = len(G.edges())
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("RdBu_r")
    edge_colors = [cmap(i) for i in np.linspace(0, 1, 256)]
    edge_color_dict = {}
    for i, (u, v, weight) in enumerate(G.edges(data="weight")):
        color_idx = int((weight + 1) * 127.5)  # Scale weight to [0, 255]
        # edge_colors[color_idx] to rgb string
        rgb_v = f"rgb({edge_colors[color_idx][0]}, {edge_colors[color_idx][1]}, {edge_colors[color_idx][2]})"
        edge_color_dict[(u, v)] = rgb_v
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]

    if node_weights_dict:
        node_sizes = np.array([node_weights_dict[node] * 25 for node in G.nodes()])
        # replace nan with 1
        node_sizes[np.isnan(node_sizes)] = 1
    else:
        node_sizes = [25] * len(G.nodes())

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers",
        hoverinfo="text",
        marker=dict(
            showscale=True,
            colorscale="YlGnBu_r",
            color=node_sizes / 25,
            size=node_sizes,
            colorbar=dict(
                thickness=15,
                title=dict(text="Mean TF Expression"),
                xanchor="left",
            ),
        ),
    )

    label_trace = go.Scatter(
        x=node_x,
        y=np.array(node_y),
        mode="text",
        hoverinfo="none",
    )

    # Create node labels
    node_text = [f"{node}" for node in G.nodes()]
    label_trace.text = node_text

    if hoverinfo_dict is not None:
        node_text = [
            f"{node} TF expressions:<br />{hoverinfo_dict[node]}" for node in G.nodes()
        ]
    node_trace.text = node_text
    traces = []

    if edge_count < 20:
        for u, v, weight in G.edges(data="weight"):
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_trace = go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                line=dict(width=np.abs(weight) * 10, color=edge_color_dict[(u, v)]),
                mode="lines",
                hoverinfo="text",
                text=f"{u} -> {v} weight: {str(weight)}",
            )
            traces.append(edge_trace)
    else:
        edge_x, edge_y = [], []
        for u, v in G.edges():
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=0.5, color="#888"),
            hoverinfo="none",
            mode="lines",
        )
        traces.append(edge_trace)

    traces.append(node_trace)
    traces.append(label_trace)
    fig = go.Figure(
        data=traces,
        layout=go.Layout(
            showlegend=False,
            hovermode="closest",
            margin=dict(b=0, l=0, r=0, t=0),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="white",
        ),
    )

    return fig


def get_parents_subnet(G, n):
    """Get the subnet of the parents of a node.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to get the subnet from.
    n: str
        The node to get the subnet of.
    """
    return G.subgraph(list(G.predecessors(n)) + [n])


def get_children_subnet(G, n):
    """Get the subnet of the children of a node.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to get the subnet from.
    n: str
        The node to get the subnet of.
    """
    return G.subgraph(list(G.successors(n)) + [n])


def get_neighbors_subnet(G, n):
    """Get the subnet of the neighbors of a node.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to get the subnet from.
    n: str
        The node to get the subnet of.
    """
    return G.subgraph(list(nx.all_neighbors(G, n)) + [n])


def get_subnet(G, n, type="neighbors"):
    """Get the subnet of a node.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to get the subnet from.
    n: str
        The node to get the subnet of.
    type: str, optional
        The type of subnet to get, defaults to "neighbors".
    """
    if type == "parents":
        return get_parents_subnet(G, n)
    elif type == "children":
        return get_children_subnet(G, n)
    elif type == "neighbors":
        return get_neighbors_subnet(G, n)
    else:
        raise ValueError("type must be one of parents, children, or neighbors")


def preprocess_net(G, threshold=0.0, remove_nodes=True, detect_communities=True):
    """Preprocess the network.

    Parameters
    ----------
    G: nx.DiGraph
        The graph to preprocess.
    threshold: float, optional
        The threshold for the edge weights, defaults to 0.0.
    remove_nodes: bool, optional
        Whether to remove the nodes, defaults to True.
    detect_communities: bool, optional
        Whether to detect the communities, defaults to True.
    """
    G.remove_edges_from(
        [
            (n1, n2)
            for n1, n2, w in G.edges(data="weight")
            if (w < threshold and w > -threshold)
        ]
    )
    G.remove_edges_from(nx.selfloop_edges(G))
    if remove_nodes:
        G.remove_nodes_from(list(nx.isolates(G)))
    if detect_communities:
        communities = sorted(
            nxcom.greedy_modularity_communities(G), key=len, reverse=True
        )
        # Set node and edge communities
        set_node_community(G, communities)
        set_edge_community(G)
    return G
