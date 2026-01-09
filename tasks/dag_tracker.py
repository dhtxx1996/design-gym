"""
DAG Tracker for Agent Execution Traces

Tracks artifacts (data/files) and functions (tool calls) as nodes in a DAG,
with edges representing data flow between them.

Node Types:
- artifact: Data or file (input or output)
- function: Tool/function call

Edge Types:
- input: artifact -> function (function consumes artifact)
- output: function -> artifact (function produces artifact)
"""

import json
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Literal, Optional
from datetime import datetime


@dataclass
class Node:
    """A node in the execution DAG."""
    id: str
    type: Literal["artifact", "function"]
    label: str
    metadata: dict = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        return asdict(self)


@dataclass  
class Edge:
    """An edge connecting nodes in the DAG."""
    source: str  # node id
    target: str  # node id
    label: str = ""
    metadata: dict = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        return asdict(self)


class DAGTracker:
    """
    Tracks execution as a DAG with artifacts and functions as nodes.
    
    The LLM can use add_node() and add_edge() to explicitly build the graph
    as it executes tools. This provides:
    - Full provenance tracking
    - Visualization of data flow
    - Reproducibility information
    """
    
    def __init__(self):
        self.nodes: dict[str, Node] = {}
        self.edges: list[Edge] = []
        self._node_counter = 0
        self._edge_counter = 0
        
    def _generate_id(self, prefix: str) -> str:
        """Generate a unique node ID."""
        self._node_counter += 1
        return f"{prefix}_{self._node_counter}"
    
    def add_artifact(self, label: str, artifact_id: str = None, **metadata) -> str:
        """
        Add an artifact node (data/file).
        
        Args:
            label: Human-readable name (e.g., "scaffold.pdb", "sequence_data")
            artifact_id: Optional explicit ID, auto-generated if not provided
            **metadata: Additional info (path, size, type, etc.)
            
        Returns:
            The node ID
        """
        node_id = artifact_id or self._generate_id("artifact")
        self.nodes[node_id] = Node(
            id=node_id,
            type="artifact", 
            label=label,
            metadata=metadata
        )
        return node_id
    
    def add_function(self, label: str, function_id: str = None, **metadata) -> str:
        """
        Add a function node (tool call).
        
        Args:
            label: Function/tool name (e.g., "run_python", "proteinmpnn")
            function_id: Optional explicit ID, auto-generated if not provided
            **metadata: Additional info (arguments, result, etc.)
            
        Returns:
            The node ID
        """
        node_id = function_id or self._generate_id("func")
        self.nodes[node_id] = Node(
            id=node_id,
            type="function",
            label=label,
            metadata=metadata
        )
        return node_id
    
    def add_edge(self, source: str, target: str, label: str = "", **metadata) -> None:
        """
        Add an edge from source to target node.
        
        Args:
            source: Source node ID
            target: Target node ID  
            label: Optional edge label
            **metadata: Additional edge info
        """
        if source not in self.nodes:
            raise ValueError(f"Source node '{source}' not found")
        if target not in self.nodes:
            raise ValueError(f"Target node '{target}' not found")
            
        self.edges.append(Edge(
            source=source,
            target=target,
            label=label,
            metadata=metadata
        ))
    
    def connect(self, inputs: list[str], function_id: str, outputs: list[str]) -> None:
        """
        Convenience method to connect inputs -> function -> outputs.
        
        Args:
            inputs: List of input artifact node IDs
            function_id: The function node ID
            outputs: List of output artifact node IDs
        """
        for inp in inputs:
            self.add_edge(inp, function_id)
        for out in outputs:
            self.add_edge(function_id, out)
    
    def get_node(self, node_id: str) -> Optional[Node]:
        """Get a node by ID."""
        return self.nodes.get(node_id)
    
    def list_nodes(self, node_type: str = None) -> list[Node]:
        """List all nodes, optionally filtered by type."""
        nodes = list(self.nodes.values())
        if node_type:
            nodes = [n for n in nodes if n.type == node_type]
        return nodes
    
    def list_edges(self) -> list[Edge]:
        """List all edges."""
        return self.edges.copy()
    
    def to_dict(self) -> dict:
        """Serialize the DAG to a dictionary."""
        return {
            "nodes": [n.to_dict() for n in self.nodes.values()],
            "edges": [e.to_dict() for e in self.edges],
            "metadata": {
                "node_count": len(self.nodes),
                "edge_count": len(self.edges),
                "artifact_count": len([n for n in self.nodes.values() if n.type == "artifact"]),
                "function_count": len([n for n in self.nodes.values() if n.type == "function"]),
            }
        }
    
    def save(self, path: str | Path) -> None:
        """Save the DAG to a JSON file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, path: str | Path) -> "DAGTracker":
        """Load a DAG from a JSON file."""
        with open(path) as f:
            data = json.load(f)
        
        dag = cls()
        for node_data in data.get("nodes", []):
            node = Node(**node_data)
            dag.nodes[node.id] = node
            
        for edge_data in data.get("edges", []):
            dag.edges.append(Edge(**edge_data))
            
        return dag
    
    def _find_connected_components(self) -> list[set[str]]:
        """Find connected components in the DAG (treating as undirected graph)."""
        # Build undirected adjacency
        adj = {nid: set() for nid in self.nodes}
        for edge in self.edges:
            adj[edge.source].add(edge.target)
            adj[edge.target].add(edge.source)
        
        visited = set()
        components = []
        
        def bfs(start: str) -> set[str]:
            component = set()
            queue = [start]
            while queue:
                node = queue.pop(0)
                if node in visited:
                    continue
                visited.add(node)
                component.add(node)
                for neighbor in adj[node]:
                    if neighbor not in visited:
                        queue.append(neighbor)
            return component
        
        for node_id in self.nodes:
            if node_id not in visited:
                component = bfs(node_id)
                if component:
                    components.append(component)
        
        return components
    
    def visualize(self, output_path: str | Path = None, title: str = "Execution DAG", 
                  figsize: tuple = (14, 10), dpi: int = 150,
                  min_cluster_size: int = 3, show_all: bool = False) -> None:
        """
        Visualize the DAG with a pleasant hierarchical layout.
        
        Uses Sugiyama-style layered layout for left-to-right flow.
        
        Args:
            output_path: Path to save the figure (optional)
            title: Figure title
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
            min_cluster_size: Minimum number of nodes in a cluster to display (default: 3)
            show_all: If True, show all nodes regardless of cluster size (default: False)
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyBboxPatch
        import numpy as np
        
        if not self.nodes:
            print("No nodes to visualize")
            return
        
        # Filter nodes by cluster size if not showing all
        if show_all:
            nodes_to_show = set(self.nodes.keys())
        else:
            components = self._find_connected_components()
            nodes_to_show = set()
            for component in components:
                if len(component) >= min_cluster_size:
                    nodes_to_show.update(component)
        
        if not nodes_to_show:
            print(f"No clusters with at least {min_cluster_size} nodes to visualize")
            return
        
        # Filter edges to only include those between shown nodes
        edges_to_show = [e for e in self.edges if e.source in nodes_to_show and e.target in nodes_to_show]
        
        # Build adjacency for layout computation (using filtered nodes)
        adj = {nid: [] for nid in nodes_to_show}
        in_degree = {nid: 0 for nid in nodes_to_show}
        for edge in edges_to_show:
            adj[edge.source].append(edge.target)
            in_degree[edge.target] += 1
        
        # Compute layers using topological sort (Kahn's algorithm)
        layers = self._compute_layers_filtered(adj, in_degree, nodes_to_show)
        
        # Compute positions with layer-based layout
        positions = self._compute_positions_filtered(layers, edges_to_show)
        
        # Create figure - compute optimal size based on graph dimensions
        max_layer_size = max(len(layer) for layer in layers) if layers else 1
        num_layers = len(layers)
        
        # Dynamic figure size: wider for more layers, taller for bigger layers
        fig_width = max(figsize[0], num_layers * 4.0)  # More width per layer
        fig_height = max(figsize[1], max_layer_size * 0.8)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        # Don't force equal aspect - let it stretch naturally
        
        # Color scheme (matching example image style)
        colors = {
            "artifact": {
                "file": "#4A90D9",      # Blue - files
                "data": "#50C878",      # Green - data  
                "result": "#FFD700",    # Yellow/Gold - results
                "input": "#DDA0DD",     # Pink/Plum - inputs
                "default": "#87CEEB",   # Light blue
            },
            "function": {
                "run_python": "#FF8C00",    # Orange
                "tamarind": "#32CD32",      # Lime green
                "read": "#6495ED",          # Cornflower blue
                "write": "#FFB347",         # Pastel orange
                "default": "#FF7F50",       # Coral
            }
        }
        
        # Draw edges first (behind nodes)
        self._draw_edges_filtered(ax, positions, edges_to_show)
        
        # Draw nodes
        self._draw_nodes_filtered(ax, positions, colors, nodes_to_show)
        
        # Set limits with padding
        all_x = [p[0] for p in positions.values()]
        all_y = [p[1] for p in positions.values()]
        
        x_margin = 1.5
        y_margin = 1.0
        ax.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)
        ax.set_ylim(min(all_y) - y_margin, max(all_y) + y_margin)
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        
        if output_path:
            path = Path(output_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(path, dpi=dpi, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            print(f"DAG visualization saved to: {path}")
        
        plt.close()
    
    def _compute_layers(self, adj: dict, in_degree: dict) -> list[list[str]]:
        """Compute node layers using Kahn's algorithm with balanced distribution."""
        return self._compute_layers_filtered(adj, in_degree, set(self.nodes.keys()))
    
    def _compute_layers_filtered(self, adj: dict, in_degree: dict, nodes_to_show: set[str]) -> list[list[str]]:
        """Compute node layers using Kahn's algorithm with balanced distribution."""
        layers = []
        remaining = set(nodes_to_show)
        current_in_degree = in_degree.copy()
        
        while remaining:
            # Find all nodes with in-degree 0 among remaining
            layer = [n for n in remaining if current_in_degree[n] == 0]
            
            if not layer:
                # Handle cycles - just take remaining nodes
                layer = list(remaining)
            
            layers.append(layer)
            
            # Remove these nodes and update in-degrees
            for node in layer:
                remaining.discard(node)
                for neighbor in adj.get(node, []):
                    current_in_degree[neighbor] -= 1
        
        return layers
    
    def _compute_positions(self, layers: list[list[str]]) -> dict[str, tuple]:
        """Compute node positions from layers with multi-pass crossing reduction."""
        return self._compute_positions_filtered(layers, self.edges)
    
    def _compute_positions_filtered(self, layers: list[list[str]], edges: list[Edge]) -> dict[str, tuple]:
        """Compute node positions from layers with multi-pass crossing reduction."""
        x_spacing = 4.5  # Horizontal spacing between layers
        
        # Dynamic y_spacing based on max layer size for more compact layout
        max_layer_size = max(len(layer) for layer in layers) if layers else 1
        y_spacing = max(1.0, min(1.4, 20 / max_layer_size))  # Tighter for large layers
        
        # Collect all nodes in layers
        all_layer_nodes = set()
        for layer in layers:
            all_layer_nodes.update(layer)
        
        # Build adjacency structures
        predecessors = {nid: [] for nid in all_layer_nodes}
        successors = {nid: [] for nid in all_layer_nodes}
        for edge in edges:
            if edge.source in all_layer_nodes and edge.target in all_layer_nodes:
                predecessors[edge.target].append(edge.source)
                successors[edge.source].append(edge.target)
        
        # Create mutable layer lists and initial ordering
        layers = [list(layer) for layer in layers]
        
        # Assign initial y-positions (needed for barycenter)
        node_y = {}
        for layer in layers:
            for i, node in enumerate(layer):
                node_y[node] = i
        
        def median(values):
            """Get median of values list."""
            if not values:
                return 0
            s = sorted(values)
            n = len(s)
            if n % 2 == 1:
                return s[n // 2]
            return (s[n // 2 - 1] + s[n // 2]) / 2
        
        # Multi-pass barycenter/median heuristic
        num_passes = 4
        for pass_num in range(num_passes):
            # Forward pass (left to right)
            for layer_idx in range(1, len(layers)):
                layer = layers[layer_idx]
                
                def forward_score(node):
                    preds = predecessors[node]
                    if not preds:
                        return node_y.get(node, 0)
                    return median([node_y[p] for p in preds if p in node_y])
                
                layers[layer_idx] = sorted(layer, key=forward_score)
                # Update y positions
                for i, node in enumerate(layers[layer_idx]):
                    node_y[node] = i
            
            # Backward pass (right to left)
            for layer_idx in range(len(layers) - 2, -1, -1):
                layer = layers[layer_idx]
                
                def backward_score(node):
                    succs = successors[node]
                    if not succs:
                        return node_y.get(node, 0)
                    return median([node_y[s] for s in succs if s in node_y])
                
                layers[layer_idx] = sorted(layer, key=backward_score)
                # Update y positions
                for i, node in enumerate(layers[layer_idx]):
                    node_y[node] = i
        
        # Compute final positions
        positions = {}
        for layer_idx, layer in enumerate(layers):
            x = layer_idx * x_spacing
            
            # Center the layer vertically
            layer_height = (len(layer) - 1) * y_spacing
            start_y = -layer_height / 2
            
            for i, node_id in enumerate(layer):
                y = start_y + i * y_spacing
                positions[node_id] = (x, y)
        
        return positions
    
    def _draw_edges(self, ax, positions: dict) -> None:
        """Draw curved bezier edges between nodes with smart routing."""
        self._draw_edges_filtered(ax, positions, self.edges)
    
    def _draw_edges_filtered(self, ax, positions: dict, edges: list[Edge]) -> None:
        """Draw curved bezier edges between nodes with smart routing."""
        from matplotlib.patches import FancyArrowPatch
        from matplotlib.path import Path as MPath
        import matplotlib.patches as mpatches
        
        for edge in edges:
            if edge.source not in positions or edge.target not in positions:
                continue
                
            x1, y1 = positions[edge.source]
            x2, y2 = positions[edge.target]
            
            # Compute node widths based on label lengths
            src_label = self.nodes[edge.source].label
            tgt_label = self.nodes[edge.target].label
            src_width = max(1.5, len(src_label) * 0.18)
            tgt_width = max(1.5, len(tgt_label) * 0.18)
            
            x1 += src_width / 2
            x2 -= tgt_width / 2
            
            # Compute bezier control points based on distance
            dx = x2 - x1
            dy = y2 - y1
            
            # Adjust curve based on vertical distance
            # More vertical distance = more horizontal control offset for smoother S-curves
            vertical_factor = min(abs(dy) / 3.0, 1.0)  # Normalize
            base_offset = dx * 0.35
            extra_offset = dx * 0.15 * vertical_factor
            ctrl_offset = base_offset + extra_offset
            
            # For edges going backwards (shouldn't happen in DAG, but handle it)
            if dx <= 0:
                ctrl_offset = max(1.0, abs(dy) * 0.5)
            
            # Path with bezier curve - smooth S-curve
            verts = [
                (x1, y1),
                (x1 + ctrl_offset, y1),
                (x2 - ctrl_offset, y2),
                (x2, y2),
            ]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
            path = MPath(verts, codes)
            
            # Implicit edges (from reasoning) get dashed lines
            linestyle = '--' if edge.label == "implicit" else '-'
            edgecolor = '#AAAAAA' if edge.label == "implicit" else '#888888'

            patch = mpatches.PathPatch(
                path,
                facecolor='none',
                edgecolor=edgecolor,
                linewidth=1.2,
                alpha=0.6,
                linestyle=linestyle,
            )
            ax.add_patch(patch)
            
            # No arrowhead - clean line only
    
    def _draw_nodes(self, ax, positions: dict, colors: dict) -> None:
        """Draw nodes as rounded rectangles with labels."""
        self._draw_nodes_filtered(ax, positions, colors, set(self.nodes.keys()))
    
    def _draw_nodes_filtered(self, ax, positions: dict, colors: dict, nodes_to_show: set[str]) -> None:
        """Draw nodes as rounded rectangles with labels."""
        from matplotlib.patches import FancyBboxPatch
        
        for node_id, (x, y) in positions.items():
            if node_id not in nodes_to_show:
                continue
            node = self.nodes[node_id]
            
            # Determine color based on node type and metadata
            if node.type == "artifact":
                color = self._get_artifact_color(node, colors["artifact"])
            else:
                color = self._get_function_color(node, colors["function"])
            
            # Node dimensions - width based on label length
            label = node.label
            width = max(1.5, len(label) * 0.18)
            height = 0.55  # 1.1x taller than original 0.5
            
            # Draw rounded rectangle
            bbox = FancyBboxPatch(
                (x - width/2, y - height/2),
                width, height,
                boxstyle="round,pad=0.03,rounding_size=0.1",
                facecolor=color,
                edgecolor='#333333',
                linewidth=1.5,
            )
            ax.add_patch(bbox)
            
            # Add full label (no truncation) - monospace font
            ax.text(x, y, label, ha='center', va='center',
                   fontsize=7.5, fontweight='medium', fontfamily='monospace', color='#222222')
    
    def _get_artifact_color(self, node: Node, palette: dict) -> str:
        """Get color for artifact node based on its properties."""
        label_lower = node.label.lower()
        
        if any(ext in label_lower for ext in ['.pdb', '.fa', '.fasta', '.csv', '.json']):
            return palette["file"]
        elif 'input' in label_lower or 'source' in label_lower:
            return palette["input"]
        elif 'result' in label_lower or 'output' in label_lower or 'answer' in label_lower:
            return palette["result"]
        elif 'data' in label_lower:
            return palette["data"]
        else:
            return palette["default"]
    
    def _get_function_color(self, node: Node, palette: dict) -> str:
        """Get color for function node based on its properties."""
        label_lower = node.label.lower()
        
        if 'python' in label_lower:
            return palette["run_python"]
        elif 'tamarind' in label_lower or any(t in label_lower for t in ['proteinmpnn', 'esmfold', 'alphafold']):
            return palette["tamarind"]
        elif 'read' in label_lower:
            return palette["read"]
        elif 'write' in label_lower:
            return palette["write"]
        else:
            return palette["default"]


def create_toy_dag() -> DAGTracker:
    """Create a toy DAG for testing visualization layout."""
    dag = DAGTracker()
    
    # Layer 0: Input data
    input_pdb = dag.add_artifact("scaffold.pdb", artifact_id="input_pdb")
    input_seq = dag.add_artifact("input_seq.fa", artifact_id="input_seq")
    
    # Layer 1: Initial processing
    func1 = dag.add_function("run_python", function_id="parse_pdb")
    func2 = dag.add_function("read_file", function_id="read_seq")
    
    dag.add_edge(input_pdb, func1)
    dag.add_edge(input_seq, func2)
    
    # Intermediate artifacts
    parsed_data = dag.add_artifact("parsed_data", artifact_id="parsed")
    seq_data = dag.add_artifact("seq_data", artifact_id="seq")
    
    dag.add_edge(func1, parsed_data)
    dag.add_edge(func2, seq_data)
    
    # Layer 2: ML tool calls
    func3 = dag.add_function("proteinmpnn", function_id="design")
    func4 = dag.add_function("esmfold", function_id="fold")
    
    dag.add_edge(parsed_data, func3)
    dag.add_edge(seq_data, func3)
    dag.add_edge(seq_data, func4)
    
    # More artifacts
    designs = dag.add_artifact("designs.fa", artifact_id="designs")
    structure = dag.add_artifact("folded.pdb", artifact_id="structure")
    
    dag.add_edge(func3, designs)
    dag.add_edge(func4, structure)
    
    # Layer 3: Analysis
    func5 = dag.add_function("run_python", function_id="analyze")
    
    dag.add_edge(designs, func5)
    dag.add_edge(structure, func5)
    
    # Final outputs
    results = dag.add_artifact("results.json", artifact_id="results")
    answer = dag.add_artifact("answer.txt", artifact_id="answer")
    
    dag.add_edge(func5, results)
    dag.add_edge(func5, answer)
    
    # Add some additional branches for complexity
    func6 = dag.add_function("write_file", function_id="save_log")
    log = dag.add_artifact("log.txt", artifact_id="log")
    dag.add_edge(parsed_data, func6)
    dag.add_edge(func6, log)
    
    return dag


if __name__ == "__main__":
    # Test with toy DAG
    print("Creating toy DAG for visualization test...")
    dag = create_toy_dag()
    
    print(f"Nodes: {len(dag.nodes)}")
    print(f"Edges: {len(dag.edges)}")
    
    # Save and visualize
    dag.save("toy_dag.json")
    print("Saved to toy_dag.json")
    
    dag.visualize(output_path="toy_dag.png", title="Toy Execution DAG")

