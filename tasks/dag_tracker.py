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
    
    def visualize(self, output_path: str | Path = None, title: str = "Execution DAG", 
                  figsize: tuple = (14, 10), dpi: int = 150) -> None:
        """
        Visualize the DAG with a pleasant hierarchical layout.
        
        Uses Sugiyama-style layered layout for left-to-right flow.
        
        Args:
            output_path: Path to save the figure (optional)
            title: Figure title
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyBboxPatch
        import numpy as np
        
        if not self.nodes:
            print("No nodes to visualize")
            return
        
        # Build adjacency for layout computation
        adj = {nid: [] for nid in self.nodes}
        in_degree = {nid: 0 for nid in self.nodes}
        for edge in self.edges:
            adj[edge.source].append(edge.target)
            in_degree[edge.target] += 1
        
        # Compute layers using topological sort (Kahn's algorithm)
        layers = self._compute_layers(adj, in_degree)
        
        # Compute positions with layer-based layout
        positions = self._compute_positions(layers)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_aspect('equal')
        
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
        self._draw_edges(ax, positions)
        
        # Draw nodes
        self._draw_nodes(ax, positions, colors)
        
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
        """Compute node layers using modified Kahn's algorithm."""
        layers = []
        remaining = set(self.nodes.keys())
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
        """Compute node positions from layers with optimized ordering."""
        positions = {}
        
        x_spacing = 3.5  # Horizontal spacing between layers (wider for full labels)
        y_spacing = 1.2  # Vertical spacing between nodes in same layer
        
        # Build reverse adjacency for ordering optimization
        predecessors = {nid: [] for nid in self.nodes}
        for edge in self.edges:
            predecessors[edge.target].append(edge.source)
        
        for layer_idx, layer in enumerate(layers):
            x = layer_idx * x_spacing
            
            # Order nodes in layer to minimize edge crossings (barycenter heuristic)
            if layer_idx > 0:
                def barycenter(node):
                    preds = predecessors[node]
                    if not preds:
                        return 0
                    return sum(positions[p][1] for p in preds if p in positions) / len(preds)
                layer = sorted(layer, key=barycenter)
            
            # Center the layer vertically
            layer_height = (len(layer) - 1) * y_spacing
            start_y = -layer_height / 2
            
            for i, node_id in enumerate(layer):
                y = start_y + i * y_spacing
                positions[node_id] = (x, y)
        
        return positions
    
    def _draw_edges(self, ax, positions: dict) -> None:
        """Draw curved bezier edges between nodes."""
        from matplotlib.patches import FancyArrowPatch
        from matplotlib.path import Path as MPath
        import matplotlib.patches as mpatches
        
        for edge in self.edges:
            if edge.source not in positions or edge.target not in positions:
                continue
                
            x1, y1 = positions[edge.source]
            x2, y2 = positions[edge.target]
            
            # Compute node widths based on label lengths
            src_label = self.nodes[edge.source].label
            tgt_label = self.nodes[edge.target].label
            src_width = max(1.5, len(src_label) * 0.12)
            tgt_width = max(1.5, len(tgt_label) * 0.12)
            
            x1 += src_width / 2
            x2 -= tgt_width / 2
            
            # Bezier control points for smooth curves
            dx = x2 - x1
            ctrl_offset = dx * 0.4
            
            # Path with bezier curve
            verts = [
                (x1, y1),
                (x1 + ctrl_offset, y1),
                (x2 - ctrl_offset, y2),
                (x2, y2),
            ]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
            path = MPath(verts, codes)
            
            patch = mpatches.PathPatch(
                path,
                facecolor='none',
                edgecolor='#888888',
                linewidth=1.5,
                alpha=0.7,
            )
            ax.add_patch(patch)
            
            # Add arrowhead
            arrow = FancyArrowPatch(
                (x2 - 0.15, y2), (x2, y2),
                arrowstyle='-|>',
                mutation_scale=12,
                color='#888888',
                linewidth=1.5,
            )
            ax.add_patch(arrow)
    
    def _draw_nodes(self, ax, positions: dict, colors: dict) -> None:
        """Draw nodes as rounded rectangles with labels."""
        from matplotlib.patches import FancyBboxPatch
        
        for node_id, (x, y) in positions.items():
            node = self.nodes[node_id]
            
            # Determine color based on node type and metadata
            if node.type == "artifact":
                color = self._get_artifact_color(node, colors["artifact"])
            else:
                color = self._get_function_color(node, colors["function"])
            
            # Node dimensions - width based on label length
            label = node.label
            width = max(1.5, len(label) * 0.12)
            height = 0.5
            
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
            
            # Add full label (no truncation)
            ax.text(x, y, label, ha='center', va='center',
                   fontsize=8, fontweight='bold', color='#222222')
    
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

